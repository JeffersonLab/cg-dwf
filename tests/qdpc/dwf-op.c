#define QDP_Nc 3
#define QDP_Precision 'D'

#include <stdarg.h>
#include <stdio.h>
#include <qdp.h>
#include <qla.h>
#include <qmp.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#ifndef PREFIX
#define PREFIX MIT_altivecf_
#endif

#include "../../dwf-ssef.h"
#include "../../dwf-ssed.h"
#include "../../dwf-altivecf.h"
#include "../../dwf-bluelightf.h"
#include "../../dwf-bluelightd.h"

#define __L3(a,b)  a##b
#define _L3(a,b)   __L3(a,b)
#define L3(n)      _L3(PREFIX,n)

static int src_pos[5];
static int src_c;
static int src_d;
static int src_ri;


static void
error(const char *fmt, ...)
{
    if (QDP_this_node == 0) {
	va_list va;

	fflush(stdout);
	fprintf(stderr, "FATAL ERROR:");
	va_start(va, fmt);
	vfprintf(stderr, fmt, va);
	va_end(va);
	fflush(stderr);
    }
    QDP_abort();
}

static int lattice_size[5];

static void
initial_iseed(QLA_Int *li, int coords[])
{
  int i,t;

  t = coords[0];
  for(i=1; i<4; ++i) {
    t = t*lattice_size[i] + coords[i];
  }

  *li = t;
}

static void
print(const char *fmt, ...)
{
    va_list va;

    if (QDP_this_node != 0)
	return;
    
    va_start(va, fmt);
    vprintf(fmt, va);
    va_end(va);
}

/*****************************************************************************/
static double
gauge_reader(const void *ptr, void *env, 
	     const int latt_coord[], int mu, int row, int col, int reim)
{
    QLA_ColorMatrix **u = (QLA_ColorMatrix **)ptr;
    QLA_Complex       c;
    QLA_Real          v;
    int               node = QDP_node_number((int *)latt_coord);
    int               linear = QDP_index((int *)latt_coord);

    if (node != QDP_this_node) {
	error("Coordinates [%d %d %d %d] are on node %d, not on %d\n",
	      latt_coord[0], latt_coord[1], latt_coord[2], latt_coord[3],
	      node, QDP_this_node);
    }

    QLA_C_eq_elem_M(&c, &u[mu][linear], row, col);
    if (reim == 0)
	QLA_R_eq_re_C(&v, &c);
    else
	QLA_R_eq_im_C(&v, &c);

    return v;
}

static double
fermion_reader(const void *ptr, void *env, 
	       const int latt_coord[5], int color, int spin, int reim)
{
    int i;

    for (i = 0; i < 5; i++) {
	if (latt_coord[i] != src_pos[i])
	    return 0;
    }
    if (color != src_c) 
	return 0;
    if (spin != src_d)
	return 0;
    if (reim != src_ri)
	return 0;
    return 1;
}

static void
fermion_writer(void *pre, void *env,
	       const int latt_coord[5], int color, int spin, int reim,
	       double value)
{
    if (value == 0)
	return;
    printf("[%05d]: w [%2d %2d %2d %2d %2d] (%d,%d).%s=%g\n",
	   QDP_this_node,
	   latt_coord[0],
	   latt_coord[1],
	   latt_coord[2],
	   latt_coord[3],
	   latt_coord[4],
	   color,
	   spin,
	   reim?"im":"re",
	   value);
}

static void
operator(QDP_ColorMatrix *xU[],
	 double M5,
	 double m_f)
{
    QDP_ColorMatrix *xV[4];
    QLA_ColorMatrix *U[4];
    QLA_ColorMatrix *V[4];
    L3(DWF_Gauge)   *g;
    L3(DWF_Fermion) *x;
    L3(DWF_Fermion) *Ax;
    int              i;
    double           M_0;
    
    for (i = 0; i < 4; i++) {
	xV[i] = QDP_create_M();
	QDP_M_eq_sM(xV[i], xU[i], QDP_neighbor[i], QDP_forward, QDP_all);
    }

    QDP_suspend_comm();

    for (i = 0; i < 4; i++) {
	U[i] = QDP_expose_M(xU[i]);
	V[i] = QDP_expose_M(xV[i]);
    }

    if (L3(DWF_init)(lattice_size, NULL, NULL)) {
	error("DWF init() failed");
    }

    g   = L3(DWF_load_gauge)(U, V, NULL, gauge_reader);
    x = L3(DWF_load_fermion)(NULL, NULL, fermion_reader);
    Ax = L3(DWF_allocate_fermion)();

    M_0 = -2 * (5.0 - M5);
    
    L3(DWF_Dirac_Operator)(Ax, g, M_0, m_f, x);

    L3(DWF_save_fermion)(NULL, NULL, fermion_writer, Ax);

    L3(DWF_fini)();
    
    QDP_resume_comm();
}

/*****************************************************************************/

int
main(int argc, char *argv[])
{
    QLA_Complex cone = { 1, 0 };
    QDP_Int           *iseed;
    QDP_RandomState   *rs;
    int                i;
    QDP_ColorMatrix   *U[4];
    double             M5;
    double             m_f;

    QDP_initialize(&argc, &argv);
    
    if (argc != 16) {
	print("usage: "
	      "%s Lx Ly Lz Lt Ls  M5 m_f  px py pz pt ps color dirac re/im\n",
	      argv[0]);
	QDP_finalize();
	return 1;
    }

    for (i = 0; i < 5; i++)
	lattice_size[i] = atoi(argv[1 + i]);

    M5       = atof(argv[6]);
    m_f      = atof(argv[7]);

    for (i = 0; i < 5; i++)
	src_pos[i] = atoi(argv[8 + i]);
    src_c = atoi(argv[13]);
    src_d = atoi(argv[14]);
    src_ri = atoi(argv[15]);

    print("lattice %d %d %d %d %d\n",
	  lattice_size[0],
	  lattice_size[1],
	  lattice_size[2],
	  lattice_size[3],
	  lattice_size[4]);


    print("M5         = %g\n", M5);
    print("m_f        = %g\n", m_f);

    print("source: %d %d %d %d %d  %d %d %s\n",
	  src_pos[0],
	  src_pos[1],
	  src_pos[2],
	  src_pos[3],
	  src_pos[4],
	  src_c, src_d, src_ri?"im": "re");

    QDP_set_latsize(4, lattice_size);
    QDP_create_layout();
    QDP_check_comm(1);

    print("sites per node %d\n", QDP_sites_on_node);
    {
      int d = QMP_get_logical_number_of_dimensions();
      const int *s = QMP_get_logical_dimensions();

      if (s)
	print("block [%d %d %d %d]\n",
	      d < 0? -1: (int)s[0],
	      d < 1? -1: (int)s[1],
	      d < 2? -1: (int)s[2],
	      d < 3? -1: (int)s[3]);
    }



    iseed = QDP_create_I();
    rs = QDP_create_S();
    QDP_I_eq_func(iseed, initial_iseed, QDP_all);
    QDP_S_eq_seed_i_I(rs, 1235470343, iseed, QDP_all);

    for (i = 0; i < 4; i++) {
	U[i] = QDP_create_M();
        QDP_M_eq_c( U[i], &cone, QDP_all );
    }

    operator(U, M5, m_f);
    

    QDP_finalize();
    return 0;
}
