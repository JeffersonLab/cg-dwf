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

#if 0
    printf("gauge_reader: lat[%d %d %d %d], mu %d, r %d, c %d, ri %d, v %g\n",
	   latt_coord[0],
	   latt_coord[1],
	   latt_coord[2],
	   latt_coord[3],
	   mu,
	   row,
	   col,
	   reim,
	   v);
#endif
    return v;
}

static double
fermion_reader_guess(const void *ptr, void *env, 
		     const int latt_coord[5], int color, int spin, int reim)
{
    QLA_DiracFermion **psi    = (QLA_DiracFermion **)ptr;
    QLA_Complex        c;
    QLA_Real           v;
    int                Ls1    = lattice_size[4] - 1;
    int                s      = latt_coord[4];
    int                node   = QDP_node_number((int *)latt_coord);
    int                linear = QDP_index((int *)latt_coord);

    if (node != QDP_this_node) {
	error("Coordinates [%d %d %d %d %d] are on node %d, not on %d\n",
	      latt_coord[0], latt_coord[1], latt_coord[2], latt_coord[3],
	      latt_coord[4],
	      node, QDP_this_node);
    }
    QLA_C_eq_elem_D(&c, &psi[Ls1-s][linear], color, spin);
    if (reim == 0)
	QLA_R_eq_re_C(&v, &c);
    else
	QLA_R_eq_im_C(&v, &c);

    if (spin >= 2)
	v = -v;

#if 0
    printf("fermion_guess: lat[%d %d %d %d %d], c %d, s %d, ri %d, v %g\n",
	   latt_coord[0],
	   latt_coord[1],
	   latt_coord[2],
	   latt_coord[3],
	   latt_coord[4],
	   color,
	   spin,
	   reim,
	   v);
#endif

    return -0.5 * v;
}


static void
test(QDP_ColorMatrix *xU[],
     double M5,
     double m_f,
     QDP_DiracFermion **a,
     QDP_DiracFermion **b)
{
    int i;
    double M_0;
    QDP_ColorMatrix *xV[4];
    QLA_ColorMatrix *U[4];
    QLA_ColorMatrix *V[4];
    QLA_DiracFermion **xa;
    QLA_DiracFermion **xb;
    L3(DWF_Gauge)   *g;
    L3(DWF_Fermion) *A;
    L3(DWF_Fermion) *B;
    L3(DWF_Fermion) *DA;
    L3(DWF_Fermion) *BD;
    double bD_a_re, bD_a_im;
    double b_Da_re, b_Da_im;
    
    for (i = 0; i < 4; i++) {
	xV[i] = QDP_create_M();
	QDP_M_eq_sM(xV[i], xU[i], QDP_neighbor[i], QDP_backward, QDP_all);
    }

    QDP_suspend_comm();

    for (i = 0; i < 4; i++) {
	U[i] = QDP_expose_M(xU[i]);
	V[i] = QDP_expose_M(xV[i]);
    }

    xa = malloc(lattice_size[4] * sizeof (QLA_DiracFermion *));
    xb  = malloc(lattice_size[4] * sizeof (QLA_DiracFermion *));
    for (i = 0; i < lattice_size[4]; i++) {
	xa[i] = QDP_expose_D(a[i]);
	xb[i] = QDP_expose_D(b[i]);
    }


    if (L3(DWF_init)(lattice_size, NULL, NULL)) {
	error("DWF init() failed");
    }

    g   = L3(DWF_load_gauge)(U, V, NULL, gauge_reader);
    A   = L3(DWF_load_fermion)(xa, NULL, fermion_reader_guess);
    B   = L3(DWF_load_fermion)(xb, NULL, fermion_reader_guess);
    DA  = L3(DWF_allocate_fermion)();
    BD  = L3(DWF_allocate_fermion)();

    M_0 = -2*(5.0-M5) ;  
    L3(DWF_Dirac_Operator)(DA, g, M_0, m_f, A);
    L3(DWF_Dirac_Operator_conjugate)(BD, g, M_0, m_f, B);
    L3(DWF_Fermion_Dot_Product)(&bD_a_re, &bD_a_im, BD, A);
    L3(DWF_Fermion_Dot_Product)(&b_Da_re, &b_Da_im, B, DA);

    print("D+ vs. D bD_a: %g %g\n", bD_a_re, bD_a_im);
    print("D+ vs. D b_Da: %g %g\n", b_Da_re, b_Da_im);

    L3(DWF_fini)();
    
    QDP_resume_comm();
    
    for (i = 0; i < lattice_size[4]; i++) {
	QDP_reset_D(a[i]);
	QDP_reset_D(b[i]);
    }

    for (i = 0; i < 4; i++) {
	QDP_reset_M(xU[i]);
	QDP_reset_M(xV[i]);
	QDP_destroy_M(xV[i]);
    }
}

/*****************************************************************************/

int
main(int argc, char *argv[])
{
    int                i;
    QDP_Int           *iseed;
    QDP_RandomState   *rs;
    QDP_ColorMatrix   *U[4];
    QDP_DiracFermion **a;
    QDP_DiracFermion **b;
    double             M5;
    double             m_f;

    QDP_initialize(&argc, &argv);
    
    if (argc != 8) {
	print("usage: %s Lx Ly Lz Lt Ls M5 m_f\n", argv[0]);
	QDP_finalize();
	return 1;
    }

    for (i = 0; i < 5; i++)
	lattice_size[i] = atoi(argv[1 + i]);

    M5       = atof(argv[6]);
    m_f      = atof(argv[7]);

    print("Test of D and D^+ on random U\n");

    print("lattice %d %d %d %d %d\n",
	  lattice_size[0],
	  lattice_size[1],
	  lattice_size[2],
	  lattice_size[3],
	  lattice_size[4]);

    print("M5         = %g\n", M5);
    print("m_f        = %g\n", m_f);
    
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
	QDP_M_eq_gaussian_S(U[i], rs, QDP_all);
    }

    a = malloc(lattice_size[4] * sizeof (QDP_DiracFermion *));
    b = malloc(lattice_size[4] * sizeof (QDP_DiracFermion *));

    for (i = 0; i < lattice_size[4]; i++) {
	a[i] = QDP_create_D();
	QDP_D_eq_gaussian_S(a[i], rs, QDP_all);
	b[i] = QDP_create_D();
	QDP_D_eq_gaussian_S(b[i], rs, QDP_all);
    }

    /* Do what it takes */
    test(U, M5, m_f, a, b);

    QDP_finalize();
    return 0;
}
