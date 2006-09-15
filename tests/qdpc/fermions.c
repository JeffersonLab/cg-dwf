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
static void
the_dot(const char *name, 
	L3(DWF_Fermion) *x,
	L3(DWF_Fermion) *y)
{
  double d_re, d_im;

  L3(DWF_Fermion_Dot_Product)(&d_re, &d_im, x, y);
  print("dot product (%s): %g %g\n", name, d_re, d_im);
}

static double
fermion_reader_rhs(const void *ptr, void *env, 
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
    return v;
}

static void
test(QDP_DiracFermion **a,
     QDP_DiracFermion **b,
     QDP_DiracFermion **c,
     double factor)
{
    int i;
    QLA_DiracFermion **xA;
    QLA_DiracFermion **xB;
    QLA_DiracFermion **xC;
    L3(DWF_Fermion) *A;
    L3(DWF_Fermion) *B;
    L3(DWF_Fermion) *C;

    QDP_suspend_comm();

    xA = malloc(lattice_size[4] * sizeof (QLA_DiracFermion *));
    xB = malloc(lattice_size[4] * sizeof (QLA_DiracFermion *));
    xC = malloc(lattice_size[4] * sizeof (QLA_DiracFermion *));
    for (i = 0; i < lattice_size[4]; i++) {
	xA[i] = QDP_expose_D(a[i]);
	xB[i] = QDP_expose_D(b[i]);
	xC[i] = QDP_expose_D(c[i]);
    }

    if (L3(DWF_init)(lattice_size, NULL, NULL)) {
	error("DWF init() failed");
    }

    A = L3(DWF_load_fermion)(xA, NULL, fermion_reader_rhs);
    B = L3(DWF_load_fermion)(xB, NULL, fermion_reader_rhs);
    C = L3(DWF_load_fermion)(xC, NULL, fermion_reader_rhs);

    L3(DWF_Add_Fermion)(C, A, factor, B);
    the_dot("a2", A, A);
    the_dot("ab", A, B);
    the_dot("ba", B, A);
    the_dot("ac", A, C);
    the_dot("ca", C, A);
    the_dot("b2", B, B);
    the_dot("bc", B, C);
    the_dot("cb", C, B);
    the_dot("c2", C, C);

    L3(DWF_fini)();
    
    QDP_resume_comm();
    
    for (i = 0; i < lattice_size[4]; i++) {
	QDP_reset_D(a[i]);
	QDP_reset_D(b[i]);
	QDP_reset_D(c[i]);
    }
}

/*****************************************************************************/

int
main(int argc, char *argv[])
{
    int                i;
    QDP_Int           *iseed;
    QDP_RandomState   *rs;
    QDP_DiracFermion **a;
    QDP_DiracFermion **b;
    QDP_DiracFermion **c;
    double             factor;

    QDP_initialize(&argc, &argv);
    
    if (argc != 7) {
	print("usage: %s Lx Ly Lz Lt Ls factor\n", argv[0]);
	QDP_finalize();
	return 1;
    }

    for (i = 0; i < 5; i++)
	lattice_size[i] = atoi(argv[1 + i]);

    factor = atof(argv[6]);

    print("Fermion helper functions test\n");
    print("lattice %d %d %d %d %d\n",
	  lattice_size[0],
	  lattice_size[1],
	  lattice_size[2],
	  lattice_size[3],
	  lattice_size[4]);

    print("factor     = %g\n", factor);
    
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

    a      = malloc(lattice_size[4] * sizeof (QDP_DiracFermion *));
    b      = malloc(lattice_size[4] * sizeof (QDP_DiracFermion *));
    c      = malloc(lattice_size[4] * sizeof (QDP_DiracFermion *));

    for (i = 0; i < lattice_size[4]; i++) {
	a[i] = QDP_create_D();
	QDP_D_eq_gaussian_S(a[i], rs, QDP_all);

	b[i] = QDP_create_D();
	QDP_D_eq_gaussian_S(b[i], rs, QDP_all);

	c[i] = QDP_create_D();
	QDP_D_eq_gaussian_S(c[i], rs, QDP_all);
    }

    /* Do what it takes */
    test(a, b, c, factor);

    QDP_finalize();
    return 0;
}
