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

void
xxx(void)
{
  fflush(stdout);
}

static double
diff_time(const struct rusage *r0,
	  const struct rusage *r1)
{
  return (r1->ru_utime.tv_sec - r0->ru_utime.tv_sec +
	  1e-6 * (r1->ru_utime.tv_usec - r0->ru_utime.tv_usec));
}

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
#if 0
    printf("fermion_rhs: lat[%d %d %d %d %d], c %d, s %d, ri %d, v %g\n",
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
solve(QDP_DiracFermion **solution,
      QDP_ColorMatrix *xU[],
      double M5,
      double m_f,
      QDP_DiracFermion **rhs,
      QDP_DiracFermion **x0,
      double eps,
      int iter)
{
    int status;
    int i;
    int out_iter;
    double out_eps;
    double M_0;
    QDP_ColorMatrix *xV[4];
    QLA_ColorMatrix *U[4];
    QLA_ColorMatrix *V[4];
    QLA_DiracFermion **xRHS;
    QLA_DiracFermion **xX0;
    L3(DWF_Gauge)   *g;
    L3(DWF_Fermion) *eta;
    L3(DWF_Fermion) *X0;
    L3(DWF_Fermion) *res;
    L3(DWF_Fermion) *Mr;
    L3(DWF_Fermion) *dlt;
    struct rusage r0, r1;
    double d_r, d_i;
    double r_r, r_i;
    double flops;
    double v4, ls;
    
    for (i = 0; i < 4; i++) {
	xV[i] = QDP_create_M();
	QDP_M_eq_sM(xV[i], xU[i], QDP_neighbor[i], QDP_forward, QDP_all);
    }

    QDP_suspend_comm();

    for (i = 0; i < 4; i++) {
	U[i] = QDP_expose_M(xU[i]);
	V[i] = QDP_expose_M(xV[i]);
    }

    xRHS = malloc(lattice_size[4] * sizeof (QLA_DiracFermion *));
    xX0  = malloc(lattice_size[4] * sizeof (QLA_DiracFermion *));
    for (i = 0; i < lattice_size[4]; i++) {
	xRHS[i] = QDP_expose_D(rhs[i]);
	if (x0 != rhs)
	  xX0[i] = QDP_expose_D(x0[i]);
	else
	  xX0[i] = xRHS[i];
    }


    if (L3(DWF_init)(lattice_size, NULL, NULL)) {
	error("DWF init() failed");
    }

    g   = L3(DWF_load_gauge)(U, V, NULL, gauge_reader);
    eta = L3(DWF_load_fermion)(xRHS, NULL, fermion_reader_rhs);
    X0  = L3(DWF_load_fermion)(xX0, NULL, fermion_reader_guess);
    res = L3(DWF_allocate_fermion)();
    Mr  = L3(DWF_allocate_fermion)();
    dlt = L3(DWF_allocate_fermion)();

    M_0 = -2*(5.0-M5) ;  
    getrusage(RUSAGE_SELF, &r0);
    status = L3(DWF_cg_solver)(res, &out_eps, &out_iter,
			       g, M_0, m_f, X0, eta, eps, 1, iter);
    getrusage(RUSAGE_SELF, &r1);

    ls = lattice_size[4];
    v4 = lattice_size[0] * lattice_size[1] * lattice_size[2] * lattice_size[3];
    flops = v4 * (out_iter * (ls * 3312. + 48.) + ls * 9708. + 144.);
    
    L3(DWF_Dirac_Operator)(Mr, g, M_0, m_f, res);
    L3(DWF_Add_Fermion)(dlt, Mr, -1.0, eta);
    L3(DWF_Fermion_Dot_Product)(&d_r, &d_i, dlt, dlt);
    L3(DWF_Fermion_Dot_Product)(&r_r, &r_i, eta, eta);
    
    print("L3 DWF solver: status = %d, iterations=%d, epsilon=%g\n",
	  status, out_iter, out_eps);
    print("L3 DWF time: %g sec, %d iterations, %d Ls\n",
	  diff_time(&r0, &r1), out_iter, lattice_size[4]);
    print("L3 DWF true residue squared: %g\n", d_r);
    print("L3 DWF rhs normalization: %g\n", r_r);
    print("L3 DWF residue normalized: %g\n", d_r / r_r);
    {
	int d = QMP_get_logical_number_of_dimensions();
	const int *s = QMP_get_logical_dimensions();

	print("L3 DWF perf: %g Mflops, l %d %d %d %d %d  n %d %d %d %d\n",
	      flops * 1e-6 / diff_time(&r0, &r1),
	      lattice_size[0],
	      lattice_size[1],
	      lattice_size[2],
	      lattice_size[3],
	      lattice_size[4],
	      d < 0? 1: s[0],
	      d < 1? 1: s[1],
	      d < 2? 1: s[2],
	      d < 3? 1: s[3]);
    }

    /* Free all allocated L3 fields here */

    L3(DWF_fini)();
    
    QDP_resume_comm();
    
    for (i = 0; i < lattice_size[4]; i++) {
        if (x0 != rhs)
	    QDP_reset_D(x0[i]);
	QDP_reset_D(rhs[i]);
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
    QLA_Complex cone = { 1, 0 };
    QDP_Int           *iseed;
    QDP_RandomState   *rs;
    int                i;
    QDP_ColorMatrix   *U[4];
    QDP_DiracFermion **rhs;
    QDP_DiracFermion **solution;
    double             M5;
    double             m_f;
    double             epsilon;
    int                iter;

    QDP_initialize(&argc, &argv);
    
    if (argc != 10) {
	print("usage: %s Lx Ly Lz Lt Ls M5 m_f epsilon iter\n", argv[0]);
	QDP_finalize();
	return 1;
    }

    for (i = 0; i < 5; i++)
	lattice_size[i] = atoi(argv[1 + i]);

    M5       = atof(argv[6]);
    m_f      = atof(argv[7]);
    epsilon  = atof(argv[8]);
    iter     = atoi(argv[9]);


    print("lattice %d %d %d %d %d\n",
	  lattice_size[0],
	  lattice_size[1],
	  lattice_size[2],
	  lattice_size[3],
	  lattice_size[4]);


    print("M5         = %g\n", M5);
    print("m_f        = %g\n", m_f);
    print("epsilon    = %g\n", epsilon);
    print("iterations = %d\n", iter);
    
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

    rhs      = malloc(lattice_size[4] * sizeof (QDP_DiracFermion *));
    solution = malloc(lattice_size[4] * sizeof (QDP_DiracFermion *));

    for (i = 0; i < lattice_size[4]; i++) {
	solution[i] = QDP_create_D();
	rhs[i]      = QDP_create_D();
	QDP_D_eq_gaussian_S(rhs[i], rs, QDP_all);
    }

    /* Do what it takes */
    solve(solution, U, M5, m_f, rhs, rhs, epsilon, iter);
    

    QDP_finalize();
    return 0;
}
