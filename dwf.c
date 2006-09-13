#line 641 "dwf.nw"
#include <string.h>
#include <stdlib.h>
#line 1938 "dwf.nw"
#include <qmp.h>
#line 3003 "dwf.nw"
#define COMPUTE_YA(j,t,c) {                \
   REAL *v_re = (REAL *)&rs->f[j][c].re;   \
   REAL *v_im = (REAL *)&rs->f[j][c].im;   \
   BLOCKOF_YA(j,t,c,re)                    \
   BLOCKOF_YA(j,t,c,im) }
#line 3011 "dwf.nw"
#define COMPUTE_YB(j,t,c) {                \
   REAL *v_re = (REAL *)&rs->f[j][c].re;   \
   REAL *v_im = (REAL *)&rs->f[j][c].im;   \
   BLOCKOF_YB(j,t,c,re)                    \
   BLOCKOF_YB(j,t,c,im) }
#line 3023 "dwf.nw"
#define BLOCKOF2_YA(j,t,c,ri)                                             \
  v_##ri[1] = inv_a * v_##ri[1] + b_over_a * yOut[t][c].ri;               \
  yOut[t][c].ri = v_##ri[0] = inv_a * v_##ri[0] + b_over_a * v_##ri[1];
#line 3029 "dwf.nw"
#define BLOCKOF4_YA(j,t,c,ri)                                             \
  v_##ri[3] = inv_a * v_##ri[3] + b_over_a * yOut[t][c].ri;               \
  v_##ri[2] = inv_a * v_##ri[2] + b_over_a * v_##ri[3];                   \
  v_##ri[1] = inv_a * v_##ri[1] + b_over_a * v_##ri[2];                   \
  yOut[t][c].ri = v_##ri[0] = inv_a * v_##ri[0] + b_over_a * v_##ri[1];
#line 3043 "dwf.nw"
#define BLOCKOF2_YB(j,t,c,ri)                                             \
  v_##ri[0] = inv_a * v_##ri[0] + b_over_a * yOut[t][c].ri;               \
  yOut[t][c].ri = v_##ri[1] = inv_a * v_##ri[1] + b_over_a * v_##ri[0];
#line 3049 "dwf.nw"
#define BLOCKOF4_YB(j,t,c,ri)                                             \
  v_##ri[0] = inv_a * v_##ri[0] + b_over_a * yOut[t][c].ri;               \
  v_##ri[1] = inv_a * v_##ri[1] + b_over_a * v_##ri[0];                   \
  v_##ri[2] = inv_a * v_##ri[2] + b_over_a * v_##ri[1];                   \
  yOut[t][c].ri = v_##ri[3] = inv_a * v_##ri[3] + b_over_a * v_##ri[2];
#line 4165 "dwf.nw"
#ifdef DEBUG_QMP
#undef DEBUG_QMP
#define DEBUG_QMP(msg,a ...) do \
     printf("[%05d] %s:%d:QMP/%s(): " msg, QMP_get_node_number(), \
                             __FILE__,__LINE__,__FUNCTION__,##a); \
  while(0);
#else /* !defined(DEBUG_QMP) */
#define DEBUG_QMP(msg,a ...)
#define NO_DEBUG_QMP
#endif /* defined(DEBUG_QMP) */
#line 4178 "dwf.nw"
#ifdef DEBUG_DWF
#undef DEBUG_DWF
#define DEBUG_DWF(msg,a ...) do \
     printf("[%05d] %s: %d: DWF/%s(): " msg, QMP_get_node_number(), \
                               __FILE__,__LINE__,__FUNCTION__,##a); \
  while(0);
#else /* !defined(DEBUG_DWF) */
#define DEBUG_DWF(msg, a ...)
#define NO_DEBUG_DWF
#endif /* defined(DEBUG_DWF) */
#line 801 "dwf.nw"
struct L3(DWF_Gauge) {
    scalar_complex v[Nc][Nc];
};
#line 1024 "dwf.nw"
typedef struct {
    vector_complex f[Fd][Nc];
} vFermion;
#line 1031 "dwf.nw"
typedef struct {
    vector_complex f[Fd/2][Nc];
} vHalfFermion;
#line 1038 "dwf.nw"
typedef struct {
    scalar_complex v[Nc][Nc];
} SU3;
#line 1045 "dwf.nw"
typedef struct {
    vector_complex v[Nc][Nc];
} vSU3;
#line 1052 "dwf.nw"
typedef struct {
    vFermion f;
} vEvenFermion;

typedef struct {
    vFermion f;
} vOddFermion;
#line 1063 "dwf.nw"
struct L3(DWF_Fermion) {
    vEvenFermion *even;
    vOddFermion *odd;
};
#line 1128 "dwf.nw"
struct memblock {
    struct memblock *next;
    struct memblock *prev;
    void *data;
    size_t size;
};
#line 1434 "dwf.nw"
struct bounds {
    int lo[DIM];
    int hi[DIM];
};
#line 1550 "dwf.nw"
struct neighbor {
    int              size;            /* size of site table */
    int              inside_size;     /* number of inside sites */
    int              boundary_size;   /* number of boundary sites */
    int              snd_size[2*DIM]; /* size of send buffers in 8 dirs */
    int              rcv_size[2*DIM]; /* size of receive buffers */
    int             *snd[2*DIM];      /* i->x translation for send buffers */
    int             *inside;          /* i->x translation for inside sites */
    struct boundary *boundary;        /* i->x,mask translation for boundary */
    struct site     *site;            /* x->site translation for sites */
    vHalfFermion    *snd_buf[2*DIM];  /* Send buffers */
    vHalfFermion    *rcv_buf[2*DIM];  /* Receive buffers */

    int              qmp_size[4*DIM]; /* sizes of QMP buffers */
    void            *qmp_xbuf[4*DIM]; /* QMP snd/rcv buffer addresses */
    vHalfFermion    *qmp_buf[4*DIM];  /* send and receive buffers for QMP */
    QMP_msgmem_t     qmp_mm[4*DIM];   /* msgmem's for send and receive */
    int              Nx;              /* number of msegs */

    int              qmp_smask;       /* send flags for qmp_sh[] */
    QMP_msghandle_t  qmp_handle;      /* common send & receive handle */
};
#line 1576 "dwf.nw"
struct boundary {
    int   index;   /* x-index of this boundary site */
    int   mask;    /* bitmask of the borders */
};
#line 1585 "dwf.nw"
struct site {
  int Uup;           /* up-links are Uup, Uup+1, Uup+2, Uup+3 */
  int Udown[DIM];    /* four down-links */
  int F[2*DIM];      /* eight neighboring fermions on the other sublattice */
};
#line 513 "dwf.nw"
static vReal vcbn;
static vReal vbnc;
static vReal vb;
static vReal va;
#line 588 "dwf.nw"
static int inited_p = 0;
#line 636 "dwf.nw"
static void *(*tmalloc)(size_t size);
static void (*tfree)(void *ptr);
#line 665 "dwf.nw"
static int tlattice[DIM+1];
#line 850 "dwf.nw"
static SU3 *U;
#line 1119 "dwf.nw"
static struct memblock memblock = {
    &memblock,
    &memblock,
    NULL,
    0
};
#line 1201 "dwf.nw"
static int network[DIM];
static int coord[DIM];
#line 1414 "dwf.nw"
static struct bounds bounds;
static int gauge_XYZT;
static int Sv, Sv_1;
#line 1444 "dwf.nw"
static struct neighbor neighbor;
#ifndef DEBUG_QMP
static struct neighbor odd_even;
static struct neighbor even_odd;
#else
struct neighbor odd_even;
struct neighbor even_odd;
#endif
#line 2192 "dwf.nw"
static vOddFermion *auxA_o, *auxB_o, *Phi_o;
static vEvenFermion *auxA_e;
#line 2270 "dwf.nw"
static vOddFermion *r_o, *p_o, *q_o;
#line 2291 "dwf.nw"
static vEvenFermion *auxB_e;
#line 2653 "dwf.nw"
static REAL inv_a;
static REAL b_over_a;
#line 2728 "dwf.nw"
static vReal vfx_A;
static vReal vfx_B;
static vReal vab;
static REAL c0;
#line 913 "dwf.nw"
static inline void
collect_add(vFermion *r,
            const vFermion *x,
	    double A,
            const vFermion *y,
	    int n)
{
  int i;
  vReal a = vmk_1(A);

  for (i = 0; i < n; i++) {
     r->f[0][0].re = x->f[0][0].re + a * y->f[0][0].re;
     r->f[0][0].im = x->f[0][0].im + a * y->f[0][0].im;
     r->f[0][1].re = x->f[0][1].re + a * y->f[0][1].re;
     r->f[0][1].im = x->f[0][1].im + a * y->f[0][1].im;
     r->f[0][2].re = x->f[0][2].re + a * y->f[0][2].re;
     r->f[0][2].im = x->f[0][2].im + a * y->f[0][2].im;

     r->f[1][0].re = x->f[1][0].re + a * y->f[1][0].re;
     r->f[1][0].im = x->f[1][0].im + a * y->f[1][0].im;
     r->f[1][1].re = x->f[1][1].re + a * y->f[1][1].re;
     r->f[1][1].im = x->f[1][1].im + a * y->f[1][1].im;
     r->f[1][2].re = x->f[1][2].re + a * y->f[1][2].re;
     r->f[1][2].im = x->f[1][2].im + a * y->f[1][2].im;

     r->f[2][0].re = x->f[2][0].re + a * y->f[2][0].re;
     r->f[2][0].im = x->f[2][0].im + a * y->f[2][0].im;
     r->f[2][1].re = x->f[2][1].re + a * y->f[2][1].re;
     r->f[2][1].im = x->f[2][1].im + a * y->f[2][1].im;
     r->f[2][2].re = x->f[2][2].re + a * y->f[2][2].re;
     r->f[2][2].im = x->f[2][2].im + a * y->f[2][2].im;

     r->f[3][0].re = x->f[3][0].re + a * y->f[3][0].re;
     r->f[3][0].im = x->f[3][0].im + a * y->f[3][0].im;
     r->f[3][1].re = x->f[3][1].re + a * y->f[3][1].re;
     r->f[3][1].im = x->f[3][1].im + a * y->f[3][1].im;
     r->f[3][2].re = x->f[3][2].re + a * y->f[3][2].re;
     r->f[3][2].im = x->f[3][2].im + a * y->f[3][2].im;
  }
}
#line 974 "dwf.nw"
static inline void
collect_dot(double *r_re,
            double *r_im,
	    const vFermion *a,
	    const vFermion *b,
	    int n)
{
  int i;
  vReal c0_re, c1_re, c2_re;
  vReal c0_im, c1_im, c2_im;

  for (i = 0; i < n; i++, a++, b++) {
    c0_re = a->f[0][0].re*b->f[0][0].re; c0_re += a->f[0][0].im*b->f[0][0].im;
    c0_im = a->f[0][0].im*b->f[0][0].re; c0_im += a->f[0][0].re*b->f[0][0].im;
    c1_re = a->f[0][1].re*b->f[0][1].re; c1_re += a->f[0][1].im*b->f[0][1].im;
    c1_im = a->f[0][1].im*b->f[0][1].re; c1_im += a->f[0][1].re*b->f[0][1].im;
    c2_re = a->f[0][2].re*b->f[0][2].re; c2_re += a->f[0][2].im*b->f[0][2].im;
    c2_im = a->f[0][2].im*b->f[0][2].re; c2_im += a->f[0][2].re*b->f[0][2].im;

    c0_re += a->f[1][0].re*b->f[1][0].re; c0_re += a->f[1][0].im*b->f[1][0].im;
    c0_im += a->f[1][0].im*b->f[1][0].re; c0_im += a->f[1][0].re*b->f[1][0].im;
    c1_re += a->f[1][1].re*b->f[1][1].re; c1_re += a->f[1][1].im*b->f[1][1].im;
    c1_im += a->f[1][1].im*b->f[1][1].re; c1_im += a->f[1][1].re*b->f[1][1].im;
    c2_re += a->f[1][2].re*b->f[1][2].re; c2_re += a->f[1][2].im*b->f[1][2].im;
    c2_im += a->f[1][2].im*b->f[1][2].re; c2_im += a->f[1][2].re*b->f[1][2].im;

    c0_re += a->f[2][0].re*b->f[2][0].re; c0_re += a->f[2][0].im*b->f[2][0].im;
    c0_im += a->f[2][0].im*b->f[2][0].re; c0_im += a->f[2][0].re*b->f[2][0].im;
    c1_re += a->f[2][1].re*b->f[2][1].re; c1_re += a->f[2][1].im*b->f[2][1].im;
    c1_im += a->f[2][1].im*b->f[2][1].re; c1_im += a->f[2][1].re*b->f[2][1].im;
    c2_re += a->f[2][2].re*b->f[2][2].re; c2_re += a->f[2][2].im*b->f[2][2].im;
    c2_im += a->f[2][2].im*b->f[2][2].re; c2_im += a->f[2][2].re*b->f[2][2].im;

    c0_re += a->f[3][0].re*b->f[3][0].re; c0_re += a->f[3][0].im*b->f[3][0].im;
    c0_im += a->f[3][0].im*b->f[3][0].re; c0_im += a->f[3][0].re*b->f[3][0].im;
    c1_re += a->f[3][1].re*b->f[3][1].re; c1_re += a->f[3][1].im*b->f[3][1].im;
    c1_im += a->f[3][1].im*b->f[3][1].re; c1_im += a->f[3][1].re*b->f[3][1].im;
    c2_re += a->f[3][2].re*b->f[3][2].re; c2_re += a->f[3][2].im*b->f[3][2].im;
    c2_im += a->f[3][2].im*b->f[3][2].re; c2_im += a->f[3][2].re*b->f[3][2].im;

    *r_re += vsum(c0_re + c1_re + c2_re);
    *r_im += vsum(c0_im + c1_im + c2_im);
  }
}
#line 1138 "dwf.nw"
static vEvenFermion *allocate_even_fermion(void);
static vOddFermion *allocate_odd_fermion(void);
static L3(DWF_Gauge) *allocate_gauge_field(void);
#line 1322 "dwf.nw"
static inline vReal
import_vector(const void *z, void *env, L3(DWF_fermion_reader) reader,
              int x[DIM+1], int c, int d, int re_im)
{
    vReal f;
    REAL *v = (REAL *)&f;
    int i, xs;
    
    for (xs = x[DIM], i = 0; i < Vs; i++, x[DIM]++) {
        *v++ = reader(z, env, x, c, d, re_im);
    }
    x[DIM] = xs;
    return f;
}
#line 1371 "dwf.nw"
static inline void
save_vector(void *z, void *env, L3(DWF_fermion_writer) writer,
            int x[DIM+1], int c, int d, int re_im, vReal *f)
{
    REAL *v = (REAL *)f;
    int i, xs;
    
    for (xs = x[DIM], i = 0; i < Vs; i++, x[DIM]++) {
        writer(z, env, x, c, d, re_im, *v++);
    }
    x[DIM] = xs;
}
#line 1455 "dwf.nw"
static inline int
lattice_start(int lat, int net, int coord)
{
    int q = lat / net;
    int r = lat % net;

    return coord * q + ((coord < r)? coord: r);
}

static inline void
mk_sublattice(struct bounds *bounds,
              int coord[])
{
    int i;

    for (i = 0; i < DIM; i++) {
        bounds->lo[i] = lattice_start(tlattice[i], network[i], coord[i]);
        bounds->hi[i] = lattice_start(tlattice[i], network[i], coord[i] + 1);
    }
}
#line 1479 "dwf.nw"
static void
init_neighbor(struct bounds *bounds, struct neighbor *neighbor);
#line 1608 "dwf.nw"
static void build_neighbor(struct neighbor *out,
	                   int              parity,
	                   struct neighbor *in);
#line 1716 "dwf.nw"
static void construct_rec(struct neighbor *out,
                          int par,
                          struct bounds *bounds,
                          int dir,
                          int step);
#line 1780 "dwf.nw"
static void construct_snd(struct neighbor *out,
                          struct neighbor *in,
                          int par,
                          struct bounds *bounds,
                          int dir,
                          int step);
#line 1867 "dwf.nw"
static int
to_HFlinear(int p[],
            struct bounds *b,
            int q,
            int z)
{
    int x, d;
    for (x = 0, d = DIM; d--;) {
        int v = p[d] + ((d == q)?z:0);
        int s = b->hi[d] - b->lo[d];
        if (v < 0)
           v += tlattice[d];
        if (v >= tlattice[d])
           v -= tlattice[d];
       x = x * s + v - b->lo[d];
    }
    return x / 2;
}
#line 1890 "dwf.nw"
static int
to_Ulinear(int p[],
           struct bounds *b,
           int q)
{
    int x, d;

    if ((q < 0) || (p[q] > b->lo[q]) || (network[q] < 2)) {
        
#line 1906 "dwf.nw"
for (x = 0, d = DIM; d--;) {
    int s = b->hi[d] - b->lo[d];
    int v = p[d] - ((q == d)?1:0);
    if (v < 0)
        v += tlattice[d];
    x = x * s + v - b->lo[d];
}
return DIM * x + ((q < 0)?0:q);
#line 1899 "dwf.nw"
    } else {
        
#line 1918 "dwf.nw"
int s0, v0;
for (d = 0, v0 = 1; d < DIM; d++)
     v0 *= b->hi[d] - b->lo[d];
for (d = 0, s0 = DIM * v0; d < q; d++) {
     if (network[d] < 2)
         continue;
     s0 += v0 / (b->hi[d] - b->lo[d]);
}
for (d = DIM, x = 0; d--;) {
    int s = b->hi[d] - b->lo[d];
    int v = p[d];

    if (d == q)
        continue;
    x = x * s + v - b->lo[d];
}
return s0 + x;
#line 1901 "dwf.nw"
    }
}
#line 1949 "dwf.nw"
static int build_buffers(struct neighbor *nb);
#line 2027 "dwf.nw"
static int make_buffer(struct neighbor *nb, int size);
#line 2048 "dwf.nw"
static int make_send(struct neighbor *nb, int k, int i, int d,
                     QMP_msghandle_t SRh[4*DIM], int Nsr);
#line 2070 "dwf.nw"
static int make_receive(struct neighbor *nb, int k, int i, int d,
                        QMP_msghandle_t SRh[4*DIM], int Nsr);
#line 2107 "dwf.nw"
static void sse_aligned_buffer(struct neighbor *nb, int k, int size);
#line 2137 "dwf.nw"
static void free_buffers(struct neighbor *nb);
#line 2214 "dwf.nw"
static int cg(vOddFermion *psi,
              const vOddFermion *b,
              const vOddFermion *x0,
              double epsilon,
              int min_iter,
              int max_iter,
              double *out_eps,
              int *out_iter);
#line 2304 "dwf.nw"
static void copy_o(vOddFermion *dst, const vOddFermion *src);
#line 2323 "dwf.nw"
static void compute_sum2_o(vOddFermion *dst, double alpha, const vOddFermion *src);
#line 2341 "dwf.nw"
static void compute_sum2x_o(vOddFermion *dst, const vOddFermion *src, double alpha);
#line 2360 "dwf.nw"
static void compute_sum_e(vEvenFermion *d,
                          const vEvenFermion *x, double alpha, const vEvenFermion *y);
static void compute_sum_o(vOddFermion *d,
                          const vOddFermion *x, double alpha, const vOddFermion *y);
#line 2402 "dwf.nw"
static void compute_sum_oN(vOddFermion *d, double *norm,
                           const vOddFermion *x, double alpha, const vOddFermion *y);
#line 2455 "dwf.nw"
static void compute_sum2_oN(vOddFermion *d, double *norm,
                            double alpha, const vOddFermion *y);
#line 2511 "dwf.nw"
static void compute_MxM(vOddFermion *eta, double *norm,
                        const vOddFermion *psi);
static void compute_M(vOddFermion *eta, double *norm,
                      const vOddFermion *psi);
static void compute_Mx(vOddFermion *eta,
                       const vOddFermion *psi);
#line 2552 "dwf.nw"
static void compute_Qxx1(vFermion *eta, const vFermion *psi, int xyzt);
static void inline compute_Qee1(vEvenFermion *eta, const vEvenFermion *psi)
{
    compute_Qxx1(&eta->f, &psi->f, even_odd.size);
}
#line 2560 "dwf.nw"
static void inline compute_Qoo1(vOddFermion *eta, const vOddFermion *psi)
{
    compute_Qxx1(&eta->f, &psi->f, odd_even.size);
}
#line 2584 "dwf.nw"
static void compute_Soo1(vOddFermion *eta, const vOddFermion *psi);
#line 3066 "dwf.nw"
static void compute_Qxy(vFermion *d, const vFermion *s, struct neighbor *nb);
#line 3070 "dwf.nw"
static void inline compute_Qoe(vOddFermion *d, const vEvenFermion *s)
{
    compute_Qxy(&d->f, &s->f, &odd_even);
}
#line 3077 "dwf.nw"
static void inline compute_Qeo(vEvenFermion *d, const vOddFermion *s)
{
    compute_Qxy(&d->f, &s->f, &even_odd);
}
#line 3084 "dwf.nw"
static void compute_1Sxy(vFermion *d,
                         const vFermion *q,
			 const vFermion *s,
			 struct neighbor *nb);
static void inline compute_1Soe(vOddFermion *d,
                                const vOddFermion *q,
                                const vEvenFermion *s)
{
    compute_1Sxy(&d->f, &q->f, &s->f, &odd_even);
}
#line 3269 "dwf.nw"
static void compute_Dx(vFermion *chi,
                       const vFermion *eta,
                       const vFermion *psi,
                       struct neighbor *nb);
static void inline compute_De(vEvenFermion *chi,
                              const vEvenFermion *eta,
                              const vOddFermion *psi)
{
    compute_Dx(&chi->f, &eta->f, &psi->f, &even_odd);
}
#line 3283 "dwf.nw"
static void inline compute_Do(vOddFermion *chi,
                              const vOddFermion *eta,
                              const vEvenFermion *psi)
{
    compute_Dx(&chi->f, &eta->f, &psi->f, &odd_even);
}
#line 3293 "dwf.nw"
static void compute_Dcx(vFermion *chi,
                        const vFermion *eta,
                        const vFermion *psi,
                        struct neighbor *nb);
static void inline compute_Dce(vEvenFermion *chi,
                               const vEvenFermion *eta,
                               const vOddFermion *psi)
{
    compute_Dcx(&chi->f, &eta->f, &psi->f, &even_odd);
}
#line 3307 "dwf.nw"
static void inline compute_Dco(vOddFermion *chi,
                               const vOddFermion *eta,
                               const vEvenFermion *psi)
{
    compute_Dcx(&chi->f, &eta->f, &psi->f, &odd_even);
}
#line 3984 "dwf.nw"
static void compute_Qxx1Qxy(vFermion *d,
                            const vFermion *s,
                            struct neighbor *nb);
static void inline compute_Qee1Qeo(vEvenFermion *d, const vOddFermion *s)
{
  compute_Qxx1Qxy(&d->f, &s->f, &even_odd);
}

static void compute_Sxx1Sxy(vFermion *d,
                            const vFermion *s,
                            struct neighbor *nb);
static void inline compute_See1Seo(vEvenFermion *d, const vOddFermion *s)
{
  compute_Sxx1Sxy(&d->f, &s->f, &even_odd);
}

static void compute_1Qxx1Qxy(vFermion *d,
                             double *norm,
                             const vFermion *q,
                             const vFermion *s,
                             struct neighbor *nb);
static void inline compute_1Qoo1Qoe(vOddFermion *d,
                                    double *norm,
                                    const vOddFermion *q,
                                    const vEvenFermion *s)
{
  compute_1Qxx1Qxy(&d->f, norm, &q->f, &s->f, &odd_even);
}
#line 4119 "dwf.nw"
static inline int
parity(const int x[DIM])
{
    int i, v;
    for (i = v = 0; i < DIM; i++)
        v += x[i];
    return v & 1;
} 
#line 4131 "dwf.nw"
static double
d_pow(double x, unsigned int n)
{
     double v = 1;

     while (n) {
        if (n & 1)
            v *= x;
        x *= x;
        n /= 2;
     }
     return v;
}
#line 4148 "dwf.nw"
static inline void 
vhfzero(vHalfFermion *v)
{
   vReal z = vmk_1(0.0);

   v->f[0][0].re = v->f[0][0].im = 
   v->f[0][1].re = v->f[0][1].im = 
   v->f[0][2].re = v->f[0][2].im = 
   v->f[1][0].re = v->f[1][0].im = 
   v->f[1][1].re = v->f[1][1].im = 
   v->f[1][2].re = v->f[1][2].im = z;
}
#line 1075 "dwf.nw"
static void *
alloc16(int size)
{
    int xsize = PAD16(size + sizeof (struct memblock));
    struct memblock *p = tmalloc(xsize);
	
    if (p == 0)
        return p;
	
    p->data = ALIGN16(&p[1]);
    p->size = size;
    p->next = memblock.next;
    p->prev = &memblock;
    p->next->prev = p;
    p->prev->next = p;
	
    return p->data;
}
#line 1096 "dwf.nw"
static void
free16(void *ptr)
{
    struct memblock *p;

    if (ptr == 0)
        return;
		
    for (p = memblock.next; p != &memblock; p = p->next) {
        if (p->data != ptr)
            continue;
        p->next->prev = p->prev;
        p->prev->next = p->next;
        tfree(p);
        return;
    }
    /* this is BAD: control should not reach here! */
}
#line 1144 "dwf.nw"
vEvenFermion *
allocate_even_fermion(void)
{
    return alloc16(even_odd.size * Sv * sizeof (vFermion));
}

vOddFermion *
allocate_odd_fermion(void)
{
    return alloc16(odd_even.size * Sv * sizeof (vFermion));
}

L3(DWF_Gauge) *
allocate_gauge_field(void)
{
    return alloc16(gauge_XYZT * sizeof (L3(DWF_Gauge)));
}
#line 1397 "dwf.nw"
static int
init_tables(void)
{
    struct neighbor tmp;
    int i, v;

    init_neighbor(&bounds, &neighbor);
    
#line 1419 "dwf.nw"
Sv = tlattice[DIM] / Vs;
Sv_1 = Sv - 1;
for (v = 1, i = 0; i < DIM; i++) {
    v *= bounds.hi[i] - bounds.lo[i];
}
gauge_XYZT = DIM * v;
for (i = 0; i < DIM; i++) {
    if (network[i] < 2)
        continue;
    gauge_XYZT += v / (bounds.hi[i] - bounds.lo[i]);
}
#line 1405 "dwf.nw"
    tmp = neighbor;
    build_neighbor(&even_odd, 0, &tmp);
    build_neighbor(&odd_even, 1, &tmp);

    return 0;
}
#line 1483 "dwf.nw"
static void
init_neighbor(struct bounds *bounds, struct neighbor *neighbor)
{
    int i;

    mk_sublattice(bounds, coord);
#ifndef NO_DEBUG_QMP
    for (i = 0; i < DIM; i++)
       DEBUG_QMP("local: bounds[%d]: lo %d, hi %d\n",
                 i, bounds->lo[i], bounds->hi[i])
#endif /* defined(DEBUG_QMP) */
    neighbor->qmp_smask = 0;
    
#line 1502 "dwf.nw"
for (neighbor->size = 1, neighbor->inside_size = 1, i = 0; i < DIM; i++) {
    int ext = bounds->hi[i] - bounds->lo[i];

    neighbor->size *= ext;
    if (network[i] > 1)
       neighbor->inside_size *= ext - 2;
    else
       neighbor->inside_size *= ext;
}
neighbor->boundary_size = neighbor->size - neighbor->inside_size;
neighbor->site = tmalloc(neighbor->size * sizeof (struct site));
#ifdef DEBUG_QMP
memset(neighbor->site, -1, neighbor->size * sizeof (struct site));
#endif
#line 1496 "dwf.nw"
    
#line 1518 "dwf.nw"
if (neighbor->inside_size)
    neighbor->inside = tmalloc(neighbor->inside_size * sizeof (int));
else
    neighbor->inside = 0;
#line 1497 "dwf.nw"
    
#line 1524 "dwf.nw"
if (neighbor->boundary_size)
    neighbor->boundary = tmalloc(neighbor->boundary_size * sizeof (struct boundary));
else
    neighbor->boundary = 0;
#line 1498 "dwf.nw"
    
#line 1530 "dwf.nw"
for (i = 0; i < 2 * DIM; i++) {
    int d = i / 2;

    if (network[d] > 1) {
        neighbor->snd_size[i] = neighbor->size / (bounds->hi[d] - bounds->lo[d]);
        neighbor->snd[i] = tmalloc(neighbor->snd_size[i] * sizeof (int));
#ifdef DEBUG_QMP
	memset(neighbor->snd[i], -1, neighbor->snd_size[i] * sizeof (int));
#endif
    } else {
        neighbor->snd_size[i] = 0;
        neighbor->snd[i] = 0;
    }
    DEBUG_QMP("Compute send sizes... snd_size[%d]=%d\n",
              i, neighbor->snd_size[i])
}
#line 1499 "dwf.nw"
}
#line 1594 "dwf.nw"
static void
build_neighbor(struct neighbor *out,
	       int              par,
	       struct neighbor *in)
{
   int i,d, s, p, m;
   int x[DIM];

   
#line 1615 "dwf.nw"
*out = *in;
out->size = 0;
out->inside_size = 0;
out->boundary_size = 0;
for (d = 0; d < DIM; d++) {
  out->rcv_size[2*d] = out->snd_size[2*d] = 0;
  out->rcv_size[2*d+1] = out->snd_size[2*d+1] = 0;
}
#line 1603 "dwf.nw"
   
#line 1235 "dwf.nw"
for (i = 0; i < DIM; i++)
    x[i] = bounds.lo[i];
for (i = 0; i < DIM;) {
#line 1627 "dwf.nw"
    
#line 1646 "dwf.nw"
s = parity(x);
if (s != par)
    goto next;
#line 1628 "dwf.nw"
    
#line 1654 "dwf.nw"
p = to_HFlinear(x, &bounds, -1, 0);
for (m = 0, d = 0; d < DIM; d++) {
    if (network[d] > 1) {
        if (x[d] == bounds.lo[d])
            m |= 1 << (2 * d);
        if (x[d] + 1 == bounds.hi[d])
            m |= 1 << (2 * d + 1);
    }
}
#line 1629 "dwf.nw"
    DEBUG_QMP("A: x[%d %d %d %d], (s,par)=(%d,%d), p=%d\n",
               x[0], x[1], x[2], x[3], s, par, p)
    
#line 1667 "dwf.nw"
if (m) {
    
#line 1683 "dwf.nw"
in->boundary->index = p;
in->boundary->mask = m;
in->boundary++;
out->boundary_size++;
#line 1669 "dwf.nw"
} else {
    
#line 1676 "dwf.nw"
*in->inside++ = p;
out->inside_size++;
#line 1671 "dwf.nw"
}
#line 1632 "dwf.nw"
    
#line 1691 "dwf.nw"
in->site->Uup = to_Ulinear(x, &bounds, -1);
for (d = 0; d < DIM; d++) {
    in->site->Udown[d] = to_Ulinear(x, &bounds, d);
    if ((m & (1 << (2 * d))) == 0)
        in->site->F[2*d] = Sv * to_HFlinear(x, &bounds, d, -1);
    if ((m & (1 << (2 * d + 1))) == 0)
        in->site->F[2*d + 1] = Sv * to_HFlinear(x, &bounds, d, +1);
}
#line 1633 "dwf.nw"
    out->size++;
    in->site++;
  next:
#line 1241 "dwf.nw"
    for (i = 0; i < DIM; i++) {
        
#line 1258 "dwf.nw"
if (++x[i] == bounds.hi[i])
    x[i] = bounds.lo[i];
else
    break;
#line 1243 "dwf.nw"
    }
}
#line 1604 "dwf.nw"
   
#line 1704 "dwf.nw"
for (d = 0; d < DIM; d++) {
    if (network[d] < 2)
        continue;
    construct_rec(out, par, &bounds, d, +1);
    construct_snd(out, in, par, &bounds, d, +1);
    construct_rec(out, par, &bounds, d, -1);
    construct_snd(out, in, par, &bounds, d, -1);
}
#line 1605 "dwf.nw"
}
#line 1723 "dwf.nw"
static void
construct_rec(struct neighbor *out,
              int par,
              struct bounds *bounds,
              int dir,
              int step)
{
     struct bounds xb;
     int xc[DIM], x[DIM];
     int s, d, p, k;
     int dz = dir * 2 + ((step>0)?1:0);

     
#line 1749 "dwf.nw"
for (d = 0; d < DIM; d++) {
    int v = coord[d] + ((d==dir)?step:0);

    if (v < 0)
        v += network[d];
    if (v >= network[d])
        v -= network[d];
    xc[d] = v;
}
mk_sublattice(&xb, xc);
#ifndef NO_DEBUG_QMP
   DEBUG_QMP("par=%d, dir=%d, step=%d\n", par, dir, step)
   for (d = 0; d < DIM; d++)
      DEBUG_QMP("neighbor: xb[%d] lo %d, di %d\n", d, xb.lo[d], xb.hi[d])
#endif /* defined(DEBUG_QMP) */
#line 1736 "dwf.nw"
     
#line 1767 "dwf.nw"
for (d = 0; d < DIM; d++)
    x[d] = ((d == dir) && (step < 0))? (xb.hi[d] - 1): xb.lo[d];
#line 1737 "dwf.nw"
     
#line 1815 "dwf.nw"
for (k = 0, d = 0; d < DIM; ) {
#line 1738 "dwf.nw"
         
#line 1641 "dwf.nw"
s = parity(x);
if (s == par)
    goto next;
#line 1739 "dwf.nw"
         
#line 1771 "dwf.nw"
p = to_HFlinear(x, bounds, dir, -step);
#line 1740 "dwf.nw"
         DEBUG_QMP("B: x[%d %d %d %d], (s,par)=(%d,%d), p=%d, k=%d\n",
                   x[0], x[1], x[2], x[3], s, par, p, k)
         
#line 1774 "dwf.nw"
out->site[p].F[dz] = Sv * k++;
#line 1743 "dwf.nw"
     
#line 1818 "dwf.nw"
  next:
    for (d = 0; d < DIM; d++) {
        if (d == dir)
            continue;
        if (++x[d] == xb.hi[d])
            x[d] = xb.lo[d];
        else
            break;
    }
}
#line 1744 "dwf.nw"
     out->rcv_size[dz^1] = k;
}
#line 1788 "dwf.nw"
static void
construct_snd(struct neighbor *out,
              struct neighbor *in,
              int par,
              struct bounds *bounds,
              int dir,
              int step)
{
     struct bounds xb;
     int xc[DIM], x[DIM];
     int s, d, p, k;
     int dz = dir * 2 + ((step>0)?1:0);

     
#line 1749 "dwf.nw"
for (d = 0; d < DIM; d++) {
    int v = coord[d] + ((d==dir)?step:0);

    if (v < 0)
        v += network[d];
    if (v >= network[d])
        v -= network[d];
    xc[d] = v;
}
mk_sublattice(&xb, xc);
#ifndef NO_DEBUG_QMP
   DEBUG_QMP("par=%d, dir=%d, step=%d\n", par, dir, step)
   for (d = 0; d < DIM; d++)
      DEBUG_QMP("neighbor: xb[%d] lo %d, di %d\n", d, xb.lo[d], xb.hi[d])
#endif /* defined(DEBUG_QMP) */
#line 1802 "dwf.nw"
     
#line 1767 "dwf.nw"
for (d = 0; d < DIM; d++)
    x[d] = ((d == dir) && (step < 0))? (xb.hi[d] - 1): xb.lo[d];
#line 1803 "dwf.nw"
     
#line 1815 "dwf.nw"
for (k = 0, d = 0; d < DIM; ) {
#line 1804 "dwf.nw"
         
#line 1646 "dwf.nw"
s = parity(x);
if (s != par)
    goto next;
#line 1805 "dwf.nw"
         
#line 1771 "dwf.nw"
p = to_HFlinear(x, bounds, dir, -step);
#line 1806 "dwf.nw"
         DEBUG_QMP("C: x[%d %d %d %d], (s,par)=(%d,%d), p=%d, k=%d\n",
                   x[0], x[1], x[2], x[3], s, par, p, k)
         *in->snd[dz]++ = p * Sv;
         k++;
     
#line 1818 "dwf.nw"
  next:
    for (d = 0; d < DIM; d++) {
        if (d == dir)
            continue;
        if (++x[d] == xb.hi[d])
            x[d] = xb.lo[d];
        else
            break;
    }
}
#line 1811 "dwf.nw"
     out->snd_size[dz] = k;
}
#line 1952 "dwf.nw"
static int
build_buffers(struct neighbor *nb)
{
    int i, k, Nh;
    QMP_msghandle_t SRh[4*DIM];

    DEBUG_QMP("--------------------------------------------\n")
    DEBUG_QMP("build buffers [%s]\n", nb == &even_odd? "even": "odd")
#ifndef NO_DEBUG_QMP
    {
      int i;

      DEBUG_QMP("nb->size        = %d\n", nb->size)
      DEBUG_QMP("nb->inside_size = %d\n", nb->inside_size)
      DEBUG_QMP("nb->boundary_size = %d\n", nb->boundary_size)
      for (i = 0; i < 2 * DIM; i++)
         DEBUG_QMP("[%d]: snd=%d, rcv=%d\n",
	            i, nb->snd_size[i], nb->rcv_size[i])
    }
#endif /* defined(DEBUG_QMP) */
    Nh = nb->Nx = 0;
    for (i = 0; i < DIM; i++) {
        switch (network[i]) {
            case 1: break;
            case 2:
                 
#line 1993 "dwf.nw"
DEBUG_QMP("Allocate up and down buffers, i=%d\n", i)
k = make_buffer(nb, nb->snd_size[2*i] + nb->snd_size[2*i+1]);
nb->snd_buf[2*i] = nb->qmp_buf[k];
nb->snd_buf[2*i+1] = nb->snd_buf[2*i] + Sv * nb->snd_size[2*i];
Nh = make_send(nb, k, i, +1, SRh, Nh);

k = make_buffer(nb, nb->rcv_size[2*i] + nb->rcv_size[2*i+1]);
nb->rcv_buf[2*i+1] = nb->qmp_buf[k]; /* should be opposite to snd_buf[] */
nb->rcv_buf[2*i] = nb->rcv_buf[2*i+1] + Sv * nb->rcv_size[2*i+1];
Nh = make_receive(nb, k, i, -1, SRh, Nh); /* -1 fixes a bug in GigE QMP */
#line 1978 "dwf.nw"
                 break;
            default:
	         /* Order here is important */
                 
#line 2016 "dwf.nw"
DEBUG_QMP("Allocate down buffers, i=%d\n", i)
k = make_buffer(nb, nb->snd_size[2*i]);
nb->snd_buf[2*i] = nb->qmp_buf[k];
Nh = make_send(nb, k, i, -1, SRh, Nh);

k = make_buffer(nb, nb->rcv_size[2*i]);
nb->rcv_buf[2*i] = nb->qmp_buf[k];
Nh = make_receive(nb, k, i, -1, SRh, Nh);
#line 1982 "dwf.nw"
                 
#line 2006 "dwf.nw"
DEBUG_QMP("Allocate up buffers, i=%d\n", i)
k = make_buffer(nb, nb->snd_size[2*i+1]);
nb->snd_buf[2*i+1] = nb->qmp_buf[k];
Nh = make_send(nb, k, i, +1, SRh, Nh);

k = make_buffer(nb, nb->rcv_size[2*i+1]);
nb->rcv_buf[2*i+1] = nb->qmp_buf[k];
Nh = make_receive(nb, k, i, +1, SRh, Nh);
#line 1983 "dwf.nw"
                 break;
        }
    }
    
#line 2088 "dwf.nw"
if (nb->qmp_smask) {
    nb->qmp_handle = QMP_declare_multiple(SRh, Nh);
#ifndef NO_DEBUG_QMP
    {
        int i;
        for (i = 0; i < Nh; i++) {
            DEBUG_QMP("declare_multiple([%d]=0x%x)\n", i, (int)SRh[i])
        }
        DEBUG_QMP("declare_multiple(..., %d)=0x%x\n", Nh, (int)nb->qmp_handle)
    }
#endif
}
#line 1987 "dwf.nw"
    return 0;
}
#line 2030 "dwf.nw"
static int
make_buffer(struct neighbor *nb, int size)
{
    int bcount = size * Sv * sizeof (vHalfFermion);
    int N = nb->Nx;

    nb->qmp_size[N] = size;
    sse_aligned_buffer(nb, N, bcount);
    nb->qmp_mm[N] = QMP_declare_msgmem(nb->qmp_buf[N], bcount);
    nb->Nx = N + 1;
    DEBUG_QMP("declare_msgmem(%p,%d)=0x%x\n",
              nb->qmp_buf[N], bcount, (int)nb->qmp_mm[N])
    return N;  
}
#line 2052 "dwf.nw"
static int
make_send(struct neighbor *nb, int k, int i, int d,
          QMP_msghandle_t SRh[4*DIM], int Nsr)
{
    int j = 2 * i + ((d < 0)? 0: 1);

    nb->qmp_smask |= (1 << j);
    SRh[Nsr] = QMP_declare_send_relative(nb->qmp_mm[k], i, d, 1);

    DEBUG_QMP("declare_send_relative(0x%x,%d,%d,1)=0x%x\n",
              (int)nb->qmp_mm[k], i, d, (int)SRh[Nsr])

    return Nsr+1;
}
#line 2074 "dwf.nw"
static int
make_receive(struct neighbor *nb, int k, int i, int d,
             QMP_msghandle_t SRh[4*DIM], int Nsr)
{
    SRh[Nsr] = QMP_declare_receive_relative(nb->qmp_mm[k], i, d, 1);

    DEBUG_QMP("declare_receive_relative(0x%x,%d,%d,1)=0x%x\n",
              (int)nb->qmp_mm[k], i, d, (int)SRh[Nsr])

    return Nsr+1;
}
#line 2110 "dwf.nw"
static void
sse_aligned_buffer(struct neighbor *nb, int k, int size)
{
#ifdef USE_QMP2
  nb->qmp_xbuf[k] = QMP_allocate_aligned_memory(size, 128, 0);
  nb->qmp_buf[k] = QMP_get_memory_pointer(nb->qmp_xbuf[k]);
#else
  int xcount = size + 15;
  char *ptr = QMP_allocate_aligned_memory(xcount);
  nb->qmp_buf[k] = (void *)(~15 & (15 + (unsigned long)(ptr)));
  nb->qmp_xbuf[k] = ptr;
#endif
  DEBUG_QMP("(%p,%d,%d): allocate: 0x%x\n",
	    nb, k, size, (int)nb->qmp_xbuf[k])
  DEBUG_QMP("(%p,%d,%d): ptr: %p\n",
	    nb, k, size, (void *)nb->qmp_buf[k])
}
#line 2140 "dwf.nw"
static void
free_buffers(struct neighbor *nb)
{
   int i;

   
#line 2155 "dwf.nw"
if (nb->qmp_handle) {
   QMP_free_msghandle(nb->qmp_handle);
   DEBUG_QMP("free_msghandle(0x%x) / common receive handle\n",
             (int)nb->qmp_handle)
}
#line 2146 "dwf.nw"
   
#line 2163 "dwf.nw"
for (i = nb->Nx; i--;) {
       if (nb->qmp_mm[i])
          QMP_free_msgmem(nb->qmp_mm[i]);
       DEBUG_QMP("free_msgmem(0x%x)\n", (int)nb->qmp_mm[i]);       
#ifdef USE_QMP2
       if (nb->qmp_xbuf[i])
          QMP_free_memory(nb->qmp_xbuf[i]);
#else /* QMP 1.3 */
       if (nb->qmp_xbug[i])
          QMP_free_aligned_memory(nb->qmp_xbuf[i]);
#endif
       DEBUG_QMP("free_memory(0x%x)\n", (int)nb->qmp_xbuf[i]);       
}
#line 2147 "dwf.nw"
}
#line 2224 "dwf.nw"
static int
cg(vOddFermion *x_o,
   const vOddFermion *b,
   const vOddFermion *x0,
   double epsilon,
   int N0,
   int N,
   double *out_eps,
   int *out_N)
{
    double rho, alpha, beta, gamma, norm_z;
    int status = 1;
    int k;

    copy_o(x_o, x0);
    compute_MxM(p_o, &norm_z, x_o);
    compute_sum_oN(r_o, &rho, b, -1, p_o);
    copy_o(p_o, r_o);
    
#line 4067 "dwf.nw"
/* relax, QMP does not support split reductions yet. */

#line 2244 "dwf.nw"
    for (k = 0; (k < N0) || ((rho > epsilon) && (k < N)); k++) {
        compute_MxM(q_o, &norm_z, p_o);
        
#line 4067 "dwf.nw"
/* relax, QMP does not support split reductions yet. */
#line 2247 "dwf.nw"
        alpha = rho / norm_z;
        compute_sum2_oN(r_o, &gamma, -alpha, q_o);
        compute_sum2_o(x_o, alpha, p_o);
        
#line 4067 "dwf.nw"
/* relax, QMP does not support split reductions yet. */
#line 2251 "dwf.nw"
        DEBUG_DWF("cg loop: k=%d, rho=%g, norm_z=%g, alpha=%g, gamma=%g\n",
                   k, rho, norm_z, alpha, gamma)
        if (gamma <= epsilon) {
            rho = gamma;
            status = 0;
            break;
        }
	beta = gamma / rho;
	rho = gamma;
	compute_sum2x_o(p_o, r_o, beta);
    }
    *out_N = k;
    *out_eps = rho;

    return status;
}
#line 2307 "dwf.nw"
static void
copy_o(vOddFermion *dst, const vOddFermion *src)
{
  int i = odd_even.size * Sv * sizeof (vOddFermion) / sizeof (vReal);
  vReal *d = (vReal *)dst;
  const vReal *s = (const vReal *)src;

  for ( ;i--;)
      *d++ = *s++;
}
#line 2326 "dwf.nw"
static void
compute_sum2_o(vOddFermion *dst, double alpha, const vOddFermion *src)
{
  int i = odd_even.size * Sv * sizeof (vOddFermion) / sizeof (vReal);
  vReal a = vmk_1(alpha);
  vReal *d = (vReal *)dst;
  const vReal *s = (const vReal *)src;

  for ( ;i--;)
      *d++ += a * *s++;
}
#line 2344 "dwf.nw"
static void
compute_sum2x_o(vOddFermion *dst, const vOddFermion *src, double alpha)
{
  int i = odd_even.size * Sv * sizeof (vOddFermion) / sizeof (vReal);
  vReal a = vmk_1(alpha);
  vReal *d = (vReal *)dst;
  const vReal *s = (const vReal *)src;

  for ( ;i--; d++)
      *d = a * *d + *s++;
}
#line 2366 "dwf.nw"
static void
compute_sum_e(vEvenFermion *d,
              const vEvenFermion *x, double alpha, const vEvenFermion *y)
{
  const vReal *X = (const vReal *)x;
  const vReal *Y = (const vReal *)y;
  vReal *D = (vReal *)d;
  vReal a = vmk_1(alpha);
  int i = even_odd.size * Sv * sizeof (vEvenFermion) / sizeof (vReal);

  for (;i--;)
     *D++ = *X++ + a * *Y++;
}
#line 2381 "dwf.nw"
static void
compute_sum_o(vOddFermion *d,
              const vOddFermion *x, double alpha, const vOddFermion *y)
{
  const vReal *X = (const vReal *)x;
  const vReal *Y = (const vReal *)y;
  vReal *D = (vReal *)d;
  vReal a = vmk_1(alpha);
  int i = odd_even.size * Sv * sizeof (vOddFermion) / sizeof (vReal);

  for (;i--;)
     *D++ = *X++ + a * *Y++;
}
#line 2405 "dwf.nw"
static void
compute_sum_oN(vOddFermion *d, double *norm,
               const vOddFermion *x, double alpha, const vOddFermion *y)
{
  vFermion *D = &d->f;
  const vFermion *X = &x->f;
  const vFermion *Y = &y->f;
  int n = odd_even.size * Sv;
  vReal a = vmk_1(alpha);
  vReal q0r, q1r, q2r, q0i, q1i, q2i;
  vReal s0r, s1r, s2r, s0i, s1i, s2i;
  int i;

  *norm = 0;
  for (i = 0; i < n; i++, X++, Y++, D++) {
    q0r = D->f[0][0].re = X->f[0][0].re + a * Y->f[0][0].re; s0r = q0r * q0r;
    q0i = D->f[0][0].im = X->f[0][0].im + a * Y->f[0][0].im; s0i = q0i * q0i;
    q1r = D->f[0][1].re = X->f[0][1].re + a * Y->f[0][1].re; s1r = q1r * q1r;
    q1i = D->f[0][1].im = X->f[0][1].im + a * Y->f[0][1].im; s1i = q1i * q1i;
    q2r = D->f[0][2].re = X->f[0][2].re + a * Y->f[0][2].re; s2r = q2r * q2r;
    q2i = D->f[0][2].im = X->f[0][2].im + a * Y->f[0][2].im; s2i = q2i * q2i;

    q0r = D->f[1][0].re = X->f[1][0].re + a * Y->f[1][0].re; s0r += q0r * q0r;
    q0i = D->f[1][0].im = X->f[1][0].im + a * Y->f[1][0].im; s0i += q0i * q0i;
    q1r = D->f[1][1].re = X->f[1][1].re + a * Y->f[1][1].re; s1r += q1r * q1r;
    q1i = D->f[1][1].im = X->f[1][1].im + a * Y->f[1][1].im; s1i += q1i * q1i;
    q2r = D->f[1][2].re = X->f[1][2].re + a * Y->f[1][2].re; s2r += q2r * q2r;
    q2i = D->f[1][2].im = X->f[1][2].im + a * Y->f[1][2].im; s2i += q2i * q2i;

    q0r = D->f[2][0].re = X->f[2][0].re + a * Y->f[2][0].re; s0r += q0r * q0r;
    q0i = D->f[2][0].im = X->f[2][0].im + a * Y->f[2][0].im; s0i += q0i * q0i;
    q1r = D->f[2][1].re = X->f[2][1].re + a * Y->f[2][1].re; s1r += q1r * q1r;
    q1i = D->f[2][1].im = X->f[2][1].im + a * Y->f[2][1].im; s1i += q1i * q1i;
    q2r = D->f[2][2].re = X->f[2][2].re + a * Y->f[2][2].re; s2r += q2r * q2r;
    q2i = D->f[2][2].im = X->f[2][2].im + a * Y->f[2][2].im; s2i += q2i * q2i;

    q0r = D->f[3][0].re = X->f[3][0].re + a * Y->f[3][0].re; s0r += q0r * q0r;
    q0i = D->f[3][0].im = X->f[3][0].im + a * Y->f[3][0].im; s0i += q0i * q0i;
    q1r = D->f[3][1].re = X->f[3][1].re + a * Y->f[3][1].re; s1r += q1r * q1r;
    q1i = D->f[3][1].im = X->f[3][1].im + a * Y->f[3][1].im; s1i += q1i * q1i;
    q2r = D->f[3][2].re = X->f[3][2].re + a * Y->f[3][2].re; s2r += q2r * q2r;
    q2i = D->f[3][2].im = X->f[3][2].im + a * Y->f[3][2].im; s2i += q2i * q2i;

    *norm += vsum(s0r + s0i + s1r + s1i + s2r + s2i);
    }
  
#line 4112 "dwf.nw"
QMP_sum_double(norm);
DEBUG_QMP("sum_double(%p)\n", norm)
#line 2451 "dwf.nw"
}
#line 2459 "dwf.nw"
static void
compute_sum2_oN(vOddFermion *d, double *norm,
                double alpha, const vOddFermion *y)
{
  vFermion *D = &d->f;
  const vFermion *Y = &y->f;
  int n = odd_even.size * Sv;
  vReal a = vmk_1(alpha);
  vReal q0r, q1r, q2r, q0i, q1i, q2i;
  vReal s0r, s1r, s2r, s0i, s1i, s2i;
  int i;

  *norm = 0;
  for (i = 0; i < n; i++, Y++, D++) {
    q0r = D->f[0][0].re += a * Y->f[0][0].re; s0r = q0r * q0r;
    q0i = D->f[0][0].im += a * Y->f[0][0].im; s0i = q0i * q0i;
    q1r = D->f[0][1].re += a * Y->f[0][1].re; s1r = q1r * q1r;
    q1i = D->f[0][1].im += a * Y->f[0][1].im; s1i = q1i * q1i;
    q2r = D->f[0][2].re += a * Y->f[0][2].re; s2r = q2r * q2r;
    q2i = D->f[0][2].im += a * Y->f[0][2].im; s2i = q2i * q2i;

    q0r = D->f[1][0].re += a * Y->f[1][0].re; s0r += q0r * q0r;
    q0i = D->f[1][0].im += a * Y->f[1][0].im; s0i += q0i * q0i;
    q1r = D->f[1][1].re += a * Y->f[1][1].re; s1r += q1r * q1r;
    q1i = D->f[1][1].im += a * Y->f[1][1].im; s1i += q1i * q1i;
    q2r = D->f[1][2].re += a * Y->f[1][2].re; s2r += q2r * q2r;
    q2i = D->f[1][2].im += a * Y->f[1][2].im; s2i += q2i * q2i;

    q0r = D->f[2][0].re += a * Y->f[2][0].re; s0r += q0r * q0r;
    q0i = D->f[2][0].im += a * Y->f[2][0].im; s0i += q0i * q0i;
    q1r = D->f[2][1].re += a * Y->f[2][1].re; s1r += q1r * q1r;
    q1i = D->f[2][1].im += a * Y->f[2][1].im; s1i += q1i * q1i;
    q2r = D->f[2][2].re += a * Y->f[2][2].re; s2r += q2r * q2r;
    q2i = D->f[2][2].im += a * Y->f[2][2].im; s2i += q2i * q2i;

    q0r = D->f[3][0].re += a * Y->f[3][0].re; s0r += q0r * q0r;
    q0i = D->f[3][0].im += a * Y->f[3][0].im; s0i += q0i * q0i;
    q1r = D->f[3][1].re += a * Y->f[3][1].re; s1r += q1r * q1r;
    q1i = D->f[3][1].im += a * Y->f[3][1].im; s1i += q1i * q1i;
    q2r = D->f[3][2].re += a * Y->f[3][2].re; s2r += q2r * q2r;
    q2i = D->f[3][2].im += a * Y->f[3][2].im; s2i += q2i * q2i;

    *norm += vsum(s0r + s0i + s1r + s1i + s2r + s2i);
    }
  
#line 4112 "dwf.nw"
QMP_sum_double(norm);
DEBUG_QMP("sum_double(%p)\n", norm)
#line 2504 "dwf.nw"
}
#line 2519 "dwf.nw"
static void
compute_MxM(vOddFermion *eta, double *norm,
            const vOddFermion *psi)
{
     compute_M(auxB_o, norm, psi);
     compute_Mx(eta, auxB_o);
}
#line 2530 "dwf.nw"
static void compute_M(vOddFermion *eta, double *norm,
                      const vOddFermion *psi)
{
   compute_Qee1Qeo(auxA_e, psi);
   compute_1Qoo1Qoe(eta, norm, psi, auxA_e);
}
#line 2540 "dwf.nw"
static void compute_Mx(vOddFermion *eta,
                       const vOddFermion *psi)
{
   compute_Soo1(auxA_o, psi);
   compute_See1Seo(auxA_e, auxA_o);
   compute_1Soe(eta, psi, auxA_e);
}
#line 2567 "dwf.nw"
static void
compute_Qxx1(vFermion *chi, const vFermion *psi, int size)
{
    const vFermion *qs, *qx5;
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 2572 "dwf.nw"
    
#line 2767 "dwf.nw"
vReal fx;
vHalfFermion zV;
vector_complex zn;
scalar_complex zX[Fd/2][Nc];
#line 2955 "dwf.nw"
scalar_complex yOut[Fd/2][Nc];

#line 2574 "dwf.nw"
    for (i = 0; i < size; i++) {
        xyzt5 = i * Sv;
        
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 2577 "dwf.nw"
	
#line 4045 "dwf.nw"
qx5 = &psi[xyzt5];
#line 2578 "dwf.nw"
        
#line 2774 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2777 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2781 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2786 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2789 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2922 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2969 "dwf.nw"
  COMPUTE_YA(0,0,0) COMPUTE_YA(0,0,1) COMPUTE_YA(0,0,2)
#line 2925 "dwf.nw"
    
#line 2973 "dwf.nw"
  COMPUTE_YA(1,1,0) COMPUTE_YA(1,1,1) COMPUTE_YA(1,1,2)
#line 2926 "dwf.nw"
}
#line 2881 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2884 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
  rs = &rx5[s];
  QSETUP(s)
  
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2888 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2893 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2895 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2946 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2993 "dwf.nw"
  COMPUTE_YB(2,0,0) COMPUTE_YB(2,0,1) COMPUTE_YB(2,0,2)
#line 2949 "dwf.nw"
    
#line 2997 "dwf.nw"
  COMPUTE_YB(3,1,0) COMPUTE_YB(3,1,1) COMPUTE_YB(3,1,2)
#line 2950 "dwf.nw"
}
#line 2579 "dwf.nw"
    }
}
#line 2587 "dwf.nw"
static void
compute_Soo1(vOddFermion *Chi, const vOddFermion *Psi)
{
    vFermion *chi = &Chi->f;
    const vFermion *psi = &Psi->f;
    int size = odd_even.size;
    const vFermion *qs, *qx5;
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 2595 "dwf.nw"
    
#line 2767 "dwf.nw"
vReal fx;
vHalfFermion zV;
vector_complex zn;
scalar_complex zX[Fd/2][Nc];
#line 2955 "dwf.nw"
scalar_complex yOut[Fd/2][Nc];

#line 2597 "dwf.nw"
    for (i = 0; i < size; i++) {
        xyzt5 = i * Sv;
        
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 2600 "dwf.nw"
	
#line 4045 "dwf.nw"
qx5 = &psi[xyzt5];
#line 2601 "dwf.nw"
        
#line 2815 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2818 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2822 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2827 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2829 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2930 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2977 "dwf.nw"
  COMPUTE_YA(2,0,0) COMPUTE_YA(2,0,1) COMPUTE_YA(2,0,2)
#line 2933 "dwf.nw"
    
#line 2981 "dwf.nw"
  COMPUTE_YA(3,1,0) COMPUTE_YA(3,1,1) COMPUTE_YA(3,1,2)
#line 2934 "dwf.nw"
}
#line 2853 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2856 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2860 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2865 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2868 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2938 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2985 "dwf.nw"
  COMPUTE_YB(0,0,0) COMPUTE_YB(0,0,1) COMPUTE_YB(0,0,2)
#line 2941 "dwf.nw"
    
#line 2989 "dwf.nw"
  COMPUTE_YB(1,1,0) COMPUTE_YB(1,1,1) COMPUTE_YB(1,1,2)
#line 2942 "dwf.nw"
}
#line 2602 "dwf.nw"
    }
}
#line 3097 "dwf.nw"
static void
compute_Qxy(vFermion *chi,
            const vFermion *psi,
            struct neighbor *nb)
{
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 3103 "dwf.nw"
    
#line 4022 "dwf.nw"
int xyzt, k, d;
const vFermion *f;
vHalfFermion *g;
vHalfFermion gg[2*DIM], hh[2*DIM];
vSU3 V[2*DIM];
int ps[2*DIM], p5[2*DIM];
#line 4030 "dwf.nw"
const SU3 *Uup, *Udown;
int c1, c2;

#line 3105 "dwf.nw"
    
#line 3259 "dwf.nw"
#define qx5 rx5
#define qs rs
#line 3106 "dwf.nw"
    
#line 3326 "dwf.nw"
{
   int k, i, s, c, *src;
   const vFermion *f;
   vHalfFermion *g;

   k = 0; 
#line 3359 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3363 "dwf.nw"
        }
    }
}
#line 3332 "dwf.nw"
   k = 1; 
#line 3368 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3372 "dwf.nw"
        }
    }
}
#line 3333 "dwf.nw"
   k = 2; 
#line 3377 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3381 "dwf.nw"
        }
    }
}
#line 3334 "dwf.nw"
   k = 3; 
#line 3386 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3390 "dwf.nw"
        }
    }
}
#line 3335 "dwf.nw"
   k = 4; 
#line 3395 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3399 "dwf.nw"
        }
    }
}
#line 3336 "dwf.nw"
   k = 5; 
#line 3404 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
     for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3408 "dwf.nw"
        }
    }
}
#line 3337 "dwf.nw"
   k = 6; 
#line 3413 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3417 "dwf.nw"
        }
    }
}
#line 3338 "dwf.nw"
   k = 7; 
#line 3422 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3426 "dwf.nw"
        }
    }
}
#line 3339 "dwf.nw"
}
#line 3107 "dwf.nw"
    
#line 4079 "dwf.nw"
if (nb->qmp_smask) {
    QMP_start(nb->qmp_handle);
    DEBUG_QMP("start sends and receives (0x%x)\n", (int)nb->qmp_handle)
}
#line 3108 "dwf.nw"
    
#line 3433 "dwf.nw"
for (i = 0; i < nb->inside_size; i++) {
    xyzt = nb->inside[i];
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3437 "dwf.nw"
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3438 "dwf.nw"
    
#line 3455 "dwf.nw"
for (s = 0; s < Sv; s++) {
    
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3472 "dwf.nw"
for (c = 0; c < Nc; c++) {
    k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3474 "dwf.nw"
    k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3475 "dwf.nw"
    k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3476 "dwf.nw"
    k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3477 "dwf.nw"
    k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3478 "dwf.nw"
    k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3479 "dwf.nw"
    k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3480 "dwf.nw"
    k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3481 "dwf.nw"
}
#line 3457 "dwf.nw"
    
#line 3528 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3533 "dwf.nw"
}
#line 3458 "dwf.nw"
    
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3459 "dwf.nw"
}
#line 3439 "dwf.nw"
}
#line 3109 "dwf.nw"
    
#line 4088 "dwf.nw"
if (nb->qmp_smask) {
    QMP_wait(nb->qmp_handle);
    DEBUG_QMP("waiting for sends and receives (0x%x)\n",
               (int)nb->qmp_handle)
}
#line 3110 "dwf.nw"
    
#line 3443 "dwf.nw"
for (i = 0; i < nb->boundary_size; i++) {
    int m = nb->boundary[i].mask;

    xyzt = nb->boundary[i].index;
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3449 "dwf.nw"
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3450 "dwf.nw"
    
#line 3463 "dwf.nw"
for (s = 0; s < Sv; s++) {
  
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3486 "dwf.nw"
for (c = 0; c < 3; c++) {
    if ((m & 0x01) == 0) {
        k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3489 "dwf.nw"
    }
    if ((m & 0x02) == 0) {
        k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3492 "dwf.nw"
    }
    if ((m & 0x04) == 0) {
        k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3495 "dwf.nw"
    }
    if ((m & 0x08) == 0) {
        k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3498 "dwf.nw"
    }
    if ((m & 0x10) == 0) {
        k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3501 "dwf.nw"
    }
    if ((m & 0x20) == 0) {
        k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3504 "dwf.nw"
    }
    if ((m & 0x40) == 0) {
        k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3507 "dwf.nw"
    }
    if ((m & 0x80) == 0) {
        k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3510 "dwf.nw"
    }
}
#line 3465 "dwf.nw"
  
#line 3537 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = (m & (1 << d))? &nb->rcv_buf[d][ps[d]]: &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3542 "dwf.nw"
}
#line 3466 "dwf.nw"
  
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3467 "dwf.nw"
}
#line 3451 "dwf.nw"
}
#line 3111 "dwf.nw"
    
#line 3263 "dwf.nw"
#undef qs
#undef qx5
#line 3112 "dwf.nw"
}
#line 3119 "dwf.nw"
static void
compute_1Sxy(vFermion *chi,
             const vFermion *eta,
             const vFermion *psi,
             struct neighbor *nb)
{
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 3126 "dwf.nw"
    
#line 4022 "dwf.nw"
int xyzt, k, d;
const vFermion *f;
vHalfFermion *g;
vHalfFermion gg[2*DIM], hh[2*DIM];
vSU3 V[2*DIM];
int ps[2*DIM], p5[2*DIM];
#line 4030 "dwf.nw"
const SU3 *Uup, *Udown;
int c1, c2;

#line 3128 "dwf.nw"
    
#line 3259 "dwf.nw"
#define qx5 rx5
#define qs rs
#line 3129 "dwf.nw"
    
#line 3343 "dwf.nw"
{
   int k, i, s, c, *src;
   const vFermion *f;
   vHalfFermion *g;

   k = 0; 
#line 3368 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3372 "dwf.nw"
        }
    }
}
#line 3349 "dwf.nw"
   k = 1; 
#line 3359 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3363 "dwf.nw"
        }
    }
}
#line 3350 "dwf.nw"
   k = 2; 
#line 3386 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3390 "dwf.nw"
        }
    }
}
#line 3351 "dwf.nw"
   k = 3; 
#line 3377 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3381 "dwf.nw"
        }
    }
}
#line 3352 "dwf.nw"
   k = 4; 
#line 3404 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
     for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3408 "dwf.nw"
        }
    }
}
#line 3353 "dwf.nw"
   k = 5; 
#line 3395 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3399 "dwf.nw"
        }
    }
}
#line 3354 "dwf.nw"
   k = 6; 
#line 3422 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3426 "dwf.nw"
        }
    }
}
#line 3355 "dwf.nw"
   k = 7; 
#line 3413 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3417 "dwf.nw"
        }
    }
}
#line 3356 "dwf.nw"
}
#line 3130 "dwf.nw"
    
#line 4079 "dwf.nw"
if (nb->qmp_smask) {
    QMP_start(nb->qmp_handle);
    DEBUG_QMP("start sends and receives (0x%x)\n", (int)nb->qmp_handle)
}
#line 3131 "dwf.nw"
    
#line 3563 "dwf.nw"
for (i = 0; i < nb->inside_size; i++) {
    const vFermion *ex5, *es;

    xyzt = nb->inside[i];
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3569 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3571 "dwf.nw"
    
#line 3590 "dwf.nw"
for (s = 0; s < Sv; s++) {
    
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3607 "dwf.nw"
for (c = 0; c < Nc; c++) {
    k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3609 "dwf.nw"
    k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3610 "dwf.nw"
    k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3611 "dwf.nw"
    k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3612 "dwf.nw"
    k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3613 "dwf.nw"
    k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3614 "dwf.nw"
    k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3615 "dwf.nw"
    k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3616 "dwf.nw"
}
#line 3592 "dwf.nw"
    
#line 3528 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3533 "dwf.nw"
}
#line 3593 "dwf.nw"
    
#line 3652 "dwf.nw"
rs = &rx5[s];
es = &ex5[s];
for (c = 0; c < Nc; c++) {
    k = 6; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3656 "dwf.nw"
    k = 7; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3657 "dwf.nw"
    k = 2; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3658 "dwf.nw"
    k = 3; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3659 "dwf.nw"
    k = 1; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3660 "dwf.nw"
    k = 0; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3661 "dwf.nw"
    k = 5; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3662 "dwf.nw"
    k = 4; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3663 "dwf.nw"
    
#line 3668 "dwf.nw"
rs->f[0][c].re = es->f[0][c].re - rs->f[0][c].re;
rs->f[0][c].im = es->f[0][c].im - rs->f[0][c].im;
rs->f[1][c].re = es->f[1][c].re - rs->f[1][c].re;
rs->f[1][c].im = es->f[1][c].im - rs->f[1][c].im;
rs->f[2][c].re = es->f[2][c].re - rs->f[2][c].re;
rs->f[2][c].im = es->f[2][c].im - rs->f[2][c].im;
rs->f[3][c].re = es->f[3][c].re - rs->f[3][c].re;
rs->f[3][c].im = es->f[3][c].im - rs->f[3][c].im;
#line 3664 "dwf.nw"
}
#line 3594 "dwf.nw"
}
#line 3572 "dwf.nw"
}
#line 3132 "dwf.nw"
    
#line 4088 "dwf.nw"
if (nb->qmp_smask) {
    QMP_wait(nb->qmp_handle);
    DEBUG_QMP("waiting for sends and receives (0x%x)\n",
               (int)nb->qmp_handle)
}
#line 3133 "dwf.nw"
    
#line 3576 "dwf.nw"
for (i = 0; i < nb->boundary_size; i++) {
    const vFermion *ex5, *es;
    int m = nb->boundary[i].mask;

    xyzt = nb->boundary[i].index;
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3583 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3585 "dwf.nw"
    
#line 3598 "dwf.nw"
for (s = 0; s < Sv; s++) {
  

#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3623 "dwf.nw"
for (c = 0; c < Nc; c++) {
    if ((m & 0x01) == 0) {
        k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3626 "dwf.nw"
    }
    if ((m & 0x02) == 0) {
        k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3629 "dwf.nw"
    }
    if ((m & 0x04) == 0) {
        k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3632 "dwf.nw"
    }
    if ((m & 0x08) == 0) {
        k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3635 "dwf.nw"
    }
    if ((m & 0x10) == 0) {
        k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3638 "dwf.nw"
    }
    if ((m & 0x20) == 0) {
        k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3641 "dwf.nw"
    }
    if ((m & 0x40) == 0) {
        k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3644 "dwf.nw"
    }
    if ((m & 0x80) == 0) {
        k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3647 "dwf.nw"
    }
}
#line 3600 "dwf.nw"
  
#line 3537 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = (m & (1 << d))? &nb->rcv_buf[d][ps[d]]: &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3542 "dwf.nw"
}
#line 3601 "dwf.nw"
  
#line 3652 "dwf.nw"
rs = &rx5[s];
es = &ex5[s];
for (c = 0; c < Nc; c++) {
    k = 6; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3656 "dwf.nw"
    k = 7; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3657 "dwf.nw"
    k = 2; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3658 "dwf.nw"
    k = 3; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3659 "dwf.nw"
    k = 1; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3660 "dwf.nw"
    k = 0; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3661 "dwf.nw"
    k = 5; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3662 "dwf.nw"
    k = 4; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3663 "dwf.nw"
    
#line 3668 "dwf.nw"
rs->f[0][c].re = es->f[0][c].re - rs->f[0][c].re;
rs->f[0][c].im = es->f[0][c].im - rs->f[0][c].im;
rs->f[1][c].re = es->f[1][c].re - rs->f[1][c].re;
rs->f[1][c].im = es->f[1][c].im - rs->f[1][c].im;
rs->f[2][c].re = es->f[2][c].re - rs->f[2][c].re;
rs->f[2][c].im = es->f[2][c].im - rs->f[2][c].im;
rs->f[3][c].re = es->f[3][c].re - rs->f[3][c].re;
rs->f[3][c].im = es->f[3][c].im - rs->f[3][c].im;
#line 3664 "dwf.nw"
}
#line 3602 "dwf.nw"
}
#line 3586 "dwf.nw"
}
#line 3134 "dwf.nw"
    
#line 3263 "dwf.nw"
#undef qs
#undef qx5
#line 3135 "dwf.nw"
}
#line 3140 "dwf.nw"
static void
compute_Qxx1Qxy(vFermion *chi,
                const vFermion *psi,
                struct neighbor *nb)
{
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 3146 "dwf.nw"
    
#line 4022 "dwf.nw"
int xyzt, k, d;
const vFermion *f;
vHalfFermion *g;
vHalfFermion gg[2*DIM], hh[2*DIM];
vSU3 V[2*DIM];
int ps[2*DIM], p5[2*DIM];
#line 4030 "dwf.nw"
const SU3 *Uup, *Udown;
int c1, c2;
#line 3147 "dwf.nw"
    
#line 2767 "dwf.nw"
vReal fx;
vHalfFermion zV;
vector_complex zn;
scalar_complex zX[Fd/2][Nc];
#line 2955 "dwf.nw"
scalar_complex yOut[Fd/2][Nc];

#line 3149 "dwf.nw"
    
#line 3259 "dwf.nw"
#define qx5 rx5
#define qs rs
#line 3150 "dwf.nw"
    
#line 3326 "dwf.nw"
{
   int k, i, s, c, *src;
   const vFermion *f;
   vHalfFermion *g;

   k = 0; 
#line 3359 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3363 "dwf.nw"
        }
    }
}
#line 3332 "dwf.nw"
   k = 1; 
#line 3368 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3372 "dwf.nw"
        }
    }
}
#line 3333 "dwf.nw"
   k = 2; 
#line 3377 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3381 "dwf.nw"
        }
    }
}
#line 3334 "dwf.nw"
   k = 3; 
#line 3386 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3390 "dwf.nw"
        }
    }
}
#line 3335 "dwf.nw"
   k = 4; 
#line 3395 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3399 "dwf.nw"
        }
    }
}
#line 3336 "dwf.nw"
   k = 5; 
#line 3404 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
     for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3408 "dwf.nw"
        }
    }
}
#line 3337 "dwf.nw"
   k = 6; 
#line 3413 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3417 "dwf.nw"
        }
    }
}
#line 3338 "dwf.nw"
   k = 7; 
#line 3422 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3426 "dwf.nw"
        }
    }
}
#line 3339 "dwf.nw"
}
#line 3151 "dwf.nw"
    
#line 4079 "dwf.nw"
if (nb->qmp_smask) {
    QMP_start(nb->qmp_handle);
    DEBUG_QMP("start sends and receives (0x%x)\n", (int)nb->qmp_handle)
}
#line 3152 "dwf.nw"
    
#line 3680 "dwf.nw"
for (i = 0; i < nb->inside_size; i++) {
    xyzt = nb->inside[i];
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3684 "dwf.nw"
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3685 "dwf.nw"
    
#line 3455 "dwf.nw"
for (s = 0; s < Sv; s++) {
    
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3472 "dwf.nw"
for (c = 0; c < Nc; c++) {
    k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3474 "dwf.nw"
    k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3475 "dwf.nw"
    k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3476 "dwf.nw"
    k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3477 "dwf.nw"
    k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3478 "dwf.nw"
    k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3479 "dwf.nw"
    k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3480 "dwf.nw"
    k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3481 "dwf.nw"
}
#line 3457 "dwf.nw"
    
#line 3528 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3533 "dwf.nw"
}
#line 3458 "dwf.nw"
    
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3459 "dwf.nw"
}
#line 3686 "dwf.nw"
    
#line 2774 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2777 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2781 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2786 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2789 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2922 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2969 "dwf.nw"
  COMPUTE_YA(0,0,0) COMPUTE_YA(0,0,1) COMPUTE_YA(0,0,2)
#line 2925 "dwf.nw"
    
#line 2973 "dwf.nw"
  COMPUTE_YA(1,1,0) COMPUTE_YA(1,1,1) COMPUTE_YA(1,1,2)
#line 2926 "dwf.nw"
}
#line 2881 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2884 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
  rs = &rx5[s];
  QSETUP(s)
  
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2888 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2893 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2895 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2946 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2993 "dwf.nw"
  COMPUTE_YB(2,0,0) COMPUTE_YB(2,0,1) COMPUTE_YB(2,0,2)
#line 2949 "dwf.nw"
    
#line 2997 "dwf.nw"
  COMPUTE_YB(3,1,0) COMPUTE_YB(3,1,1) COMPUTE_YB(3,1,2)
#line 2950 "dwf.nw"
}
#line 3687 "dwf.nw"
}
#line 3153 "dwf.nw"
    
#line 4088 "dwf.nw"
if (nb->qmp_smask) {
    QMP_wait(nb->qmp_handle);
    DEBUG_QMP("waiting for sends and receives (0x%x)\n",
               (int)nb->qmp_handle)
}
#line 3154 "dwf.nw"
    
#line 3691 "dwf.nw"
for (i = 0; i < nb->boundary_size; i++) {
    int m = nb->boundary[i].mask;

    xyzt = nb->boundary[i].index;
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3697 "dwf.nw"
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3698 "dwf.nw"
    
#line 3463 "dwf.nw"
for (s = 0; s < Sv; s++) {
  
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3486 "dwf.nw"
for (c = 0; c < 3; c++) {
    if ((m & 0x01) == 0) {
        k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3489 "dwf.nw"
    }
    if ((m & 0x02) == 0) {
        k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3492 "dwf.nw"
    }
    if ((m & 0x04) == 0) {
        k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3495 "dwf.nw"
    }
    if ((m & 0x08) == 0) {
        k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3498 "dwf.nw"
    }
    if ((m & 0x10) == 0) {
        k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3501 "dwf.nw"
    }
    if ((m & 0x20) == 0) {
        k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3504 "dwf.nw"
    }
    if ((m & 0x40) == 0) {
        k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3507 "dwf.nw"
    }
    if ((m & 0x80) == 0) {
        k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3510 "dwf.nw"
    }
}
#line 3465 "dwf.nw"
  
#line 3537 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = (m & (1 << d))? &nb->rcv_buf[d][ps[d]]: &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3542 "dwf.nw"
}
#line 3466 "dwf.nw"
  
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3467 "dwf.nw"
}
#line 3699 "dwf.nw"
    
#line 2774 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2777 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2781 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2786 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2789 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2922 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2969 "dwf.nw"
  COMPUTE_YA(0,0,0) COMPUTE_YA(0,0,1) COMPUTE_YA(0,0,2)
#line 2925 "dwf.nw"
    
#line 2973 "dwf.nw"
  COMPUTE_YA(1,1,0) COMPUTE_YA(1,1,1) COMPUTE_YA(1,1,2)
#line 2926 "dwf.nw"
}
#line 2881 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2884 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
  rs = &rx5[s];
  QSETUP(s)
  
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2888 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2893 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2895 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2946 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2993 "dwf.nw"
  COMPUTE_YB(2,0,0) COMPUTE_YB(2,0,1) COMPUTE_YB(2,0,2)
#line 2949 "dwf.nw"
    
#line 2997 "dwf.nw"
  COMPUTE_YB(3,1,0) COMPUTE_YB(3,1,1) COMPUTE_YB(3,1,2)
#line 2950 "dwf.nw"
}
#line 3700 "dwf.nw"
}
#line 3155 "dwf.nw"
    
#line 3263 "dwf.nw"
#undef qs
#undef qx5
#line 3156 "dwf.nw"
}
#line 3161 "dwf.nw"
static void
compute_Sxx1Sxy(vFermion *chi,
                const vFermion *psi,
                struct neighbor *nb)
{
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 3167 "dwf.nw"
    
#line 4022 "dwf.nw"
int xyzt, k, d;
const vFermion *f;
vHalfFermion *g;
vHalfFermion gg[2*DIM], hh[2*DIM];
vSU3 V[2*DIM];
int ps[2*DIM], p5[2*DIM];
#line 4030 "dwf.nw"
const SU3 *Uup, *Udown;
int c1, c2;
#line 3168 "dwf.nw"
    
#line 2767 "dwf.nw"
vReal fx;
vHalfFermion zV;
vector_complex zn;
scalar_complex zX[Fd/2][Nc];
#line 2955 "dwf.nw"
scalar_complex yOut[Fd/2][Nc];

#line 3170 "dwf.nw"
    
#line 3259 "dwf.nw"
#define qx5 rx5
#define qs rs
#line 3171 "dwf.nw"
    
#line 3343 "dwf.nw"
{
   int k, i, s, c, *src;
   const vFermion *f;
   vHalfFermion *g;

   k = 0; 
#line 3368 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3372 "dwf.nw"
        }
    }
}
#line 3349 "dwf.nw"
   k = 1; 
#line 3359 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3363 "dwf.nw"
        }
    }
}
#line 3350 "dwf.nw"
   k = 2; 
#line 3386 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3390 "dwf.nw"
        }
    }
}
#line 3351 "dwf.nw"
   k = 3; 
#line 3377 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3381 "dwf.nw"
        }
    }
}
#line 3352 "dwf.nw"
   k = 4; 
#line 3404 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
     for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3408 "dwf.nw"
        }
    }
}
#line 3353 "dwf.nw"
   k = 5; 
#line 3395 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3399 "dwf.nw"
        }
    }
}
#line 3354 "dwf.nw"
   k = 6; 
#line 3422 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3426 "dwf.nw"
        }
    }
}
#line 3355 "dwf.nw"
   k = 7; 
#line 3413 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3417 "dwf.nw"
        }
    }
}
#line 3356 "dwf.nw"
}
#line 3172 "dwf.nw"
    
#line 4079 "dwf.nw"
if (nb->qmp_smask) {
    QMP_start(nb->qmp_handle);
    DEBUG_QMP("start sends and receives (0x%x)\n", (int)nb->qmp_handle)
}
#line 3173 "dwf.nw"
    
#line 3705 "dwf.nw"
for (i = 0; i < nb->inside_size; i++) {
    xyzt = nb->inside[i];
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3709 "dwf.nw"
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3710 "dwf.nw"
    
#line 3729 "dwf.nw"
for (s = 0; s < Sv; s++) {
    
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3607 "dwf.nw"
for (c = 0; c < Nc; c++) {
    k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3609 "dwf.nw"
    k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3610 "dwf.nw"
    k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3611 "dwf.nw"
    k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3612 "dwf.nw"
    k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3613 "dwf.nw"
    k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3614 "dwf.nw"
    k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3615 "dwf.nw"
    k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3616 "dwf.nw"
}
#line 3731 "dwf.nw"
    
#line 3528 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3533 "dwf.nw"
}
#line 3732 "dwf.nw"
    
#line 3745 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 6; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3748 "dwf.nw"
    k = 7; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3749 "dwf.nw"
    k = 2; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3750 "dwf.nw"
    k = 3; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3751 "dwf.nw"
    k = 1; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3752 "dwf.nw"
    k = 0; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3753 "dwf.nw"
    k = 5; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3754 "dwf.nw"
    k = 4; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3755 "dwf.nw"
}
#line 3733 "dwf.nw"
}
#line 3711 "dwf.nw"
    
#line 2815 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2818 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2822 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2827 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2829 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2930 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2977 "dwf.nw"
  COMPUTE_YA(2,0,0) COMPUTE_YA(2,0,1) COMPUTE_YA(2,0,2)
#line 2933 "dwf.nw"
    
#line 2981 "dwf.nw"
  COMPUTE_YA(3,1,0) COMPUTE_YA(3,1,1) COMPUTE_YA(3,1,2)
#line 2934 "dwf.nw"
}
#line 2853 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2856 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2860 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2865 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2868 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2938 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2985 "dwf.nw"
  COMPUTE_YB(0,0,0) COMPUTE_YB(0,0,1) COMPUTE_YB(0,0,2)
#line 2941 "dwf.nw"
    
#line 2989 "dwf.nw"
  COMPUTE_YB(1,1,0) COMPUTE_YB(1,1,1) COMPUTE_YB(1,1,2)
#line 2942 "dwf.nw"
}
#line 3712 "dwf.nw"
}
#line 3174 "dwf.nw"
    
#line 4088 "dwf.nw"
if (nb->qmp_smask) {
    QMP_wait(nb->qmp_handle);
    DEBUG_QMP("waiting for sends and receives (0x%x)\n",
               (int)nb->qmp_handle)
}
#line 3175 "dwf.nw"
    
#line 3716 "dwf.nw"
for (i = 0; i < nb->boundary_size; i++) {
    int m = nb->boundary[i].mask;

    xyzt = nb->boundary[i].index;
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3722 "dwf.nw"
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3723 "dwf.nw"
    
#line 3737 "dwf.nw"
for (s = 0; s < Sv; s++) {
  

#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3623 "dwf.nw"
for (c = 0; c < Nc; c++) {
    if ((m & 0x01) == 0) {
        k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3626 "dwf.nw"
    }
    if ((m & 0x02) == 0) {
        k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3629 "dwf.nw"
    }
    if ((m & 0x04) == 0) {
        k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3632 "dwf.nw"
    }
    if ((m & 0x08) == 0) {
        k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3635 "dwf.nw"
    }
    if ((m & 0x10) == 0) {
        k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3638 "dwf.nw"
    }
    if ((m & 0x20) == 0) {
        k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3641 "dwf.nw"
    }
    if ((m & 0x40) == 0) {
        k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3644 "dwf.nw"
    }
    if ((m & 0x80) == 0) {
        k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3647 "dwf.nw"
    }
}
#line 3739 "dwf.nw"
  
#line 3537 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = (m & (1 << d))? &nb->rcv_buf[d][ps[d]]: &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3542 "dwf.nw"
}
#line 3740 "dwf.nw"
  
#line 3745 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 6; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3748 "dwf.nw"
    k = 7; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3749 "dwf.nw"
    k = 2; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3750 "dwf.nw"
    k = 3; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3751 "dwf.nw"
    k = 1; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3752 "dwf.nw"
    k = 0; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3753 "dwf.nw"
    k = 5; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3754 "dwf.nw"
    k = 4; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3755 "dwf.nw"
}
#line 3741 "dwf.nw"
}
#line 3724 "dwf.nw"
    
#line 2815 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2818 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2822 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2827 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2829 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2930 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2977 "dwf.nw"
  COMPUTE_YA(2,0,0) COMPUTE_YA(2,0,1) COMPUTE_YA(2,0,2)
#line 2933 "dwf.nw"
    
#line 2981 "dwf.nw"
  COMPUTE_YA(3,1,0) COMPUTE_YA(3,1,1) COMPUTE_YA(3,1,2)
#line 2934 "dwf.nw"
}
#line 2853 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2856 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2860 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2865 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2868 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2938 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2985 "dwf.nw"
  COMPUTE_YB(0,0,0) COMPUTE_YB(0,0,1) COMPUTE_YB(0,0,2)
#line 2941 "dwf.nw"
    
#line 2989 "dwf.nw"
  COMPUTE_YB(1,1,0) COMPUTE_YB(1,1,1) COMPUTE_YB(1,1,2)
#line 2942 "dwf.nw"
}
#line 3725 "dwf.nw"
}
#line 3176 "dwf.nw"
    
#line 3263 "dwf.nw"
#undef qs
#undef qx5
#line 3177 "dwf.nw"
}
#line 3183 "dwf.nw"
static void
compute_1Qxx1Qxy(vFermion *chi,
                 double *norm,
                 const vFermion *eta,
                 const vFermion *psi,
                 struct neighbor *nb)
{
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 3191 "dwf.nw"
    
#line 4022 "dwf.nw"
int xyzt, k, d;
const vFermion *f;
vHalfFermion *g;
vHalfFermion gg[2*DIM], hh[2*DIM];
vSU3 V[2*DIM];
int ps[2*DIM], p5[2*DIM];
#line 4030 "dwf.nw"
const SU3 *Uup, *Udown;
int c1, c2;
#line 3192 "dwf.nw"
    
#line 2767 "dwf.nw"
vReal fx;
vHalfFermion zV;
vector_complex zn;
scalar_complex zX[Fd/2][Nc];
#line 2955 "dwf.nw"
scalar_complex yOut[Fd/2][Nc];
#line 3193 "dwf.nw"
    vReal vv;
    vReal nv;
    *norm = 0;

    
#line 3259 "dwf.nw"
#define qx5 rx5
#define qs rs
#line 3198 "dwf.nw"
    
#line 3326 "dwf.nw"
{
   int k, i, s, c, *src;
   const vFermion *f;
   vHalfFermion *g;

   k = 0; 
#line 3359 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3363 "dwf.nw"
        }
    }
}
#line 3332 "dwf.nw"
   k = 1; 
#line 3368 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3372 "dwf.nw"
        }
    }
}
#line 3333 "dwf.nw"
   k = 2; 
#line 3377 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3381 "dwf.nw"
        }
    }
}
#line 3334 "dwf.nw"
   k = 3; 
#line 3386 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3390 "dwf.nw"
        }
    }
}
#line 3335 "dwf.nw"
   k = 4; 
#line 3395 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3399 "dwf.nw"
        }
    }
}
#line 3336 "dwf.nw"
   k = 5; 
#line 3404 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
     for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3408 "dwf.nw"
        }
    }
}
#line 3337 "dwf.nw"
   k = 6; 
#line 3413 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3417 "dwf.nw"
        }
    }
}
#line 3338 "dwf.nw"
   k = 7; 
#line 3422 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3426 "dwf.nw"
        }
    }
}
#line 3339 "dwf.nw"
}
#line 3199 "dwf.nw"
    
#line 4079 "dwf.nw"
if (nb->qmp_smask) {
    QMP_start(nb->qmp_handle);
    DEBUG_QMP("start sends and receives (0x%x)\n", (int)nb->qmp_handle)
}
#line 3200 "dwf.nw"
    
#line 3760 "dwf.nw"
for (i = 0; i < nb->inside_size; i++) {
    const vFermion *ex5, *es;

    xyzt = nb->inside[i];
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3766 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3768 "dwf.nw"
    
#line 3455 "dwf.nw"
for (s = 0; s < Sv; s++) {
    
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3472 "dwf.nw"
for (c = 0; c < Nc; c++) {
    k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3474 "dwf.nw"
    k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3475 "dwf.nw"
    k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3476 "dwf.nw"
    k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3477 "dwf.nw"
    k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3478 "dwf.nw"
    k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3479 "dwf.nw"
    k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3480 "dwf.nw"
    k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3481 "dwf.nw"
}
#line 3457 "dwf.nw"
    
#line 3528 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3533 "dwf.nw"
}
#line 3458 "dwf.nw"
    
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3459 "dwf.nw"
}
#line 3769 "dwf.nw"
    
#line 2774 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2777 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2781 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2786 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2789 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2922 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2969 "dwf.nw"
  COMPUTE_YA(0,0,0) COMPUTE_YA(0,0,1) COMPUTE_YA(0,0,2)
#line 2925 "dwf.nw"
    
#line 2973 "dwf.nw"
  COMPUTE_YA(1,1,0) COMPUTE_YA(1,1,1) COMPUTE_YA(1,1,2)
#line 2926 "dwf.nw"
}
#line 2881 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2884 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
  rs = &rx5[s];
  QSETUP(s)
  
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2888 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2893 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2895 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2946 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2993 "dwf.nw"
  COMPUTE_YB(2,0,0) COMPUTE_YB(2,0,1) COMPUTE_YB(2,0,2)
#line 2949 "dwf.nw"
    
#line 2997 "dwf.nw"
  COMPUTE_YB(3,1,0) COMPUTE_YB(3,1,1) COMPUTE_YB(3,1,2)
#line 2950 "dwf.nw"
}
#line 3790 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    es = &ex5[s];
    nv = vmk_1(0.0);
    for (c = 0; c < Nc; c++) {
        
#line 3802 "dwf.nw"
vv = es->f[0][c].re - rs->f[0][c].re; rs->f[0][c].re = vv; nv += vv * vv;
vv = es->f[0][c].im - rs->f[0][c].im; rs->f[0][c].im = vv; nv += vv * vv;
vv = es->f[1][c].re - rs->f[1][c].re; rs->f[1][c].re = vv; nv += vv * vv;
vv = es->f[1][c].im - rs->f[1][c].im; rs->f[1][c].im = vv; nv += vv * vv;
vv = es->f[2][c].re - rs->f[2][c].re; rs->f[2][c].re = vv; nv += vv * vv;
vv = es->f[2][c].im - rs->f[2][c].im; rs->f[2][c].im = vv; nv += vv * vv;
vv = es->f[3][c].re - rs->f[3][c].re; rs->f[3][c].re = vv; nv += vv * vv;
vv = es->f[3][c].im - rs->f[3][c].im; rs->f[3][c].im = vv; nv += vv * vv;
#line 3796 "dwf.nw"
    }
    *norm += vsum(nv);
}
#line 3770 "dwf.nw"
}
#line 3201 "dwf.nw"
    
#line 4088 "dwf.nw"
if (nb->qmp_smask) {
    QMP_wait(nb->qmp_handle);
    DEBUG_QMP("waiting for sends and receives (0x%x)\n",
               (int)nb->qmp_handle)
}
#line 3202 "dwf.nw"
    
#line 3774 "dwf.nw"
for (i = 0; i < nb->boundary_size; i++) {
    const vFermion *ex5, *es;
    int m = nb->boundary[i].mask;

    xyzt = nb->boundary[i].index;
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3781 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3783 "dwf.nw"
    
#line 3463 "dwf.nw"
for (s = 0; s < Sv; s++) {
  
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3486 "dwf.nw"
for (c = 0; c < 3; c++) {
    if ((m & 0x01) == 0) {
        k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3489 "dwf.nw"
    }
    if ((m & 0x02) == 0) {
        k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3492 "dwf.nw"
    }
    if ((m & 0x04) == 0) {
        k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3495 "dwf.nw"
    }
    if ((m & 0x08) == 0) {
        k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3498 "dwf.nw"
    }
    if ((m & 0x10) == 0) {
        k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3501 "dwf.nw"
    }
    if ((m & 0x20) == 0) {
        k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3504 "dwf.nw"
    }
    if ((m & 0x40) == 0) {
        k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3507 "dwf.nw"
    }
    if ((m & 0x80) == 0) {
        k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3510 "dwf.nw"
    }
}
#line 3465 "dwf.nw"
  
#line 3537 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = (m & (1 << d))? &nb->rcv_buf[d][ps[d]]: &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3542 "dwf.nw"
}
#line 3466 "dwf.nw"
  
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3467 "dwf.nw"
}
#line 3784 "dwf.nw"
    
#line 2774 "dwf.nw"
vhfzero(&zV);
fx = vfx_A;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2777 "dwf.nw"
for (s = 0; s < Sv_1; s++, fx = fx * vab) {
    rs = &rx5[s];
    QSETUP(s)
    
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2781 "dwf.nw"
}
rs = &rx5[Sv_1];
QSETUP(Sv_1)
fx = vput_n(fx, c0);
#line 2803 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[0][c].re; Q2R(0,re)
    zV.f[0][c].im += fx * qs->f[0][c].im; Q2R(0,im)
    zV.f[1][c].re += fx * qs->f[1][c].re; Q2R(1,re)
    zV.f[1][c].im += fx * qs->f[1][c].im; Q2R(1,im)
}
#line 2786 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);

#line 2789 "dwf.nw"
  zn.re = qs->f[0][c].re;               zn.im = qs->f[0][c].im;
  zn.re = vput_n(zn.re, zX[0][c].re);   zn.im = vput_n(zn.im, zX[0][c].im);
  rs->f[0][c].re = zn.re;               rs->f[0][c].im = zn.im;

  zn.re = qs->f[1][c].re;               zn.im = qs->f[1][c].im;
  zn.re = vput_n(zn.re, zX[1][c].re);   zn.im = vput_n(zn.im, zX[1][c].im);
  rs->f[1][c].re = zn.re;               rs->f[1][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2922 "dwf.nw"
for (s = Sv; s--;) {
    rs = &rx5[s];
    
#line 2969 "dwf.nw"
  COMPUTE_YA(0,0,0) COMPUTE_YA(0,0,1) COMPUTE_YA(0,0,2)
#line 2925 "dwf.nw"
    
#line 2973 "dwf.nw"
  COMPUTE_YA(1,1,0) COMPUTE_YA(1,1,1) COMPUTE_YA(1,1,2)
#line 2926 "dwf.nw"
}
#line 2881 "dwf.nw"
vhfzero(&zV);
fx = vfx_B;
#line 2753 "dwf.nw"
#if defined(qs)
#define QSETUP(s)
#define Q2R(d,pt)
#else
#define QSETUP(s) qs = &qx5[s];
#define Q2R(d,pt) rs->f[d][c].pt = qs->f[d][c].pt;
#endif
#line 2884 "dwf.nw"
for (s = Sv; --s; fx = fx * vab) {
  rs = &rx5[s];
  QSETUP(s)
  
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2888 "dwf.nw"
}
rs = &rx5[0];
QSETUP(0)
fx = vput_0(fx, c0);
#line 2841 "dwf.nw"
for (c = 0; c < Nc; c++) {
    zV.f[0][c].re += fx * qs->f[2][c].re; Q2R(2,re)
    zV.f[0][c].im += fx * qs->f[2][c].im; Q2R(2,im)
    zV.f[1][c].re += fx * qs->f[3][c].re; Q2R(3,re)
    zV.f[1][c].im += fx * qs->f[3][c].im; Q2R(3,im)
}
#line 2893 "dwf.nw"
for (c = 0; c < Nc; c++) {
  
#line 2908 "dwf.nw"
zX[0][c].re = vsum(zV.f[0][c].re);
zX[0][c].im = vsum(zV.f[0][c].im);
zX[1][c].re = vsum(zV.f[1][c].re);
zX[1][c].im = vsum(zV.f[1][c].im);
#line 2895 "dwf.nw"
  
  zn.re = qs->f[2][c].re;               zn.im = qs->f[2][c].im;
  zn.re = vput_0(zn.re, zX[0][c].re);   zn.im = vput_0(zn.im, zX[0][c].im);
  rs->f[2][c].re = zn.re;               rs->f[2][c].im = zn.im;

  zn.re = qs->f[3][c].re;               zn.im = qs->f[3][c].im;
  zn.re = vput_0(zn.re, zX[1][c].re);   zn.im = vput_0(zn.im, zX[1][c].im);
  rs->f[3][c].re = zn.re;               rs->f[3][c].im = zn.im;
}
#line 2762 "dwf.nw"
#undef QSETUP
#undef Q2R
#line 2958 "dwf.nw"
yOut[0][0].re = yOut[0][0].im = 0;
yOut[0][1].re = yOut[0][1].im = 0;
yOut[0][2].re = yOut[0][2].im = 0;
yOut[1][0].re = yOut[1][0].im = 0;
yOut[1][1].re = yOut[1][1].im = 0;
yOut[1][2].re = yOut[1][2].im = 0;
#line 2946 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    
#line 2993 "dwf.nw"
  COMPUTE_YB(2,0,0) COMPUTE_YB(2,0,1) COMPUTE_YB(2,0,2)
#line 2949 "dwf.nw"
    
#line 2997 "dwf.nw"
  COMPUTE_YB(3,1,0) COMPUTE_YB(3,1,1) COMPUTE_YB(3,1,2)
#line 2950 "dwf.nw"
}
#line 3790 "dwf.nw"
for (s = 0; s < Sv; s++) {
    rs = &rx5[s];
    es = &ex5[s];
    nv = vmk_1(0.0);
    for (c = 0; c < Nc; c++) {
        
#line 3802 "dwf.nw"
vv = es->f[0][c].re - rs->f[0][c].re; rs->f[0][c].re = vv; nv += vv * vv;
vv = es->f[0][c].im - rs->f[0][c].im; rs->f[0][c].im = vv; nv += vv * vv;
vv = es->f[1][c].re - rs->f[1][c].re; rs->f[1][c].re = vv; nv += vv * vv;
vv = es->f[1][c].im - rs->f[1][c].im; rs->f[1][c].im = vv; nv += vv * vv;
vv = es->f[2][c].re - rs->f[2][c].re; rs->f[2][c].re = vv; nv += vv * vv;
vv = es->f[2][c].im - rs->f[2][c].im; rs->f[2][c].im = vv; nv += vv * vv;
vv = es->f[3][c].re - rs->f[3][c].re; rs->f[3][c].re = vv; nv += vv * vv;
vv = es->f[3][c].im - rs->f[3][c].im; rs->f[3][c].im = vv; nv += vv * vv;
#line 3796 "dwf.nw"
    }
    *norm += vsum(nv);
}
#line 3785 "dwf.nw"
}
#line 3203 "dwf.nw"
    
#line 4112 "dwf.nw"
QMP_sum_double(norm);
DEBUG_QMP("sum_double(%p)\n", norm)
#line 3204 "dwf.nw"
    
#line 3263 "dwf.nw"
#undef qs
#undef qx5
#line 3205 "dwf.nw"
}
#line 3210 "dwf.nw"
static void
compute_Dx(vFermion *chi,
           const vFermion *eta,
           const vFermion *psi,
           struct neighbor *nb)
{
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 3217 "dwf.nw"
    
#line 4022 "dwf.nw"
int xyzt, k, d;
const vFermion *f;
vHalfFermion *g;
vHalfFermion gg[2*DIM], hh[2*DIM];
vSU3 V[2*DIM];
int ps[2*DIM], p5[2*DIM];
#line 4030 "dwf.nw"
const SU3 *Uup, *Udown;
int c1, c2;
#line 3218 "dwf.nw"
    
#line 3952 "dwf.nw"
const vFermion *es1;
vReal vbc;

#line 3220 "dwf.nw"
    
#line 3259 "dwf.nw"
#define qx5 rx5
#define qs rs
#line 3221 "dwf.nw"
    
#line 3326 "dwf.nw"
{
   int k, i, s, c, *src;
   const vFermion *f;
   vHalfFermion *g;

   k = 0; 
#line 3359 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3363 "dwf.nw"
        }
    }
}
#line 3332 "dwf.nw"
   k = 1; 
#line 3368 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3372 "dwf.nw"
        }
    }
}
#line 3333 "dwf.nw"
   k = 2; 
#line 3377 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3381 "dwf.nw"
        }
    }
}
#line 3334 "dwf.nw"
   k = 3; 
#line 3386 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3390 "dwf.nw"
        }
    }
}
#line 3335 "dwf.nw"
   k = 4; 
#line 3395 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3399 "dwf.nw"
        }
    }
}
#line 3336 "dwf.nw"
   k = 5; 
#line 3404 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
     for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3408 "dwf.nw"
        }
    }
}
#line 3337 "dwf.nw"
   k = 6; 
#line 3413 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3417 "dwf.nw"
        }
    }
}
#line 3338 "dwf.nw"
   k = 7; 
#line 3422 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3426 "dwf.nw"
        }
    }
}
#line 3339 "dwf.nw"
}
#line 3222 "dwf.nw"
    
#line 4079 "dwf.nw"
if (nb->qmp_smask) {
    QMP_start(nb->qmp_handle);
    DEBUG_QMP("start sends and receives (0x%x)\n", (int)nb->qmp_handle)
}
#line 3223 "dwf.nw"
    
#line 3814 "dwf.nw"
for (i = 0; i < nb->inside_size; i++) {
    const vFermion *ex5, *es;

    xyzt = nb->inside[i];
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3820 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3822 "dwf.nw"
    
#line 3455 "dwf.nw"
for (s = 0; s < Sv; s++) {
    
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3472 "dwf.nw"
for (c = 0; c < Nc; c++) {
    k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3474 "dwf.nw"
    k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3475 "dwf.nw"
    k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3476 "dwf.nw"
    k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3477 "dwf.nw"
    k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3478 "dwf.nw"
    k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3479 "dwf.nw"
    k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3480 "dwf.nw"
    k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3481 "dwf.nw"
}
#line 3457 "dwf.nw"
    
#line 3528 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3533 "dwf.nw"
}
#line 3458 "dwf.nw"
    
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3459 "dwf.nw"
}
#line 3823 "dwf.nw"
    
#line 3883 "dwf.nw"
for (s = Sv, vbc = vbnc, es1 = &ex5[0]; s--; vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_up1(es->f[d][c].r, es1->f[d][c].r)
    QXX(0,0,re); QXX(0,0,im);
    QXX(0,1,re); QXX(0,1,im);
    QXX(0,2,re); QXX(0,2,im);
    QXX(1,0,re); QXX(1,0,im);
    QXX(1,1,re); QXX(1,1,im);
    QXX(1,2,re); QXX(1,2,im);
#undef QXX
    es1 = es;
}
#line 3934 "dwf.nw"
for (s = 0, vbc = vcbn, es1 = &ex5[Sv_1]; s < Sv; s++, vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_upN(es1->f[d][c].r, es->f[d][c].r)
    QXX(2,0,re); QXX(2,0,im);
    QXX(2,1,re); QXX(2,1,im);
    QXX(2,2,re); QXX(2,2,im);
    QXX(3,0,re); QXX(3,0,im);
    QXX(3,1,re); QXX(3,1,im);
    QXX(3,2,re); QXX(3,2,im);
#undef QXX
    es1 = es;
}
#line 3824 "dwf.nw"
}
#line 3224 "dwf.nw"
    
#line 4088 "dwf.nw"
if (nb->qmp_smask) {
    QMP_wait(nb->qmp_handle);
    DEBUG_QMP("waiting for sends and receives (0x%x)\n",
               (int)nb->qmp_handle)
}
#line 3225 "dwf.nw"
    
#line 3828 "dwf.nw"
for (i = 0; i < nb->boundary_size; i++) {
    const vFermion *ex5, *es;
    int m = nb->boundary[i].mask;

    xyzt = nb->boundary[i].index;
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3835 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3837 "dwf.nw"
    
#line 3463 "dwf.nw"
for (s = 0; s < Sv; s++) {
  
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3486 "dwf.nw"
for (c = 0; c < 3; c++) {
    if ((m & 0x01) == 0) {
        k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3489 "dwf.nw"
    }
    if ((m & 0x02) == 0) {
        k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3492 "dwf.nw"
    }
    if ((m & 0x04) == 0) {
        k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3495 "dwf.nw"
    }
    if ((m & 0x08) == 0) {
        k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3498 "dwf.nw"
    }
    if ((m & 0x10) == 0) {
        k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3501 "dwf.nw"
    }
    if ((m & 0x20) == 0) {
        k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3504 "dwf.nw"
    }
    if ((m & 0x40) == 0) {
        k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3507 "dwf.nw"
    }
    if ((m & 0x80) == 0) {
        k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3510 "dwf.nw"
    }
}
#line 3465 "dwf.nw"
  
#line 3537 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = (m & (1 << d))? &nb->rcv_buf[d][ps[d]]: &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3542 "dwf.nw"
}
#line 3466 "dwf.nw"
  
#line 3514 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 7; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3517 "dwf.nw"
    k = 6; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3518 "dwf.nw"
    k = 3; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3519 "dwf.nw"
    k = 2; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3520 "dwf.nw"
    k = 1; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3521 "dwf.nw"
    k = 0; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3522 "dwf.nw"
    k = 5; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3523 "dwf.nw"
    k = 4; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3524 "dwf.nw"
}
#line 3467 "dwf.nw"
}
#line 3838 "dwf.nw"
    
#line 3883 "dwf.nw"
for (s = Sv, vbc = vbnc, es1 = &ex5[0]; s--; vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_up1(es->f[d][c].r, es1->f[d][c].r)
    QXX(0,0,re); QXX(0,0,im);
    QXX(0,1,re); QXX(0,1,im);
    QXX(0,2,re); QXX(0,2,im);
    QXX(1,0,re); QXX(1,0,im);
    QXX(1,1,re); QXX(1,1,im);
    QXX(1,2,re); QXX(1,2,im);
#undef QXX
    es1 = es;
}
#line 3934 "dwf.nw"
for (s = 0, vbc = vcbn, es1 = &ex5[Sv_1]; s < Sv; s++, vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_upN(es1->f[d][c].r, es->f[d][c].r)
    QXX(2,0,re); QXX(2,0,im);
    QXX(2,1,re); QXX(2,1,im);
    QXX(2,2,re); QXX(2,2,im);
    QXX(3,0,re); QXX(3,0,im);
    QXX(3,1,re); QXX(3,1,im);
    QXX(3,2,re); QXX(3,2,im);
#undef QXX
    es1 = es;
}
#line 3839 "dwf.nw"
}
#line 3226 "dwf.nw"
    
#line 3263 "dwf.nw"
#undef qs
#undef qx5
#line 3227 "dwf.nw"
}
#line 3233 "dwf.nw"
static void
compute_Dcx(vFermion *chi,
            const vFermion *eta,
            const vFermion *psi,
            struct neighbor *nb)
{
    
#line 4017 "dwf.nw"
int i, xyzt5, s, c;
vFermion * __restrict__ rx5, * __restrict__ rs;
#line 3240 "dwf.nw"
    
#line 4022 "dwf.nw"
int xyzt, k, d;
const vFermion *f;
vHalfFermion *g;
vHalfFermion gg[2*DIM], hh[2*DIM];
vSU3 V[2*DIM];
int ps[2*DIM], p5[2*DIM];
#line 4030 "dwf.nw"
const SU3 *Uup, *Udown;
int c1, c2;
#line 3241 "dwf.nw"
    
#line 3952 "dwf.nw"
const vFermion *es1;
vReal vbc;

#line 3243 "dwf.nw"
    
#line 3259 "dwf.nw"
#define qx5 rx5
#define qs rs
#line 3244 "dwf.nw"
    
#line 3343 "dwf.nw"
{
   int k, i, s, c, *src;
   const vFermion *f;
   vHalfFermion *g;

   k = 0; 
#line 3368 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3372 "dwf.nw"
        }
    }
}
#line 3349 "dwf.nw"
   k = 1; 
#line 3359 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3363 "dwf.nw"
        }
    }
}
#line 3350 "dwf.nw"
   k = 2; 
#line 3386 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3390 "dwf.nw"
        }
    }
}
#line 3351 "dwf.nw"
   k = 3; 
#line 3377 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3381 "dwf.nw"
        }
    }
}
#line 3352 "dwf.nw"
   k = 4; 
#line 3404 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
     for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3408 "dwf.nw"
        }
    }
}
#line 3353 "dwf.nw"
   k = 5; 
#line 3395 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3399 "dwf.nw"
        }
    }
}
#line 3354 "dwf.nw"
   k = 6; 
#line 3422 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3426 "dwf.nw"
        }
    }
}
#line 3355 "dwf.nw"
   k = 7; 
#line 3413 "dwf.nw"
for (i = nb->snd_size[k], g = nb->snd_buf[k], src = nb->snd[k]; i--; src++) {
    for (s = Sv, f = &psi[*src]; s--; g++, f++) {
        for (c = 0; c < Nc; c++) {
            
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3417 "dwf.nw"
        }
    }
}
#line 3356 "dwf.nw"
}
#line 3245 "dwf.nw"
    
#line 4079 "dwf.nw"
if (nb->qmp_smask) {
    QMP_start(nb->qmp_handle);
    DEBUG_QMP("start sends and receives (0x%x)\n", (int)nb->qmp_handle)
}
#line 3246 "dwf.nw"
    
#line 3848 "dwf.nw"
for (i = 0; i < nb->inside_size; i++) {
    const vFermion *ex5, *es;

    xyzt = nb->inside[i];
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3854 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3856 "dwf.nw"
    
#line 3729 "dwf.nw"
for (s = 0; s < Sv; s++) {
    
#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3607 "dwf.nw"
for (c = 0; c < Nc; c++) {
    k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3609 "dwf.nw"
    k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3610 "dwf.nw"
    k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3611 "dwf.nw"
    k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3612 "dwf.nw"
    k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3613 "dwf.nw"
    k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3614 "dwf.nw"
    k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3615 "dwf.nw"
    k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3616 "dwf.nw"
}
#line 3731 "dwf.nw"
    
#line 3528 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3533 "dwf.nw"
}
#line 3732 "dwf.nw"
    
#line 3745 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 6; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3748 "dwf.nw"
    k = 7; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3749 "dwf.nw"
    k = 2; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3750 "dwf.nw"
    k = 3; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3751 "dwf.nw"
    k = 1; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3752 "dwf.nw"
    k = 0; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3753 "dwf.nw"
    k = 5; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3754 "dwf.nw"
    k = 4; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3755 "dwf.nw"
}
#line 3733 "dwf.nw"
}
#line 3857 "dwf.nw"
    
#line 3917 "dwf.nw"
for (s = 0, vbc = vcbn, es1 = &ex5[Sv_1]; s < Sv; s++, vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_upN(es1->f[d][c].r, es->f[d][c].r)
    QXX(0,0,re); QXX(0,0,im);
    QXX(0,1,re); QXX(0,1,im);
    QXX(0,2,re); QXX(0,2,im);
    QXX(1,0,re); QXX(1,0,im);
    QXX(1,1,re); QXX(1,1,im);
    QXX(1,2,re); QXX(1,2,im);
#undef QXX
    es1 = es;
}
#line 3900 "dwf.nw"
for (s = Sv, vbc = vbnc, es1 = &ex5[0]; s--; vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_up1(es->f[d][c].r, es1->f[d][c].r)
    QXX(2,0,re); QXX(2,0,im);
    QXX(2,1,re); QXX(2,1,im);
    QXX(2,2,re); QXX(2,2,im);
    QXX(3,0,re); QXX(3,0,im);
    QXX(3,1,re); QXX(3,1,im);
    QXX(3,2,re); QXX(3,2,im);
#undef QXX
    es1 = es;
}
#line 3858 "dwf.nw"
}
#line 3247 "dwf.nw"
    
#line 4088 "dwf.nw"
if (nb->qmp_smask) {
    QMP_wait(nb->qmp_handle);
    DEBUG_QMP("waiting for sends and receives (0x%x)\n",
               (int)nb->qmp_handle)
}
#line 3248 "dwf.nw"
    
#line 3862 "dwf.nw"
for (i = 0; i < nb->boundary_size; i++) {
    const vFermion *ex5, *es;
    int m = nb->boundary[i].mask;

    xyzt = nb->boundary[i].index;
    xyzt5 = xyzt * Sv;
    
#line 4037 "dwf.nw"
for (d = 0; d < 2*DIM; d++)
    p5[d] = nb->site[xyzt].F[d];
#line 4042 "dwf.nw"
rx5 = &chi[xyzt5];
#line 3869 "dwf.nw"
    ex5 = &eta[xyzt5];
    
#line 3960 "dwf.nw"
Uup = &U[nb->site[xyzt].Uup];
for (d = 0; d < DIM; d++, Uup++) {
    Udown = &U[nb->site[xyzt].Udown[d]];
    for (c1 = 0; c1 < Nc; c1++) {
        for (c2 = 0; c2 < Nc; c2++) {
	    /* conjugate down-link */
	    V[d*2+0].v[c1][c2].re = vmk_1( Udown->v[c2][c1].re);
	    V[d*2+0].v[c1][c2].im = vmk_1(-Udown->v[c2][c1].im);
	    /* normal up-link */
	    V[d*2+1].v[c1][c2].re = vmk_1(Uup->v[c1][c2].re);
	    V[d*2+1].v[c1][c2].im = vmk_1(Uup->v[c1][c2].im);
        }
    }
}
#line 3871 "dwf.nw"
    
#line 3737 "dwf.nw"
for (s = 0; s < Sv; s++) {
  

#line 3977 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
  ps[d] = p5[d] + s;
}
#line 3623 "dwf.nw"
for (c = 0; c < Nc; c++) {
    if ((m & 0x01) == 0) {
        k=0; f=&psi[ps[0]]; g=&gg[0]; 
#line 138 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].re;
#line 3626 "dwf.nw"
    }
    if ((m & 0x02) == 0) {
        k=1; f=&psi[ps[1]]; g=&gg[1]; 
#line 151 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].re;
#line 3629 "dwf.nw"
    }
    if ((m & 0x04) == 0) {
        k=2; f=&psi[ps[2]]; g=&gg[2]; 
#line 179 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[2][c].im;
#line 3632 "dwf.nw"
    }
    if ((m & 0x08) == 0) {
        k=3; f=&psi[ps[3]]; g=&gg[3]; 
#line 192 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[3][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[3][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[2][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[2][c].im;
#line 3635 "dwf.nw"
    }
    if ((m & 0x10) == 0) {
        k=4; f=&psi[ps[4]]; g=&gg[4]; 
#line 220 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].re;
#line 3638 "dwf.nw"
    }
    if ((m & 0x20) == 0) {
        k=5; f=&psi[ps[5]]; g=&gg[5]; 
#line 233 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].im;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].re;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].im;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].re;
#line 3641 "dwf.nw"
    }
    if ((m & 0x40) == 0) {
        k=6; f=&psi[ps[6]]; g=&gg[6]; 
#line 261 "dwf.nw"
g->f[0][c].re = f->f[0][c].re + f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im + f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re + f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im + f->f[3][c].im;
#line 3644 "dwf.nw"
    }
    if ((m & 0x80) == 0) {
        k=7; f=&psi[ps[7]]; g=&gg[7]; 
#line 274 "dwf.nw"
g->f[0][c].re = f->f[0][c].re - f->f[2][c].re;
g->f[0][c].im = f->f[0][c].im - f->f[2][c].im;
g->f[1][c].re = f->f[1][c].re - f->f[3][c].re;
g->f[1][c].im = f->f[1][c].im - f->f[3][c].im;
#line 3647 "dwf.nw"
    }
}
#line 3739 "dwf.nw"
  
#line 3537 "dwf.nw"
for (d = 0; d < 2*DIM; d++) {
    vHalfFermion * __restrict__ h = &hh[d];
    vSU3 *u = &V[d];
    g = (m & (1 << d))? &nb->rcv_buf[d][ps[d]]: &gg[d];
    
#line 3545 "dwf.nw"
for (c = 0; c < Nc; c++) {
    h->f[0][c].re=u->v[c][0].re*g->f[0][0].re-u->v[c][0].im*g->f[0][0].im
                 +u->v[c][1].re*g->f[0][1].re-u->v[c][1].im*g->f[0][1].im
                 +u->v[c][2].re*g->f[0][2].re-u->v[c][2].im*g->f[0][2].im;
    h->f[0][c].im=u->v[c][0].im*g->f[0][0].re+u->v[c][0].re*g->f[0][0].im
                 +u->v[c][1].im*g->f[0][1].re+u->v[c][1].re*g->f[0][1].im
                 +u->v[c][2].im*g->f[0][2].re+u->v[c][2].re*g->f[0][2].im;
    h->f[1][c].re=u->v[c][0].re*g->f[1][0].re-u->v[c][0].im*g->f[1][0].im
                 +u->v[c][1].re*g->f[1][1].re-u->v[c][1].im*g->f[1][1].im
                 +u->v[c][2].re*g->f[1][2].re-u->v[c][2].im*g->f[1][2].im;
    h->f[1][c].im=u->v[c][0].im*g->f[1][0].re+u->v[c][0].re*g->f[1][0].im
                 +u->v[c][1].im*g->f[1][1].re+u->v[c][1].re*g->f[1][1].im
                 +u->v[c][2].im*g->f[1][2].re+u->v[c][2].re*g->f[1][2].im;
}
#line 3542 "dwf.nw"
}
#line 3740 "dwf.nw"
  
#line 3745 "dwf.nw"
rs = &rx5[s];
for (c = 0; c < Nc; c++) {
    k = 6; 
#line 267 "dwf.nw"
qs->f[0][c].re = hh[k].f[0][c].re; qs->f[2][c].re = hh[k].f[0][c].re;
qs->f[0][c].im = hh[k].f[0][c].im; qs->f[2][c].im = hh[k].f[0][c].im;
qs->f[1][c].re = hh[k].f[1][c].re; qs->f[3][c].re = hh[k].f[1][c].re;
qs->f[1][c].im = hh[k].f[1][c].im; qs->f[3][c].im = hh[k].f[1][c].im;
#line 3748 "dwf.nw"
    k = 7; 
#line 280 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].im -= hh[k].f[1][c].im;
#line 3749 "dwf.nw"
    k = 2; 
#line 185 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im += hh[k].f[1][c].im;
#line 3750 "dwf.nw"
    k = 3; 
#line 198 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].re += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].im += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].re -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].im -= hh[k].f[1][c].im;
#line 3751 "dwf.nw"
    k = 1; 
#line 157 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re -= hh[k].f[1][c].im;
#line 3752 "dwf.nw"
    k = 0; 
#line 144 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[3][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[3][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[2][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[2][c].re += hh[k].f[1][c].im;
#line 3753 "dwf.nw"
    k = 5; 
#line 239 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im += hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re -= hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im -= hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re += hh[k].f[1][c].im;
#line 3754 "dwf.nw"
    k = 4; 
#line 226 "dwf.nw"
qs->f[0][c].re += hh[k].f[0][c].re; qs->f[2][c].im -= hh[k].f[0][c].re;
qs->f[0][c].im += hh[k].f[0][c].im; qs->f[2][c].re += hh[k].f[0][c].im;
qs->f[1][c].re += hh[k].f[1][c].re; qs->f[3][c].im += hh[k].f[1][c].re;
qs->f[1][c].im += hh[k].f[1][c].im; qs->f[3][c].re -= hh[k].f[1][c].im;
#line 3755 "dwf.nw"
}
#line 3741 "dwf.nw"
}
#line 3872 "dwf.nw"
    
#line 3917 "dwf.nw"
for (s = 0, vbc = vcbn, es1 = &ex5[Sv_1]; s < Sv; s++, vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_upN(es1->f[d][c].r, es->f[d][c].r)
    QXX(0,0,re); QXX(0,0,im);
    QXX(0,1,re); QXX(0,1,im);
    QXX(0,2,re); QXX(0,2,im);
    QXX(1,0,re); QXX(1,0,im);
    QXX(1,1,re); QXX(1,1,im);
    QXX(1,2,re); QXX(1,2,im);
#undef QXX
    es1 = es;
}
#line 3900 "dwf.nw"
for (s = Sv, vbc = vbnc, es1 = &ex5[0]; s--; vbc = vb) {
    es = &ex5[s];
    rs = &rx5[s];

#define QXX(d,c,r) rs->f[d][c].r += va * es->f[d][c].r \
                     + vbc * shift_up1(es->f[d][c].r, es1->f[d][c].r)
    QXX(2,0,re); QXX(2,0,im);
    QXX(2,1,re); QXX(2,1,im);
    QXX(2,2,re); QXX(2,2,im);
    QXX(3,0,re); QXX(3,0,im);
    QXX(3,1,re); QXX(3,1,im);
    QXX(3,2,re); QXX(3,2,im);
#undef QXX
    es1 = es;
}
#line 3873 "dwf.nw"
}
#line 3249 "dwf.nw"
    
#line 3263 "dwf.nw"
#undef qs
#undef qx5
#line 3250 "dwf.nw"
}
#line 592 "dwf.nw"
int
L3(DWF_init)(const int lattice[DIM+1],
             void *(*allocator)(size_t size),
             void (*deallocator)(void *))
{
    
#line 31 "dwf.nw"
static const char *version = "Version 1.3.0";

#line 599 "dwf.nw"
    if (inited_p)
        return 1; /* error: second init */

    
#line 648 "dwf.nw"
if (lattice[DIM] % Vs)
    goto error;
tlattice[DIM] = lattice[DIM];
#line 655 "dwf.nw"
{
    int i;
    for (i = 0; i < DIM; i++) {
        if (lattice[i] & 1)
            goto error;
        tlattice[i] = lattice[i];
    }
}
#line 603 "dwf.nw"
    
#line 1174 "dwf.nw"
{
    int i, dn;
    const int *xn, *xc;
    
    if (!QMP_logical_topology_is_declared())
        /* The user must have declared logical topology before */
        goto error;
    dn = QMP_get_logical_number_of_dimensions();
    if (dn > DIM)
        /* Too high dimension of the logical network */
        goto error;
        
    xn = QMP_get_logical_dimensions();
    xc = QMP_get_logical_coordinates();
    for (i = 0; i < dn; i++) {
        network[i] = xn[i];
        coord[i] = xc[i];
    }

    for (; i < DIM; i++) {
        network[i] = 1;
        coord[i] = 0;
    }
}
#line 604 "dwf.nw"
    
#line 625 "dwf.nw"
if (allocator)
    tmalloc = allocator;
else
    tmalloc = malloc;
	
if (deallocator)
    tfree = deallocator;
else
    tfree = free;
#line 605 "dwf.nw"
    
#line 1389 "dwf.nw"
if (init_tables()) {
    /* Something went wrong in the table construction */
    goto error;
}
#line 606 "dwf.nw"
    
#line 2196 "dwf.nw"
Phi_o  = allocate_odd_fermion();  if (Phi_o == 0) goto error;
auxA_o = allocate_odd_fermion();  if (auxA_o == 0) goto error;
auxB_o = allocate_odd_fermion();  if (auxB_o == 0) goto error;
auxA_e = allocate_even_fermion(); if (auxA_e == 0) goto error;
#line 2273 "dwf.nw"
r_o = allocate_odd_fermion(); if (r_o == 0) goto error;
p_o = allocate_odd_fermion(); if (p_o == 0) goto error;
q_o = allocate_odd_fermion(); if (q_o == 0) goto error;
#line 2294 "dwf.nw"
auxB_e = allocate_even_fermion();  if (auxB_e == 0) goto error;
#line 607 "dwf.nw"
    
#line 1943 "dwf.nw"
if (build_buffers(&even_odd)) goto error;
if (build_buffers(&odd_even)) goto error;
#line 608 "dwf.nw"
    
#line 4097 "dwf.nw"
{
   if (QMP_get_node_number() == 0) {
      QMP_printf("DWF init: %s (" MACHINE ")\n", version);
   }
}
#line 609 "dwf.nw"
    inited_p = 1;
    DEBUG_DWF("finished init, lattice=[%d %d %d %d %d]\n",
               lattice[0], lattice[1], lattice[2], lattice[3], lattice[4])
    return 0;
	
    
#line 619 "dwf.nw"
error:
   L3(DWF_fini)();
   return 1;
#line 615 "dwf.nw"
}
#line 671 "dwf.nw"
void
L3(DWF_fini)(void)
{
    
#line 2130 "dwf.nw"
free_buffers(&even_odd);
free_buffers(&odd_even);
#line 675 "dwf.nw"
    
#line 2202 "dwf.nw"
if (auxA_e) free16(auxA_e); auxA_e = 0;
if (auxB_o) free16(auxB_o); auxB_o = 0;
if (auxA_o) free16(auxA_o); auxA_o = 0;
if (Phi_o)  free16(Phi_o);  Phi_o = 0;
#line 2278 "dwf.nw"
if (r_o) free16(r_o); r_o = 0;
if (p_o) free16(p_o); p_o = 0;
if (q_o) free16(q_o); q_o = 0;
#line 2297 "dwf.nw"
if (auxB_e) free16(auxB_e); auxB_e = 0;
#line 676 "dwf.nw"
    
#line 1831 "dwf.nw"
{
    int i;

    if (neighbor.site) {
        tfree(neighbor.site);
        neighbor.site = 0;
    }

    if (neighbor.inside) {
       tfree(neighbor.inside);
       neighbor.inside = 0;
    }

    if (neighbor.boundary) {
       tfree(neighbor.boundary);
       neighbor.boundary = 0;
    }

    for (i = 2 * DIM; i--;) {
        if (neighbor.snd[i] == 0)
            continue;
        tfree(neighbor.snd[i]);
        neighbor.snd[i] = 0;
    }
}
#line 677 "dwf.nw"
    inited_p = 0;
    DEBUG_DWF("fini done\n")
}
#line 689 "dwf.nw"
L3(DWF_Fermion) *
L3(DWF_allocate_fermion)(void)
{
    L3(DWF_Fermion) *ptr;
    
    if (!inited_p)
        return 0;
	
    ptr = tmalloc(sizeof (*ptr));
    if (ptr == 0)
        return 0;
        
    ptr->even = allocate_even_fermion();
    if (ptr->even == 0)
        goto error1;
        
    ptr->odd  = allocate_odd_fermion();
    if (ptr->odd == 0)
        goto error2;
        
    return ptr;
  error2:
    free16(ptr->even);
  error1:
    tfree(ptr);
    return 0;
}
#line 722 "dwf.nw"
L3(DWF_Fermion) *
L3(DWF_load_fermion)(const void *OuterFermion,
                     void *env,
                     L3(DWF_fermion_reader) reader)
{
    L3(DWF_Fermion) *ptr = L3(DWF_allocate_fermion)();

    DEBUG_DWF("OuterFermion=%p, reader=%p\n", OuterFermion, reader)

    /* Handle both lack of memory and missing initialization */
    if (ptr == 0)
        return 0;

    
#line 1290 "dwf.nw"
{
    int x[DIM+1], i;
    
    
#line 1235 "dwf.nw"
for (i = 0; i < DIM; i++)
    x[i] = bounds.lo[i];
for (i = 0; i < DIM;) {
#line 1294 "dwf.nw"
        
#line 1301 "dwf.nw"
{
    int p = parity(x);
    int p1 = Sv * to_HFlinear(x, &bounds, -1, 0); /* p is taken care of! */
    vFermion *f = p? &ptr->odd[p1].f: &ptr->even[p1].f;

    for (x[DIM] = 0; x[DIM] < tlattice[DIM]; x[DIM] += Vs, f++) {
        int d;
        for (d = 0; d < Fd; d++) {
            int c;
            for (c = 0; c < Nc; c++) {
                f->f[d][c].re = import_vector(OuterFermion, env, reader,
                                              x, c, d, 0);
                f->f[d][c].im = import_vector(OuterFermion, env, reader,
                                              x, c, d, 1);                
            }
        }
    }
}
#line 1295 "dwf.nw"
    
#line 1241 "dwf.nw"
    for (i = 0; i < DIM; i++) {
        
#line 1258 "dwf.nw"
if (++x[i] == bounds.hi[i])
    x[i] = bounds.lo[i];
else
    break;
#line 1243 "dwf.nw"
    }
}
#line 1296 "dwf.nw"
}

#line 737 "dwf.nw"
    return ptr;
}
#line 743 "dwf.nw"
void
L3(DWF_save_fermion)(void *OuterFermion,
                     void *env,
                     L3(DWF_fermion_writer) writer,
                     L3(DWF_Fermion) *CGfermion)
{
    if (!inited_p)
        return;

    DEBUG_DWF("Outer fermion=%p, writer=%p\n", OuterFermion, writer)

    
#line 1340 "dwf.nw"
{
    int x[DIM+1], i;
    
    
#line 1235 "dwf.nw"
for (i = 0; i < DIM; i++)
    x[i] = bounds.lo[i];
for (i = 0; i < DIM;) {
#line 1344 "dwf.nw"
        
#line 1349 "dwf.nw"
{
    int p = parity(x);
    int p1 = Sv * to_HFlinear(x, &bounds, -1, 0); /* p is taken care of! */
    vFermion *f = p? &CGfermion->odd[p1].f: &CGfermion->even[p1].f;

    for (x[DIM] = 0; x[DIM] < tlattice[DIM]; x[DIM] += Vs, f++) {
        int d;
        for (d = 0; d < Fd; d++) {
            int c;
            for (c = 0; c < Nc; c++) {
                save_vector(OuterFermion, env, writer, x, c, d, 0,
                            &f->f[d][c].re);
                save_vector(OuterFermion, env, writer, x, c, d, 1,
                            &f->f[d][c].im);
            }
        }
    }
}
#line 1345 "dwf.nw"
    
#line 1241 "dwf.nw"
    for (i = 0; i < DIM; i++) {
        
#line 1258 "dwf.nw"
if (++x[i] == bounds.hi[i])
    x[i] = bounds.lo[i];
else
    break;
#line 1243 "dwf.nw"
    }
}
#line 1346 "dwf.nw"
}
#line 755 "dwf.nw"
}
#line 761 "dwf.nw"
void
L3(DWF_delete_fermion)(L3(DWF_Fermion) *ptr)
{
    if (!inited_p)
        return;
	
    free16(ptr->even);
    free16(ptr->odd);
    tfree(ptr);
}
#line 777 "dwf.nw"
L3(DWF_Gauge) *
L3(DWF_load_gauge)(const void *OuterGauge_U,
                   const void *OuterGauge_V,
                   void *env,
                   L3(DWF_gauge_reader) reader)
{
    L3(DWF_Gauge) *g;
    
    if (!inited_p)
        return 0;

    DEBUG_DWF("U=%p, V=%p, reader=%p\n", OuterGauge_U, OuterGauge_V, reader)
        
    g = allocate_gauge_field();
    if (g == 0)
        return 0;
 
    
#line 1210 "dwf.nw"
{
    int x[DIM], i, d, a, b, p1;
    
    
#line 1235 "dwf.nw"
for (i = 0; i < DIM; i++)
    x[i] = bounds.lo[i];
for (i = 0; i < DIM;) {
#line 1214 "dwf.nw"
        
#line 1223 "dwf.nw"
p1 = to_Ulinear(x, &bounds, -1);
for (d = 0; d < DIM; d++) {
    for (a = 0; a < Nc; a++) {
        for (b = 0; b < Nc; b++) {
            g[p1 + d].v[a][b].re = reader(OuterGauge_U, env, x, d, a, b, 0);
            g[p1 + d].v[a][b].im = reader(OuterGauge_U, env, x, d, a, b, 1);
        }
    }
}
#line 1215 "dwf.nw"
    
#line 1241 "dwf.nw"
    for (i = 0; i < DIM; i++) {
        
#line 1258 "dwf.nw"
if (++x[i] == bounds.hi[i])
    x[i] = bounds.lo[i];
else
    break;
#line 1243 "dwf.nw"
    }
}

#line 1217 "dwf.nw"
    for (d = 0; d < DIM; d++)
        
#line 1266 "dwf.nw"
{
    if (network[d] == 1)
        continue;

    
#line 1235 "dwf.nw"
for (i = 0; i < DIM; i++)
    x[i] = bounds.lo[i];
for (i = 0; i < DIM;) {
#line 1271 "dwf.nw"
        
#line 1277 "dwf.nw"
p1 = to_Ulinear(x, &bounds, d);
for (a = 0; a < Nc; a++) {
    for (b = 0; b < Nc; b++) {
        g[p1].v[a][b].re = reader(OuterGauge_V, env, x, d, a, b, 0);
        g[p1].v[a][b].im = reader(OuterGauge_V, env, x, d, a, b, 1);
    }
}
#line 1272 "dwf.nw"
    
#line 1249 "dwf.nw"
    for (i = 0; i < DIM; i++) {
        if (i == d)
            continue;
        
#line 1258 "dwf.nw"
if (++x[i] == bounds.hi[i])
    x[i] = bounds.lo[i];
else
    break;
#line 1253 "dwf.nw"
    }
}
#line 1273 "dwf.nw"
}
#line 1219 "dwf.nw"
}
#line 795 "dwf.nw"
    return g;
}
#line 809 "dwf.nw"
void
L3(DWF_delete_gauge)(L3(DWF_Gauge) *ptr)
{
    if (!inited_p)
        return;

    free16(ptr);
}
#line 822 "dwf.nw"
int
L3(DWF_cg_solver)(L3(DWF_Fermion)         *psi,
                  double                  *out_eps,
                  int                     *out_iter,
                  const L3(DWF_Gauge)     *gauge,
                  double                   M,
                  double                   m_f,
                  const L3(DWF_Fermion)   *x0,
                  const L3(DWF_Fermion)   *eta,
                  double                   eps,
                  int                      min_iter,
                  int                      max_iter)
{
    int status;

    if (!inited_p)
        return 1;

    U = (SU3 *)gauge;    
    
#line 4052 "dwf.nw"
{
    double a = M;
    double b = 2.;
    double c = -2*m_f;

    
#line 2658 "dwf.nw"
inv_a = 1.0 / a;
b_over_a = -b * inv_a;
#line 2735 "dwf.nw"
c0 = 1./(1.-c/b*d_pow(b/a, Sv*Vs));
vab = vmk_1(d_pow(b/a, Vs));
vfx_A = vmk_fn(-c0*c/a, -b/a);
vfx_B = vmk_bn(-b/a, -c0*c/a);
#line 4058 "dwf.nw"
}
#line 842 "dwf.nw"
    
#line 2185 "dwf.nw"
compute_Qee1(auxA_e, eta->even);
compute_Qoe(auxB_o, auxA_e);
compute_sum_o(auxA_o, eta->odd, -1, auxB_o);
compute_Qoo1(auxB_o, auxA_o);
compute_Mx(Phi_o, auxB_o);
#line 843 "dwf.nw"
    
#line 2211 "dwf.nw"
status = cg(psi->odd, Phi_o, x0->odd, eps, min_iter, max_iter, out_eps, out_iter);
#line 844 "dwf.nw"
    
#line 2286 "dwf.nw"
compute_Qeo(auxA_e, psi->odd);
compute_sum_e(auxB_e, eta->even, -1, auxA_e);
compute_Qee1(psi->even, auxB_e);
#line 845 "dwf.nw"
    return status;
}
#line 860 "dwf.nw"
void
L3(DWF_Dirac_Operator)(L3(DWF_Fermion)         *chi,
                       const L3(DWF_Gauge)     *gauge,
                       double                   M_0,
                       double                   m_f,
                       const L3(DWF_Fermion)   *psi)
{
    if (!inited_p)
        return;

    U = (SU3 *)gauge;
    
#line 500 "dwf.nw"
{
   double a = M_0;
   double b = 2.0;
   double c = -2 * m_f;

   va   = vmk_1(a);
   vb   = vmk_1(b);
   vcbn = vmk_1n(c, b);  /* {c, b, ...} */
   vbnc = vmk_n1(b, c);  /* {b, ..., c} */
}
#line 872 "dwf.nw"
    compute_Do(chi->odd,  psi->odd,  psi->even);
    compute_De(chi->even, psi->even, psi->odd);
}
#line 880 "dwf.nw"
void
L3(DWF_Dirac_Operator_conjugate)(L3(DWF_Fermion)        *chi,
                                 const L3(DWF_Gauge)    *gauge,
                                 double                  M_0,
                                 double                  m_f,
                                 const L3(DWF_Fermion)  *psi)
{
    if (!inited_p)
        return;

    U = (SU3 *)gauge;
    
#line 500 "dwf.nw"
{
   double a = M_0;
   double b = 2.0;
   double c = -2 * m_f;

   va   = vmk_1(a);
   vb   = vmk_1(b);
   vcbn = vmk_1n(c, b);  /* {c, b, ...} */
   vbnc = vmk_n1(b, c);  /* {b, ..., c} */
}
#line 892 "dwf.nw"
    compute_Dco(chi->odd, psi->odd, psi->even);
    compute_Dce(chi->even, psi->even, psi->odd);
}
#line 901 "dwf.nw"
void
L3(DWF_Add_Fermion)(L3(DWF_Fermion)         *psi,
                    const L3(DWF_Fermion)   *phi,
                    double                   a,
                    const L3(DWF_Fermion)   *eta)
{
    collect_add(&psi->odd->f,  &phi->odd->f,  a, &eta->odd->f,  odd_even.size * Sv);
    collect_add(&psi->even->f, &phi->even->f, a, &eta->even->f, even_odd.size * Sv);
}
#line 958 "dwf.nw"
void
L3(DWF_Fermion_Dot_Product)(double                  *r_re,
	                    double                  *r_im,
                            const L3(DWF_Fermion)   *psi,
                            const L3(DWF_Fermion)   *phi)
{
    *r_re = *r_im = 0;
    collect_dot(r_re, r_im, &psi->odd->f,  &phi->odd->f,  odd_even.size * Sv);
    collect_dot(r_re, r_im, &psi->even->f, &phi->even->f, even_odd.size * Sv);
    QMP_sum_double(r_re);
    QMP_sum_double(r_im);
}
