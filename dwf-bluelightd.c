#line 4283 "dwf.nw"
#define Nc  3      /* Number of colors */
#define DIM 4      /* number of dimensions */
#define Fd  4      /* Fermion representation dimension */
#line 4872 "dwf.nw"
#include <bluelight.h>
typedef double REAL;
typedef vector double vReal;
#define Vs 2
#line 4295 "dwf.nw"
typedef struct {
    float re, im;
} scalar_complex;

typedef struct {
   vReal re, im;
} vector_complex;
#line 4767 "dwf.nw"
static inline vReal
vmk_1(double a)                                        /* return (a a ... a) */
{
     return (vReal)(vec_mk2d(a, a));
}
#line 4776 "dwf.nw"
static inline vReal
vmk_n1(double a, double b)                             /* return (a ... a b) */
{
  return (vReal)(vec_mk2d(a, b));
}
#line 4785 "dwf.nw"
static inline vReal
vmk_1n(double a, double b)                             /* return (a b ... b) */
{
  return (vReal)(vec_mk2d(a, b));
}
#line 4795 "dwf.nw"
static inline vReal
vmk_fn(double a, double b)                  /* return (a a*b ... a*b^(Vs-1)) */
{
  return (vReal)(vec_mk2d(a, a*b));
}
#line 4803 "dwf.nw"
static inline vReal
vmk_bn(double a, double b)                  /* return (a^(Vs-1)*b ... a*b b) */
{
  return (vReal)(vec_mk2d(a*b, b));
}
#line 4812 "dwf.nw"
static inline double
vsum(vReal v)                                         /* return sum(i, [i]v) */
{
   return vec_get0((vector double)v) + vec_get1((vector double)v);
}
#line 4821 "dwf.nw"
static inline vReal
vput_0(vReal a, double b)                     /* return (b [1]a ... [Vs-1]a) */
{
   return (vReal)(vec_mk2d(b, vec_get1((vector double)a)));
}
#line 4830 "dwf.nw"
static inline vReal
vput_n(vReal a, double b)                     /* return ([0]a ... [Vs-2]a b) */
{
   return (vReal)(vec_mk2d(vec_get0((vector double)a), b));
}
#line 4839 "dwf.nw"
static inline vReal
shift_up1(vReal a, vReal b)                /* return ([1]a ... [Vs-1]a [0]b) */
{
   return (vReal)(vec_mk10((vector double)a, (vector double)b));
}
#line 4847 "dwf.nw"
static inline vReal
shift_upN(vReal a, vReal b)             /* return ([Vs-1]a [0]b ... [Vs-2]b) */
{
   return (vReal)(vec_mk10((vector double)a, (vector double)b));
}
#line 4861 "dwf.nw"
#include "dwf-bluelightd.h"
#define MACHINE "bluelight double"
#define L3(n) MIT_bluelightd_##n
#define PAD16(size) (15+(size))
#define ALIGN16(addr) ((void *)(~15 & (15 + (size_t)(addr))))
#line 3028 "dwf.nw"
#define BLOCKOF_YA(j,t,c,ri) BLOCKOF2_YA(j,t,c,ri)
#define BLOCKOF_YB(j,t,c,ri) BLOCKOF2_YB(j,t,c,ri)
#line 4867 "dwf.nw"
#include "dwf.c"
