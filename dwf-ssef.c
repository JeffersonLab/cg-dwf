#line 4215 "dwf.nw"
#define Nc  3      /* Number of colors */
#define DIM 4      /* number of dimensions */
#define Fd  4      /* Fermion representation dimension */
#line 4254 "dwf.nw"
#define Vs 4
typedef float REAL;
typedef REAL vReal __attribute__((mode(V4SF),aligned(16)));
#line 4227 "dwf.nw"
typedef struct {
    float re, im;
} scalar_complex;

typedef struct {
   vReal re, im;
} vector_complex;
#line 4263 "dwf.nw"
static inline vReal
vmk_1(double a)                                        /* return (a a ... a) */
{
     float b = a;
     vReal v = __builtin_ia32_loadss((float *)&b);
     asm("shufps\t$0,%0,%0" : "+x" (v));
     return v;
}
#line 4275 "dwf.nw"
static inline vReal
vmk_n1(double a, double b)                             /* return (a ... a b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a; r[2] = a; r[3] = b;
  return v;
}
#line 4287 "dwf.nw"
static inline vReal
vmk_1n(double a, double b)                             /* return (a b ... b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = b; r[2] = b; r[3] = b;
  return v;
}
#line 4300 "dwf.nw"
static inline vReal
vmk_fn(double a, double b)                  /* return (a a*b ... a*b^(Vs-1)) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a*b; r[2] = a*b*b; r[3] = a*b*b*b;
  return v;
}
#line 4311 "dwf.nw"
static inline vReal
vmk_bn(double a, double b)                  /* return (a^(Vs-1)*b ... a*b b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a*a*a*b; r[1] = a*a*b; r[2] = a*b; r[3] = b;
  return v;
}
#line 4323 "dwf.nw"
static inline double
vsum(vReal v)                                         /* return sum(i, [i]v) */
{
  REAL *r = (REAL *)&v;
  return r[0] + r[1] + r[2] + r[3];
}
#line 4333 "dwf.nw"
static inline vReal
vput_0(vReal a, double b)                     /* return (b [1]a ... [Vs-1]a) */
{
   REAL *v = (REAL *)&a;
   v[0] = b;
   return a;
}
#line 4344 "dwf.nw"
static inline vReal
vput_n(vReal a, double b)                     /* return ([0]a ... [Vs-2]a b) */
{
   REAL *v = (REAL *)&a;
   v[3] = b;
   return a;
}
#line 4356 "dwf.nw"
static inline vReal
shift_up1(vReal a, vReal b)                /* return ([1]a ... [Vs-1]a [0]b) */
{
   vReal x = a;
   vReal y = b;
   asm("shufps\t$0x30,%0,%1\n\t"
       "shufps\t$0x29,%1,%0"
       : "+x" (x), "+x" (y));
   return x;
}
#line 4369 "dwf.nw"
static inline vReal
shift_upN(vReal a, vReal b)             /* return ([Vs-1]a [0]b ... [Vs-2]b) */
{
   vReal x = a;
   asm("shufps\t$0x03,%1,%0\n\t"
       "shufps\t$0x9c,%1,%0"
       : "+x" (x): "x" (b));
   return x;
}
#line 4241 "dwf.nw"
#include "dwf-ssef.h"
#define MACHINE "sse float"
#define L3(n) MIT_ssef_##n
#define PAD16(size) (15+(size))
#define ALIGN16(addr) ((void *)(~15 & (15 + (size_t)(addr))))
#line 3037 "dwf.nw"
#define BLOCKOF_YA(j,t,c,ri) BLOCKOF4_YA(j,t,c,ri)
#define BLOCKOF_YB(j,t,c,ri) BLOCKOF4_YB(j,t,c,ri)
#line 4247 "dwf.nw"
#include "dwf.c"
