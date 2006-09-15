#include <stdio.h>

typedef float REAL;
typedef REAL vReal __attribute__((mode(V4SF),aligned(16)));

static inline vReal
vmk_xx(REAL a[])
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a[0]; r[1] = a[1]; r[2] = a[2]; r[3] = a[3];
  return v;
}

static inline vReal
shift_up1(vReal a, vReal b)                /* return ([1]a ... [Vs-1]a [0]b) */
{
   vReal z;
   REAL *X = (REAL *)&a;
   REAL *Y = (REAL *)&b;
   REAL *Z = (REAL *)&z;

   Z[0] = X[1]; Z[1] = X[2]; Z[2] = X[3]; Z[3] = Y[0];
   return z;
}

static inline vReal
shift_upN(vReal a, vReal b)             /* return ([Vs-1]a [0]b ... [Vs-2]b) */
{
   vReal z;
   REAL *X = (REAL *)&a;
   REAL *Y = (REAL *)&b;
   REAL *Z = (REAL *)&z;

   Z[0] = X[3]; Z[1] = Y[0]; Z[2] = Y[1]; Z[3] = Y[2];
   return z;
}


static void
dump_x(const char *name, void *v, int n)
{
    int i;
    printf("%s: ", name);
    for (i = 0; i < n; i++)
	printf("%02x", ((unsigned char *)v)[i]);
    printf("\n");
}

static void
dump_d(const char *name, REAL x) 
{
    dump_x(name, &x, sizeof (x));
}

static void
dump_v(const char *name, vReal x)
{
    dump_x(name, &x, sizeof (x));
}

void
test(void)
{
    static REAL data[] = { 1.45634575467356e+1,
			   2.25643563456345e+2,
			   3.54634574867247e+3,
			   4.74567546872445e+4,
			   5.09876542345665e+5,
			   6.32456345677741e+6,
			   7.65467456781374e+7,
			   8.03673487566586e+8 };
			     
    vReal v, w;

    dump_x("data", data, sizeof (data));
    dump_d("d[0]", data[0]);
    dump_d("d[1]", data[1]);
    dump_d("d[2]", data[2]);
    dump_d("d[3]", data[3]);
    dump_d("d[4]", data[4]);
    dump_d("d[5]", data[5]);
    dump_d("d[6]", data[6]);
    dump_d("d[7]", data[7]);
    
    v = vmk_xx(&data[0]);
    w = vmk_xx(&data[4]);

    dump_x("vxx ", &v, sizeof (v));
    dump_v("v2d ", v);
    
    dump_x("wxx ", &w, sizeof (w));
    dump_v("w2d ", w);

    dump_v("s1up", shift_up1(v, w));
    dump_v("s1dn", shift_upN(v, w));
}
