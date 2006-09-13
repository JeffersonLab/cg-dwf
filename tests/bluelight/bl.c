#include <bluelight.h>
#include <stdio.h>

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
dump_d(const char *name, double x)
{
    dump_x(name, &x, sizeof (x));
}

static void
dump_v(const char *name, vector double x)
{
    dump_x(name, &x, sizeof (x));
}

void
test(void)
{
    static double data[] = { 133./64, 295./128, 73./32, 987./256 };
    vector double v, w;

    dump_x("data", data, sizeof (data));
    dump_d("d[0]", data[0]);
    dump_d("d[1]", data[1]);
    dump_d("d[2]", data[2]);
    dump_d("d[3]", data[3]);
    
    v = vec_mk2d(data[0], data[1]);
    w = vec_mk2d(data[2], data[3]);

    dump_x("vxx", &v, sizeof (v));
    dump_v("v2d", v);
    dump_d("v/0", vec_get0(v));
    dump_d("v/1", vec_get1(v));
    
    dump_x("wxx", &w, sizeof (w));
    dump_v("w2d", w);
    dump_d("w/0", vec_get0(w));
    dump_d("w/1", vec_get1(w));

    dump_v("z00", vec_mk00(v, w));
    dump_v("z10", vec_mk10(v, w));
    dump_v("z11", vec_mk11(v, w));
}
