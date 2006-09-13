#include <qdp.h>

extern void test(void);

int
main(int argc, char *argv[])
{

    QDP_initialize(&argc, &argv);

    if (QDP_this_node == 0)
	test();

    QDP_finalize();

    return 0;
}
