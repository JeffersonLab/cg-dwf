#include <stdio.h>
#include <qmp.h>

#ifndef PREFIX
#define PREFIX MIT_altivecf_
#endif

#include "../../dwf-altivecf.h"
#include "../../dwf-bluelightf.h"
#include "../../dwf-bluelightd.h"

#define DIM 4
typedef double vHalfFermion;

#include "dwf-tables.h"

#define __L3(a,b)  a##b
#define _L3(a,b)   __L3(a,b)
#define L3(n)      _L3(PREFIX,n)

int lattice_size[5];
int network[4];
int node[4];

void
xxx(void)
{
    fflush(stdout);
}

/*****************************************************************************/
static void
dump_table(const char *name, struct neighbor *nb)
{
    int d, i;

    printf("\f\n");
    printf("%s: buffers:\n", name);
    for (i = 0; i < 2*DIM; i++) {
	printf("%s: send[%d] %08p, rcv[%d] %08p\n", name,
	       i, nb->snd_buf[i], i, nb->rcv_buf[i]);
    }
    printf("%s: qmp buffers\n", name);
    for (i = 0; i < 4 * DIM; i++) {
	printf("%s: %d qmp_size %5d, qmp_xbuf %08p qmp_buf %08p\n",
	       name, i, nb->qmp_size[i], nb->qmp_xbuf[i], nb->qmp_buf[i]);
    }
    printf("%s: boundary size %d\n", name, nb->boundary_size);
    for (i = 0; i < nb->boundary_size; i++) {
	int x = nb->boundary[i].index;
	int m = nb->boundary[i].mask;
	printf("%s: bnd[%3d] %02x %4d/ ", name, i, m, x);
	for (d = 0; d < 2*DIM; d++)
	    if (m & (1 << d))
		printf(" %3d:", nb->site[x].F[d]);
	    else
		printf("  %3d", nb->site[x].F[d]);
	printf("\n");
    }
    for (d = 0; d < 2*DIM; d++) {
	printf("\n%s: send size[%d] %d\n", name, d, nb->snd_size[d]);
	for (i = 0; i < nb->snd_size[d]; i++)
	    printf("%s: send[%d][%3d] %3d\n", name, d, i, nb->snd[d][i]);
    }
}

/*****************************************************************************/

int
main(int argc, char *argv[])
{
    QMP_thread_level_t qmp_unused_1;
    int                i;

    QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &qmp_unused_1);
    
    if (argc != 14) {
	printf("usage: %s Lx Ly Lz Lt Ls Nx Ny Nz Nt nx ny nz nt\n", argv[0]);
	return 1;
    }

    for (i = 0; i < 5; i++)
	lattice_size[i] = atoi(argv[1 + i]);

    for (i = 0; i < 4; i++)
	network[i] = atoi(argv[6 + i]);

    for (i = 0; i < 4; i++)
	node[i] = atoi(argv[10 + i]);

    printf("lattice %d %d %d %d %d\n",
	   lattice_size[0],
	   lattice_size[1],
	   lattice_size[2],
	   lattice_size[3],
	   lattice_size[4]);
    printf("network %d %d %d %d\n",
	   network[0],
	   network[1],
	   network[2],
	   network[3]);
    printf("node %d %d %d %d\n",
	   node[0],
	   node[1],
	   node[2],
	   node[3]);
    
    L3(DWF_init)(lattice_size, NULL, NULL);

    dump_table("even_odd", &even_odd);
    dump_table("odd_even", &odd_even);
    
    return 0;
}


