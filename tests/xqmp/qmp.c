#include <stdarg.h>
#include <stdio.h>

#include <qmp.h>

extern int lattice_size[5];
extern int network[4];
extern int node[4];

int
QMP_get_node_number (void)
{
    int i;
    int n;

    for (i = 1, n = node[0]; i < 4; i++)
	n = n * network[i-1] + node[i];

    return n;
}

int
QMP_get_logical_number_of_dimensions (void)
{
    return 4;
}

const int*
QMP_get_logical_coordinates (void)
{
    return node;
}

const int*
QMP_get_logical_dimensions (void)
{
    return network;
}

/* rest is fake to make the linker happy */

QMP_status_t
QMP_init_msg_passing (int* argc, char*** argv,
		      QMP_thread_level_t required,
		      QMP_thread_level_t *provided)
{
    return QMP_SUCCESS;
}

QMP_mem_t *
QMP_allocate_aligned_memory (size_t nbytes,
			     size_t alignment,
			     int flags)
{
    static int count = 1;

    printf("alloc(%ld) -> mem[%2d] ", (long)nbytes, count);
    return (QMP_mem_t *)count++;
}

static int msg_count = 1;

QMP_msgmem_t
QMP_declare_msgmem (const void* mem, size_t nbytes)
{
/*     printf("msgmem(mem[%2d]) -> msgmem[%2d]\n", (int)mem, msg_count); */
    return (QMP_msgmem_t)msg_count++;
}

static int mh_count = 1;

QMP_msghandle_t   
QMP_declare_multiple (QMP_msghandle_t msgh[], 
		      int num)
{
    int i;

    printf("multiple(%d:", num);
    for (i = 0; i < num; i++)
	printf(" mh[%2d]", (int)msgh[i]);
    printf(") -> mh[%2d]\n", mh_count);
    return (QMP_msghandle_t)mh_count++;
}

QMP_msghandle_t    
QMP_declare_receive_relative (QMP_msgmem_t m, 
			      int axis,
			      int dir,
			      int priority)
{
/*
    printf("rcv(msgmem[%2d], dir[%d], ud[%2d]) -> mh[%2d]\n",
	   (int)m, axis, dir, mh_count);
*/
    printf("rcv(dir[%d], ud[%2d]) -> mh[%2d]\n", axis, dir, mh_count);
    return (QMP_msghandle_t)mh_count++;
}

QMP_msghandle_t    
QMP_declare_send_relative (QMP_msgmem_t m, 
			      int axis,
			      int dir,
			      int priority)
{
/*
    printf("snd(msgmem[%2d], dir[%d], ud[%2d]) -> mh[%2d]\n",
	   (int)m, axis, dir, mh_count);
*/
    printf("snd(dir[%d], ud[%2d]) -> mh[%2d]\n", axis, dir, mh_count);
    return (QMP_msghandle_t)mh_count++;
}

void
QMP_free_memory (QMP_mem_t* mem)
{
    printf("free_mem(%2d)\n", (int)mem);
}

void
QMP_free_msgmem (QMP_msgmem_t m)
{
}

void
QMP_free_msghandle (QMP_msghandle_t h)
{
}

void*
QMP_get_memory_pointer (QMP_mem_t* mem)
{
    return mem;
}

QMP_bool_t
QMP_logical_topology_is_declared (void)
{
    return 1;
}

QMP_status_t
QMP_start (QMP_msghandle_t h)
{
    return 0;
}

QMP_status_t
QMP_wait (QMP_msghandle_t h)
{
    return 0;
}

QMP_status_t
QMP_sum_double(double *value)
{
    return 0;
}

int
QMP_printf(const char *format, ...)
{
    va_list v;
    int status;

    va_start(v, format);
    status = vprintf(format, v);
    va_end(v);
    return status;
}
