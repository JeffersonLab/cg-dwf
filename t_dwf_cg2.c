
#include <stdio.h>
#include <string.h>
#include<sys/time.h>
#include <time.h>
#include <errno.h>

#include "qmp.h"
#include "dwf-ssef.h"
#define L3(n) MIT_ssef_##n
#define ND 4
/* Forward Declarations */
void parse_options(int argc, char *argv[], int *length, int *maxCG);
void get_int_next(int i, int argc, char *argv[], int *ret_val);
void usage();


double
DWF_gauge_reader(const void *ptr, void *env, 
		 const int latt_coord[4], int mu, int row, int col, int reim)
{
  double val;

  if (row == col)
    val = (reim == 0) ? 1.0 : 0.0;
  else
    val = 0.0;

  return val;
}

double
DWF_fermion_reader1(const void *ptr, void *env, 
		    const int latt_coord[5], int color, int spin, int reim)
{
  return 0.0;  /* Often set to 1 for a non-trivial RHS */
}

double
DWF_fermion_reader2(const void *ptr, void *env, 
		    const int latt_coord[5], int color, int spin, int reim)
{
  return 0.0;
}

void
DWF_fermion_writer(void *ptr, void *env, 
		   const int latt_coord[5], int color, int spin, int reim,
		   double val)
{
  return;
}




static void
error(const char *msg)
{
    L3(DWF_fini)();
    QMP_fprintf(stderr, "FATAL ERROR: %s\n", msg);
    QMP_abort(1);
    QMP_finalize_msg_passing();
}


int main(int argc, char **argv)
{
  int length[ND+1];
  int i,err;
  L3(DWF_Gauge)*   g;
  L3(DWF_Fermion)* rhs;
  L3(DWF_Fermion)* x0;
  L3(DWF_Fermion)* x;
  float *u = NULL, *v=NULL, *psi=NULL, *chi=NULL;
  int vol, MaxCG;
  int subgrid[ND];
  QMP_thread_level_t prv;

  QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv);
  parse_options(argc, argv, length, &MaxCG);

  /* 
  if ( QMP_is_primary_node() ) { 
    fprintf(stderr, "LatticeSize: %d %d %d %d %d MaxCG: %d\n", 
	    length[0], length[1], length[2], length[3], length[4], MaxCG);
  }
  */

  /* Set up the grid */
  err = QMP_layout_grid(length, ND);
  {
    const int *s = QMP_get_subgrid_dimensions();
    vol = 1;
    for(i=0; i < ND; ++i)
    {
      vol *= length[i];
      subgrid[i] = s[i];
    }
  }

  err = L3(DWF_init)(length, NULL, NULL);
  if (err)
    error("error initializing SSE_DWF_init");


  /* Initialize the SSE specific fields */
  if ((g = L3(DWF_load_gauge)(u, v, NULL, DWF_gauge_reader)) == NULL)
    error("error initializing g");

  if ((rhs = L3(DWF_load_fermion)(chi, NULL, DWF_fermion_reader1)) == NULL)
    error("error initializing rhs");

  if ((x0 = L3(DWF_load_fermion)(psi, NULL, DWF_fermion_reader2)) == NULL)
    error("error initializing x0");

  if ((x = L3(DWF_allocate_fermion)()) == NULL)
    error("error initializing x");

  {
    /* Call the solver */
    double out_epsilon = 0.0;
    double M_0 = 1.0;    /* Warning check sign convention */
    double m_q = 0.01;
    double Rsd = 0;
    int    Ls = length[4];
    int iter = 0;
    clock_t myt1, myt2;
    double  mydt;
    int Nc = 3;
    int Ns = 4;

    unsigned long Ndiag  = (4*Ls+2)*Nc*Ns; /* This is my count with the blas / chiral proj ops */
    //    unsigned long NdiagInv = (10*Ls-8)*Nc*Ns;
    unsigned long Neo    = Ls*(1320+24);
    unsigned long N_mpsi = 2*Ndiag + 2*Neo + Ls*24;
    unsigned long Nflops_c = (24*Ls + 2*N_mpsi) + (48*Ls);  /* constant term */
    unsigned long Nflops_s = (2*N_mpsi + Ls*(2*48+2*24));   /* slope */
    unsigned long  Nflops_per_cbsite;
    double Nflops_total;
    double mflops;

    struct timeval t_start;
    struct timeval t_end;

   
    myt1=clock();
    gettimeofday(&t_start, NULL);
    err = L3(DWF_cg_solver_timings)(x, &out_epsilon, &iter,
				    g, M_0, m_q,
				    x0, rhs, 
				    Rsd, 0, MaxCG);
    gettimeofday(&t_end, NULL);
    myt2=clock();
    {
      long secs=0;
      long usecs=0;

      secs = t_end.tv_sec - t_start.tv_sec;
      if( t_end.tv_usec < t_start.tv_usec ) 
      {
	secs -= 1;
	usecs = 1000000;
      }
      usecs += t_end.tv_usec - t_start.tv_usec;
      mydt = (double)secs + ((double)usecs / 1e6);
    }

    if (err)
      error("error in cg solver");

    /* Flop count for inverter */
    Nflops_per_cbsite = Nflops_c + iter*Nflops_s;
    Nflops_total = Nflops_per_cbsite*((double)(vol/2));
    mflops = Nflops_total / mydt;
    mflops /= 1.0e6;

    if( QMP_is_primary_node() ) { 
      printf("Nproc: %d Global: %d %d %d %d Ls: %d Local: %d %d %d %d : iter: %d  Nflops: %g  dt: %g  floppage_total: %g flops_per_process: %g\n",
	     QMP_get_number_of_nodes(),
	     length[0], length[1], length[2], length[3], length[4],
	     subgrid[0],subgrid[1],subgrid[2],subgrid[3],
	     iter, Nflops_total, mydt, 
	     mflops, mflops/(double)QMP_get_number_of_nodes());
    }
  }

    /* Save the result */
  L3(DWF_save_fermion)(psi, NULL, DWF_fermion_writer, x);

  /* Cleanup */
  L3(DWF_delete_fermion)(x);
  L3(DWF_delete_fermion)(x0);
  L3(DWF_delete_fermion)(rhs);
  L3(DWF_delete_gauge)(g);

  /* Bolt */
  L3(DWF_fini)();
  QMP_finalize_msg_passing();

  return 0;
}


#define FALSE 0
#define TRUE 1
void parse_options(int argc, char *argv[], int *length, int *maxCG)
{
  int i; 
  unsigned char found_lx, found_ly, found_lz, found_lt, found_ls, found_maxcg, found_all;

  found_ls = found_lx = found_ly = found_lz = found_lt = found_maxcg = FALSE;
  
 

  for(i=0; i < argc; i++) { 
    
    if( strcmp(argv[i], "-Lx") == 0 ) {
      get_int_next(i, argc, argv, &length[0]);
      found_lx = TRUE;
      i++;
    }
    else if( strcmp(argv[i], "-Ly") == 0 ) { 
      get_int_next(i, argc,argv, &length[1]);
      found_ly = TRUE;
      i++;
    }
    else if( strcmp(argv[i], "-Lz") == 0 ) { 
      get_int_next(i, argc, argv, &length[2]);
      found_lz = TRUE;
      i++;
    }
    else if( strcmp(argv[i], "-Lt") == 0 ) { 
      get_int_next(i, argc, argv, &length[3]);
      found_lt = TRUE;
      i++;
    }
    else if( strcmp(argv[i], "-Ls") == 0 ) {
      get_int_next(i, argc, argv, &length[4]);
      found_ls = TRUE;
      if( (length[4] % 4) != 0 ) { 
	error("Ls must be divisible by 4\n");
      }
      i++;
    }
    else if( strcmp(argv[i], "-maxCG") == 0 ) {
      get_int_next(i, argc, argv, maxCG);
      found_maxcg = TRUE;
      i++;
    }
  }

  found_all = found_lx && found_ly && found_lz && found_lt && found_ls
    && found_maxcg;

  if ( ! found_all ) { 
    usage();
  }

}

void get_int_next(int i, int argc, char *argv[], int *ret_val)
{
  if ( i+1 < argc ) { 
    *ret_val = strtol( argv[i+1], (char **)NULL, 10 );
    if( errno == EINVAL ) { 
      usage();
    }
  }
  else { 
    usage();
  }
}

void usage()
{

  if( QMP_is_primary_node() == 1) { 
    fprintf(stderr, "Usage: t_dwf_cg2 <arguments>\n");
    fprintf(stderr, "  The following arguments are required: \n");
    fprintf(stderr, "     -Lx <integer>  global size in X direction \n");
    fprintf(stderr, "     -Ly <integer>  global size in Y direction \n");
    fprintf(stderr, "     -Lz <integer>  global size in Z direction \n");
    fprintf(stderr, "     -Lt <integer>  global size in T direction \n");
    fprintf(stderr, "     -Ls <integer>  global size in S direction \n");
    fprintf(stderr, "     -maxCG <integer>  maximum no of CG iterations \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Optionally one can add QMP arguments such as -qmp-geom\n");
    fprintf(stderr, "\n");
  }

  QMP_finalize_msg_passing();
  exit(1);

}


  
