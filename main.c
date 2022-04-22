// programme principal
#include <stdio.h>
#include <stdlib.h>

#include "check.h"
#include "type.h"
#include "io.h"
#include "darboux.h"


//POur calculer le temps d'execution 
#define DIFFTEMPS(a,b) (((b).tv_sec - (a).tv_sec) + ((b).tv_usec - (a).tv_usec)/1000000.)

int main(int argc, char **argv)
{
  mnt *m, *d;

  if(argc < 2)
  {
    fprintf(stderr, "Usage: %s <input filename> [<output filename>]\n", argv[0]);
    exit(1);
  }
  
  //INIT MPI 
  MPI_Init(&argc,&argv);

  int rank ; 
  int nbproc = 0; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

  
  
  // READ INPUT
  m = mnt_read(argv[1]);


  struct timeval tv_init, tv_begin, tv_end;
  //Mesure de performances 
  if(!rank){
	  gettimeofday( &tv_init, NULL);
    gettimeofday( &tv_begin, NULL);
  }
  // COMPUTE
  d = darboux(m);

  if(!rank){
    gettimeofday( &tv_end, NULL);
   
  }
  
  // WRITE OUTPUT
  if (!rank){
      FILE *out;
    if(argc == 3)
      out = fopen(argv[2], "w");
    else{
      out = stdout;
      //mnt_write(d, out);
    }
    if(argc == 3)
      fclose(out);
    //else
      //mnt_write_lakes(m, d, stdout);

    
  }

  if(!rank){
     printf("Init : %lfs, Compute : %lfs\n",
          DIFFTEMPS(tv_init,tv_begin),
          DIFFTEMPS(tv_begin,tv_end));
  }
      
  // free
  free(m->terrain);
  free(m);
  
  if (!rank) {
    free(d->terrain);
    free(d);
    
  }
  MPI_Finalize();
   
  return(0);
}
