// programme principal
#include <stdio.h>
#include <stdlib.h>

#include "check.h"
#include "type.h"
#include "io.h"
#include "darboux.h"

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

  // COMPUTE
  
  d = darboux(m);

  // WRITE OUTPUT
  
      FILE *out;
    if(argc == 3)
      out = fopen(argv[2], "w");
    else
      out = stdout;
    mnt_write(d, out);
    if(argc == 3)
      fclose(out);
    else
      mnt_write_lakes(m, d, stdout);
  
  

  // free
  free(m->terrain);
  free(m);
  free(d->terrain);
  free(d);
  MPI_Finalize();
  return(0);
}
