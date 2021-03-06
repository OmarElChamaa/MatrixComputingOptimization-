
#include <omp.h>
#include <mpi.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
// structure de données principale
#ifndef __TYPE_H__
#define __TYPE_H__

typedef struct mnt_t
{
  int ncols, nrows;                   // size
  float xllcorner, yllcorner, cellsize; // not used
  float no_data; // mnt value unknown
  int nrowsTemp;            //Permet de stocker la valeur initiale du nombre de ligne du mnt (pour P0)       

  float *terrain;                     // linear array (size: ncols*nrows)
}
mnt;

// access to terrain in an mnt m as a 2D array:
#define TERRAIN(m,i,j) (m->terrain[(i)*m->ncols+(j)])

#endif
