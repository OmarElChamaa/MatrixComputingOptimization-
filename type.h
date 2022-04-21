
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
  int nrowsTemp;            //POur savoir combie de rows wst attribuer au rang 0 et pas son nrows total             

  float *terrain;                     // linear array (size: ncols*nrows)
}
mnt;

// access to terrain in an mnt m as a 2D array:
#define TERRAIN(m,i,j) (m->terrain[(i)*m->ncols+(j)])

#endif
