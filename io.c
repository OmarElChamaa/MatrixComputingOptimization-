// fonctions d'entrée/sortie

#include <stdio.h>
#include "check.h"
#include "type.h"
#include "io.h"



void printTerrain (mnt *m, int r){
  for (int i = 0 ; i< m->nrows*m->ncols ; i++){
    printf(" %f ",m->terrain[i]);
    if (i%m->ncols==0 && i != 0)
      printf("\n"); 
  }
  
    printf("\n***************************************** Je suis le p %d\n",r);

  printf("\n");
}

mnt *mnt_read(char *fname)
{
  mnt *m;
  CHECK((m = malloc(sizeof(*m))) != NULL);
  int taille_chunk = 0;
  int reste = 0;
  int rank ; 
  int nbproc = 0; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
  
  if(rank==0){
    FILE *f;
    CHECK((f = fopen(fname, "r")) != NULL);
    CHECK(fscanf(f, "%d", &m->ncols) == 1);
    CHECK(fscanf(f, "%d", &m->nrows) == 1);
    CHECK(fscanf(f, "%f", &m->xllcorner) == 1);
    CHECK(fscanf(f, "%f", &m->yllcorner) == 1);
    CHECK(fscanf(f, "%f", &m->cellsize) == 1);
    CHECK(fscanf(f, "%f", &m->no_data) == 1);

    /* calcul ici */
  
  
    CHECK((m->terrain = malloc(m->ncols * m->nrows * sizeof(float))) != NULL);
    
    taille_chunk = m->nrows  / (nbproc-1) ;
    reste = m->nrows  % (nbproc-1) ;  
    for(int i = 0 ; i < m->ncols * m->nrows ; i++)
    {
      CHECK(fscanf(f, "%f", &m->terrain[i]) == 1);
    }
    CHECK(fclose(f) == 0);
  }
  // Transfert de la valeur du chunk aux autres, ainsi que no_data, ncols
  MPI_Bcast(&taille_chunk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&reste , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
  MPI_Bcast(&m->no_data , 1 , MPI_INT , 0 , MPI_COMM_WORLD);  
  MPI_Bcast(&m->ncols , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
  
  //Définition du nombre de lignes à traiter par le processus courant, autre que le processus 0
  if (rank){
    m->nrows = taille_chunk + (rank <= reste ? 1 : 0) ;
    CHECK((m->terrain = malloc(m->ncols * m->nrows * sizeof(float))) != NULL);
    //printf("%i rows dans p%i\n", m->nrows, rank);
  }


  //Préparation du scatter des données de m->terrrain
  int* count_send = malloc(nbproc*sizeof(int));
  int* displacements = malloc(nbproc*sizeof(int));
  count_send[0] = 0;
  displacements[0] = 0;
  if (!rank) {
    for (int i = 1 ; i < nbproc; i++) {
      count_send[i] = taille_chunk * m->ncols;
      count_send[i] += ((i) <= reste ? m->ncols : 0);
      displacements[i] = displacements[i-1] + ( i > 1 ? count_send[i-1] : 0);
    }
    MPI_Scatterv(m->terrain, count_send, displacements, MPI_FLOAT, NULL, 0, MPI_FLOAT, 0, MPI_COMM_WORLD);
  } else { 
    
    taille_chunk += (reste >= rank ? 1 : 0);
    MPI_Scatterv(NULL , NULL , NULL , MPI_FLOAT, m->terrain , taille_chunk * m->ncols , MPI_FLOAT , 0 , MPI_COMM_WORLD);
  }
  printTerrain(m, rank);
  return(m);
}

void mnt_write(mnt *m, FILE *f)
{
  CHECK(f != NULL);

  fprintf(f, "%d\n", m->ncols);
  fprintf(f, "%d\n", m->nrows);
  fprintf(f, "%.2f\n", m->xllcorner);
  fprintf(f, "%.2f\n", m->yllcorner);
  fprintf(f, "%.2f\n", m->cellsize);
  fprintf(f, "%.2f\n", m->no_data);

  for(int i = 0 ; i < m->nrows ; i++)
  {
    for(int j = 0 ; j < m->ncols ; j++)
    {
      fprintf(f, "%.2f ", TERRAIN(m,i,j));
    }
    fprintf(f, "\n");
  }
}

void mnt_write_lakes(mnt *m, mnt *d, FILE *f)
{
  CHECK(f != NULL);

  for(int i = 0 ; i < m->nrows ; i++)
  {
    for(int j = 0 ; j < m->ncols ; j++)
    {
      const float dif = TERRAIN(d,i,j)-TERRAIN(m,i,j);
      fprintf(f, "%c", (dif>1.)?'#':(dif>0.)?'+':'.' );
    }
    fprintf(f, "\n");
  }
}
