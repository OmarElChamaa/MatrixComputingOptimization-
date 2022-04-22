// fonctions d'entrée/sortie

#include <stdio.h>
#include "check.h"
#include "type.h"
#include "io.h"



void printTerrain (mnt *m, int r){
  printf("\n***************************************** Je suis le p %d\n",r);
  for (int i = 0 ; i< m->nrows*m->ncols ; i++){
    printf(" %f ",m->terrain[i]);
    if ((i+1)%m->ncols==0)
      printf("\n"); 
  }
  
  printf("\n***************************************** /Je suis le p %d\n",r);

  printf("\n");
}

//Décalage des valeurs de m d'une ligne, pour permettre la permission
void rowShift (mnt *m){
  for (int i = m->nrows * m->ncols; i>=0; i--) {
    m->terrain[i+m->ncols] = m->terrain[i];
  }
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
  
  if(!rank){
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
    
    taille_chunk = m->nrows  / (nbproc) ;
    reste = m->nrows  % (nbproc) ;  
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

  // m -> nrowsTemp sert de valeur temp pour stocker la taille total de la matrice de P0
  if(!rank)
    m->nrowsTemp = m->nrows;

  m->nrows = taille_chunk  + (rank < reste ? 1 : 0);
  if(rank)
    CHECK((m->terrain = malloc(m->ncols * (m->nrows + 1) * sizeof(float))) != NULL);
  //printf("%i rows dans p%i\n", m->nrows, rank);
  


  //Préparation du scatter des données de m->terrrain
  int* count_send = malloc(nbproc*sizeof(int));
  int* displacements = malloc(nbproc*sizeof(int));
  if (!rank) {
    for (int i = 0 ; i < nbproc; i++) {
      count_send[i] = taille_chunk * m->ncols;
      count_send[i] += (i < reste ? m->ncols : 0);
      displacements[i] = (i == 0 ? 0 : displacements[i-1] + count_send[i-1]);
      //printf("%i:%i, count send is %i \n",i,displacements[i],count_send[i]);
    } 

    MPI_Scatterv(m->terrain, count_send, displacements, MPI_FLOAT, MPI_IN_PLACE, m->nrows * m->ncols, MPI_FLOAT, 0, MPI_COMM_WORLD);
  } else {   
    MPI_Scatterv(NULL , NULL , NULL , NULL, m->terrain , m->nrows * m->ncols , MPI_FLOAT , 0 , MPI_COMM_WORLD);
  }
  //On décale les valeurs de m d'une ligne
  if (rank){
    rowShift(m); 
  }
   
  //On modifie le nombre de rows dans m pour le traitement de calcul_wij, comme on ajoute des lignes supplémentaires
  m->nrows += ((!rank || rank == nbproc - 1) ? 1 : 2);
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

  for(int i = 0 ; i < m->nrowsTemp ; i++)
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

  for(int i = 0 ; i < m->nrowsTemp ; i++)
  {
    for(int j = 0 ; j < m->ncols ; j++)
    {
      const float dif = TERRAIN(d,i,j)-TERRAIN(m,i,j);
      fprintf(f, "%c", (dif>1.)?'#':(dif>0.)?'+':'.' );
    }
    fprintf(f, "\n");
  }
}
