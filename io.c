// fonctions d'entrée/sortie

#include <stdio.h>
#include "check.h"
#include "type.h"
#include "io.h"



void printTerrain (mnt *m, int r){
  for (int i = 0 ; i< m->nrows*m->ncols ; i++){
    printf(" %f ",m->terrain[i]);
  }
  
    printf("\n***************************************** Je suis le p %d\n",r);

  printf("\n");
}

mnt *mnt_read(char *fname)
{
  mnt *m;
  CHECK((m = malloc(sizeof(*m))) != NULL);
  int taille_chunk= 0;
  int reste = 0;
  int rank ; 
  int nbproc = 0; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
  
 
  
  if(rank==0){
    FILE *f;
    CHECK((f = fopen(fname, "r")) != NULL);
    //if rank == 0 do scans 
    CHECK(fscanf(f, "%d", &m->ncols) == 1);
    CHECK(fscanf(f, "%d", &m->nrows) == 1);
    CHECK(fscanf(f, "%f", &m->xllcorner) == 1);
    CHECK(fscanf(f, "%f", &m->yllcorner) == 1);
    CHECK(fscanf(f, "%f", &m->cellsize) == 1);
    CHECK(fscanf(f, "%f", &m->no_data) == 1);
    printf("Nombre de processeurs : %d \n",nbproc);

    /* calcul ici */
  
  
    CHECK((m->terrain = malloc(m->ncols * m->nrows * sizeof(float))) != NULL);
    
    taille_chunk = m->ncols  / nbproc ;
    reste = m->ncols  % nbproc;
    printf("Taille Calculee %i : et le reste : %i  \n",taille_chunk,reste);
    
    printf("Taille Chunk : %d \n , Taille Chunk Reste : %d \n",taille_chunk,reste);

    for(int i = 0 ; i < m->ncols * m->nrows ; i++)
    {
      CHECK(fscanf(f, "%f", &m->terrain[i]) == 1);
    }
    CHECK(fclose(f) == 0);
  }
    // Transfert de la valeur du chunk aux autres, ainsi que no_data et nrows
  MPI_Bcast(&taille_chunk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&reste , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
  MPI_Bcast(&m->no_data , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
  MPI_Bcast(&m->nrows , 1 , MPI_INT , 0 , MPI_COMM_WORLD);

  
  //Définition du nombre de lignes à traiter par le processus courant, autre que le processus 0
  if (rank){
    m->ncols = taille_chunk;
    CHECK((m->terrain = malloc(m->ncols * m->nrows * sizeof(float))) != NULL);
    printf("terrain ok dans pcs %i \n", rank);
  }
  
//Envoi et réception des valeurs du terrain aux processus
  if (!rank) {
    float *temp = m->terrain;
    for (int i = 1; i<nbproc; i++) {
      int count_send = taille_chunk * m->nrows;
      count_send += (i <= reste ? m->nrows : 0);
      printf("count send is %i for proc %i \n",count_send,i);
      
      MPI_Send(temp, count_send, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
      temp=(temp + count_send+1);
    }
    sleep(2);
  }
  else{
    //Les processus de rang reste ou moins traiteront une ligne de plus, si nbRows % nbProc != 0
    taille_chunk += (reste >= rank ? 1 : 0); 
    printf("taille_chunk = %i, nbrows = %i \n", taille_chunk, m->nrows);
    MPI_Recv(m->terrain , taille_chunk *  m->nrows , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , NULL);
  }
  for(int i = 0 ; i < nbproc ; i++){
    if(rank == i){
      printTerrain(m,rank);
    }
    else{
      sleep(2);
    }
    
  }
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
