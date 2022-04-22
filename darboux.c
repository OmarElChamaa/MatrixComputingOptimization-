// fonction de calcul principale : algorithme de Darboux
// (remplissage des cuvettes d'un MNT)
#include <string.h>

#include "check.h"
#include "type.h"
#include "darboux.h"

// si ce define n'est pas commenté, l'exécution affiche sur stderr la hauteur
// courante en train d'être calculée (doit augmenter) et l'itération du calcul
#define DARBOUX_PPRINT

#define PRECISION_FLOTTANT 1.e-5

// pour accéder à un tableau de flotant linéarisé (ncols doit être défini) :
#define WTERRAIN(w,i,j) (w[(i)*ncols+(j)])

float max_terrain_kernel(const mnt *restrict m, int rank, int nbproc)
{
  float max = m->terrain[m->ncols];
  //Décalage selon le numéro du processus
  int stride = ((!rank || rank == nbproc -1) ? 1 : 2);
  int nrows = m->nrows - stride;
  for(int i = m->ncols ; i < m->ncols * nrows ; i++)
    if(m->terrain[i] > max)
      max = m->terrain[i];
  return(max);
}


// calcule la valeur max de hauteur sur un terrain
float max_terrain(const mnt *restrict m)
{
  float max = m->terrain[0];
  for(int i = 0 ; i < m->ncols * m->nrows ; i++)
    if(m->terrain[i] > max)
      max = m->terrain[i];
  return(max);
}



// initialise le tableau W de départ à partir d'un mnt
float *init_W(const mnt *restrict m)
{
  const int ncols = m->ncols, nrows = m->nrows;
  float *restrict W;
  CHECK((W = malloc(ncols * nrows * sizeof(float))) != NULL);

  // initialisation W
  const float max = max_terrain(m) + 10.;
  for(int i = 0 ; i < nrows ; i++)
  {
    for(int j = 0 ; j < ncols ; j++)
    {
      if(i==0 || i==nrows-1 || j==0 || j==ncols-1 || TERRAIN(m,i,j) == m->no_data)
        WTERRAIN(W,i,j) = TERRAIN(m,i,j);
      else
        WTERRAIN(W,i,j) = max;
    }
  }

  return(W);
}

//Initialise le tableau W de départ à partir d'un mnt, en ajoutant une ligne vide avant et après les données de m.
//La concaténation des W
float *init_W_kernel(const mnt *restrict m, float gmax){
  
  float *restrict W;
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int nbproc ; 
  MPI_Comm_size(MPI_COMM_WORLD, &nbproc); 

  const int ncols = m->ncols ; 
  const int nrows = m->nrows  ; 
  CHECK((W = calloc( ncols * nrows , sizeof(float))) != NULL);

  //traitement différent selon le processus, pour respecter la "bordure" de W initial
  if (rank == 0) {
      for(int i = 0 ; i < nrows - 1; i++)
      {
        for(int j = 0 ; j < ncols ; j++)
        {
          if(i==0 || j==0 ||  j==ncols-1 || TERRAIN(m,i,j) == m->no_data){
            WTERRAIN(W,i,j) = TERRAIN(m,i,j);
          } 
          else
            WTERRAIN(W,i,j) = gmax;
        }
      }
    }
    else if (rank == nbproc - 1) {
      for(int i = 1; i < nrows; i++)
      {
        for(int j = 0 ; j < ncols ; j++)
        {
          if(i == nrows-1|| j==0 ||  j==ncols-1 || TERRAIN(m,i,j) == m->no_data)
            WTERRAIN(W,i,j) = TERRAIN(m,i,j);
          else
            WTERRAIN(W,i,j) = gmax;
        }
      }
    }
    else {
     for(int i = 1; i < nrows-1; i++)
      {
        for(int j = 0 ; j < ncols ; j++)
        {
          if(j==0 || j==ncols-1 || TERRAIN(m,i,j) == m->no_data)
            WTERRAIN(W,i,j) = TERRAIN(m,i,j);
          else
            WTERRAIN(W,i,j) = gmax;
        }
      }
    }

  return(W);
}

// variables globales pour l'affichage de la progression
#ifdef DARBOUX_PPRINT
float min_darboux=9999.; // ça ira bien, c'est juste de l'affichage
int iter_darboux=0;
// fonction d'affichage de la progression
void dpprint()
{
  if(min_darboux != 9999.)
  {
    fprintf(stderr, "%.3f %d\r", min_darboux, iter_darboux++);
    fflush(stderr);
    min_darboux = 9999.;
  }
  else
    fprintf(stderr, "\n");
}
#endif


// pour parcourir les 8 voisins :
const int VOISINS[8][2] = {{-1,-1}, {-1,0}, {-1,1}, {0,-1}, {0,1}, {1,-1}, {1,0}, {1,1}};

// cette fonction calcule le nouveau W[i,j] en utilisant Wprec[i,j]
// et ses 8 cases voisines : Wprec[i +/- 1 , j +/- 1],
// ainsi que le MNT initial m en position [i,j]
// inutile de modifier cette fonction (elle est sensible...):
int calcul_Wij(float *restrict W, const float *restrict Wprec, const mnt *m, const int i, const int j)
{
  
  int rank ; 
  int nbproc; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

  const int nrows = m->nrows , ncols = m->ncols;
  int modif = 0;

  // on prend la valeur précédente...
  WTERRAIN(W,i,j) = WTERRAIN(Wprec,i,j);
  // ... sauf si :
  
  if(WTERRAIN(Wprec,i,j) > TERRAIN(m,i,j))
  {
    // parcourir les 8 voisins haut/bas + gauche/droite
    for(int v=0; v<8; v++)
    {
      const int n1 = i + VOISINS[v][0];
      const int n2 = j + VOISINS[v][1];
     // printf("i : %i j: %i v: %i\n",i,j,v);
      //printf("Rank %i Wterrain:%f , Terrain:%f\n",rank,WTERRAIN(Wprec, i, j), TERRAIN(m,i,j));
      //printf("proc %i val n1 %i n2 %i nrows %i ncols %i\n",rank , n1, n2,nrows,ncols);

      // vérifie qu'on ne sort pas de la grille.
      // ceci est théoriquement impossible, si les bords de la matrice Wprec
      // sont bien initialisés avec les valeurs des bords du mnt
    
      CHECK(n1>=0 && n1<nrows && n2>=0 && n2<ncols);

      // si le voisin est inconnu, on l'ignore et passe au suivant
      if(WTERRAIN(Wprec,n1,n2) == m->no_data)
        continue;
      //printf("P%i | i : %i j: %i v: %i valWprec: %f nodata:%f\n",rank,i,j,v, WTERRAIN(Wprec,n1,n2), m->no_data);
      CHECK(TERRAIN(m,i,j)>m->no_data);
      CHECK(WTERRAIN(Wprec,i,j)>m->no_data);
      CHECK(WTERRAIN(Wprec,n1,n2)>m->no_data);

      // il est important de mettre cette valeur dans un temporaire, sinon le
      // compilo fait des arrondis flotants divergents dans les tests ci-dessous
      const float Wn = WTERRAIN(Wprec,n1,n2) + EPSILON;
      if(TERRAIN(m,i,j) >= Wn)
      {
        //printf("p%i terrain i:%i j:%i=%f\n", rank, i, j, TERRAIN(m,i,j));
        WTERRAIN(W,i,j) = TERRAIN(m,i,j);
        modif = 1;
        #ifdef DARBOUX_PPRINT
        if(WTERRAIN(W,i,j)<min_darboux)
          min_darboux = WTERRAIN(W,i,j);
        #endif
      }
      else if(WTERRAIN(Wprec,i,j) > Wn)
      {
        //printf("p%i Wterrain i:%i j:%i=%f\n", rank, i, j, TERRAIN(m,i,j));
        WTERRAIN(W,i,j) = Wn;
        modif = 1;
        #ifdef DARBOUX_PPRINT
        if(WTERRAIN(W,i,j)<min_darboux)
          min_darboux = WTERRAIN(W,i,j);
        #endif
      }
    }
  }
  return(modif);
}




void printW(float *w,int nrows, int ncols, int rank){
  printf("*********** %i ***********\n", rank);
  for(int i = 0 ; i < nrows; i++){
      for(int j = 0 ; j < ncols ; j++){
        printf("||%f||",WTERRAIN(w,i,j));
      }
      printf("\n");  
    }
      printf("*********** /%i ************\n",rank);  
      printf("\n\n");
}


/*****************************************************************************/
/*           Fonction de calcul principale - À PARALLÉLISER                  */
/*****************************************************************************/
// applique l'algorithme de Darboux sur le MNT m, pour calculer un nouveau MNT
mnt *darboux(const mnt *restrict m)
{ 
  float* restrict Wres;
  int rank ; 
  int nbproc; 
  int localModif ;
  float maxGlobal, maxLocal; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

  MPI_Request request;

  //nrows correspond à la longueur de W, avec les 2 lignes supplémentaires inclus
  const int ncols = m->ncols ; 
  const int nrows = m->nrows ;
  // initialisation
  float *restrict W, *restrict Wprec;
  //On initialise W et Wprec avec 2 lignes supplémentaires, qui sont celles des processus suivant et précédant
  CHECK((W = malloc(ncols  * nrows * sizeof(float))) != NULL);  

  //Calcul du max local puis reduction pour obtenir le max global
  maxLocal = max_terrain_kernel(m,rank,nbproc)+10.;
  MPI_Allreduce(&maxLocal, &maxGlobal, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

  Wprec = init_W_kernel(m,maxGlobal);
  
  
  //Pointeur temporaire, utilisé pour l'envoi de ligne
  float *restrict Wp; 
  
  // calcul : boucle principale
  int modif = 1;  

  while(modif){
    
    modif = 0;
    
    Wp = Wprec;
    localModif = 0; // sera mis à 1 s'il y a une modification
    //1er : envoi de sa première ligne au précédant (sauf p = 1)
    //2eme : réception de sa ligne suivante, non calculée (sauf p = n - 1)
    //3eme : réception de sa ligne précédante, calculée (sauf p = 1)
    //4eme : calcul jusqu'a dernière ligne 
    //5eme : envoi de sa dernière ligne au suivant (sauf p = n - 1)
    
    //Envoi première ligne au précédent
    if (rank) {
      Wp = Wprec + ncols;    
      MPI_Isend(Wp, ncols, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD,&request);
      Wp = Wprec ;
    }

    //reception ligne suivante
    if (rank != nbproc - 1){
      Wp = Wprec + ncols*(nrows-1) ;
      MPI_Recv(Wp,ncols, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, NULL);
    }

    //reception ligne prec calculee  
    if(rank){
      MPI_Recv(Wprec, ncols, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, NULL);
    }


    if(!rank){
      for(int i=0 ; i < nrows-1;  i++)
      {
        for(int j=0; j < ncols; j++)
        {
          // calcule la nouvelle valeur de W[i,j]
          // en utilisant les 8 voisins de la position [i,j] du tableau Wprec
          localModif |= calcul_Wij(W, Wprec, m, i, j);
        }
      }  
    }else if(rank == nbproc-1){
      for(int i=1 ; i < nrows;  i++)
      {
        for(int j=0; j < ncols; j++)
        {
          localModif |= calcul_Wij(W, Wprec, m, i, j);
        }
      
      } 
      }else{
      for(int i=1 ; i < nrows-1;  i++)
      {
        for(int j=0; j < ncols; j++)
        {
          localModif |= calcul_Wij(W, Wprec, m, i, j);
        }
      } 
    }
    
    //printf("modif proc %i : %i \n",rank, localModif);
    //Envoi de la dernière ligne, calculée
    if (rank != nbproc - 1){
      Wp = W + ncols * (nrows - 2);
      MPI_Isend(Wp, ncols, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD,&request);
    }
      
    
    // #ifdef DARBOUX_PPRINTT
    // //printf("Je suis le proc %i est je print : \n",rank);
    // dpprint();
    // #endif

    // échange W et Wprec
    // sans faire de copie mémoire : échange les pointeurs sur les deux tableaux
    float *tmp = W;
    W = Wprec;
    Wprec = tmp;
  
    //Réduction de la variable modif et envoie aux autres processus
    MPI_Allreduce(&localModif, &modif, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  }
   // fin du while principal 

  //On envoie les valeurs des lignes de W au processus 0  
  
  if(!rank){
    int* displs = (int *) malloc(nbproc * sizeof(int));
    int* counts_recv =(int *) malloc(nbproc  * sizeof(int));

    // fin du calcul, le résultat se trouve dans W de chaque processus
    int taille_chunk = m->nrowsTemp / nbproc ;
    int reste = m->nrowsTemp  % nbproc ;

    displs[0] = 0 ;
    counts_recv[0] = taille_chunk * (ncols + (reste > 0 ? 1 : 0));
    //printf("nrows %i, taille chunk %i \n reste %i\n count %i disp is %i\n",m->nrowsTemp ,taille_chunk, reste ,counts_recv[0],displs[0]) ; 
    for(int i = 1 ; i < nbproc  ; i++){
      int stride = (reste > i ? 1 : 0) ;
      displs[i] = displs[i-1] + taille_chunk * (ncols + stride) ;
      counts_recv[i] = taille_chunk * (ncols + stride) ;
      printf(" i is %i , count %i disp is %i\n",i,counts_recv[i],displs[i]) ; 
    }
    
    CHECK((Wres = malloc(ncols * m->nrowsTemp * sizeof(float))) != NULL);  
    MPI_Gatherv(W , (m->nrows-1) * m->ncols , MPI_FLOAT , Wres , counts_recv , displs , MPI_FLOAT , 0 , MPI_COMM_WORLD);
  } else {
    //Gatherv des senders
    //Décalage  w
    W += ncols;
  
    int ndata = (rank == nbproc - 1 ? 1 : 2) ; 
    MPI_Gatherv(W , (m->nrows-ndata)* m->ncols , MPI_FLOAT , NULL,  NULL, NULL , MPI_FLOAT , 0 , MPI_COMM_WORLD);
  }

  free(Wprec);  
  mnt *res = NULL;
  // crée la structure résultat et la renvoie
  if (!rank){
    
    CHECK((res=malloc(sizeof(*res))) != NULL);
    memcpy(res, m, sizeof(*res));
    res->terrain = Wres;
  }
  return(res);
  
}