#include <stdio.h>
#include <stdlib.h>

#define taille 3

double ** initialisation(){
  double **test;
  int i,j;
  test=calloc(taille,sizeof(double *));
  for (i=0;i<taille;i++)
    test[i]=malloc(sizeof(double));

  for (i=0;i<taille;i++){
    for (j=0;j<taille;j++)
      test[i][j]=0;
  }
  return test;
}

double ** matId(){
  double ** test=initialisation();
  int i;

  for (i=0;i<taille;i++)
    test[i][i]=1;
  return test;
}

void affichage(double ** test){
  int i,j;
  for (i=0;i<taille;i++){
    for (j=0;j<taille;j++)
      printf("%f  ",test[i][j]);
    printf("\n");
  }
}

double ** creationMatrice(double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8,double a9){
  double **test=initialisation();
  test[0][0]=a1;
  test[0][1]=a2;
  test[0][2]=a3;
  test[1][0]=a4;
  test[1][1]=a5;
  test[1][2]=a6;
  test[2][0]=a7;
  test[2][1]=a8;
  test[2][2]=a9;
  return test;
}

double ** produitMatrice(double **A, double **B){
  int i,j,k;
  double ** res;
  res=initialisation();

  for (i=0;i<taille;i++){
    for (j=0;j<taille;j++){
      for (k=0;k<taille;k++)
        res[i][j]=res[i][j]+A[i][k]*B[k][j];
    }
  }
  return res;
}

double ** matriceL(double ** A, int iter){
  int i, j;
  double ** res=initialisation();

  for (i=0;i<taille;i++){
    for (j=0;j<taille;j++){
      if (i==j)
        res[i][j]=1;
      else if (i<j)
        res[i][j]=0;
      else if (j!=iter)
        res[i][j]=0;
      else
        res[i][j]=-A[i][j]/A[iter][iter];
    }
  }
  return res;
}

double ** inverseMatL(double ** A, int iter){
  int i;
  for (i=iter+1;i<taille;i++)
    A[i][iter]=-A[i][iter];
  return A;
}

double * resolutionLy_B(double ** L,double * B){
  double * res=calloc(taille,sizeof(double));
  int i,k;

  for (i=0;i<taille;i++){
    res[i]=B[i];
    for (k=0;k<i;k++)
      res[i]-=L[i][k]*res[k];
  }
  return res;
}

double * resolutionUx_Y(double ** U,double * Y){
  double * res=calloc(taille,sizeof(double));
  int i,k;

  for (i=0;i<taille;i++){
    res[taille-1-i]=Y[taille-1-i];
    for (k=0;k<i;k++)
      res[taille-1-i]-=U[taille-1-i][taille-1-k]*res[taille-1-k];
  }
  return res;
}

int main(){
  /*double ** test;
  double ** test1;
  double ** test2;*/
  double ** testL;
  double * solution;
  int i;
  testL=creationMatrice(1,1,2,1,-1,-1,1,0,1);
  /*test=creationMatrice(2,-1,0,4,1,-1,0,6,2);
  test1=creationMatrice(1,0,0,-2,1,0,0,0,1);
  test2=creationMatrice(2,-1,0,4,1,-1,0,6,2);*/
  double B[taille]={5,1,3};

  /*test[0][0]=1;
  test[0][1]=1;
  test[0][2]=2;
  test[1][0]=1;
  test[1][1]=-1;
  test[1][2]=-1;
  test[2][0]=1;
  test[2][1]=0;
  test[2][2]=1;*/

  /*Affichage de la matrice initiale*/
  printf("Matrice initiale test :\n");
  affichage(testL);

  printf("\nMéthode LU :\n");
  int iter=0;
  double ** tmp=initialisation();
  double ** matL=matId();
  double ** matLfinal=matId();
  for (i=0;i<taille-1;i++){
    tmp=matriceL(testL,iter);
                    printf("Matrice intermédiaire L à l'itération %d\n",iter+1);
                    affichage(tmp);
    testL=produitMatrice(tmp,testL);
                    printf("Matrice à l'iteration %d\n",iter+1);
                    affichage(testL);
    matLfinal=produitMatrice(matLfinal,inverseMatL(tmp,iter)); /*contiendra la matrice L*/
    iter++;
    printf("\n");
  } /*A la sortie de la boucle, testL est une matrice triangulaire supérieure
  matLfinal la matrice triangulaire inférieure*/
  printf("\nLa matrice tri supérieure à la dernière itération\n");
  affichage(testL);
  printf("\nLa matrice tri inférieure à la dernière itération\n");
  affichage(matLfinal);
  printf("\nLe produit LU à la dernière itération\n");
  affichage(produitMatrice(matLfinal,testL));

  printf("\nSolution par la méthode LU\n");
  solution=resolutionUx_Y(testL,resolutionLy_B(matLfinal,B));
  for (i=0;i<taille;i++)
    printf("%3f  ",solution[i]);
  printf("\n");
  free(tmp);
  /*free(test);
  free(test1);
  free(test2);*/
  free(matL);

  return 0;
}
