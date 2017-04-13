// vim: set sw=4 ts=4 sts=4 et tw=78 foldmarker={{{,}}} foldlevel=0 foldmethod=marker spell:

// cas {{{
//#define _1D
#define _2D
// }}}

// include {{{
#include <math.h>
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;
// }}}

// define {{{

// Dans ce cas sequentiel, pour modifier le nombre de mailles il ne faut
// modifier que _NBWORKSX et _NBWORKSY

#define _M (9)                                      //Number of conservative variable
#define _TMAX (1)
#define _GAP (0)                                    // taille du recouvrement entre les work groupe
#define _TRANSBLOCK 1                               // size of the cached transposed block
#define _NBWORKSX (1<<7)                            // number of work-items in a work-group
#define _NBLOCKSX (1)                               // number of work-groups
#define _NX (_NBLOCKSX*(_NBWORKSX-2*_GAP)+2*_GAP)   // number of volume finite
#define _NXTRANSBLOCK ( (_NX%_TRANSBLOCK==0)? _NX : _NX+_TRANSBLOCK-_NX%_TRANSBLOCK )

#define _NBWORKSY (1<<7)
#define _NBLOCKSY (1)
#define _NY (_NBLOCKSY*(_NBWORKSY-2*_GAP)+2*_GAP)
#define _NYTRANSBLOCK ( (_NY%_TRANSBLOCK==0)? _NY : _NY+_TRANSBLOCK-_NY%_TRANSBLOCK )
#define _CFL (0.5)
#define _SPLIT (1)                                  // affiche 1 maille sur _SPLIT dans le fichier .msh

#define Min(a,b) (((a) < (b)) ? (a) : (b))
#define Max(a,b) (((a) > (b)) ? (a) : (b))
#define Abs(a) ((a) > (0) ? (a) : (-a))

#ifdef _1D
    #define _LONGUEURX (10)                         //Longueur du domaine suivant x
    #define _LONGUEURY (10)                         //Longueur du domaine suivant y
    #define _XMIN (-5)
    #define _XMAX (5)
    #define _YMIN (-5)
    #define _YMAX (5)
#endif
#ifdef _2D
    // Orzag Tang
    #define _LONGUEURX (6.2831853)                  //Longueur du domaine suivant x
    #define _LONGUEURY (6.2831853)                  //Longueur du domaine suivant y
    #define _XMIN (0)
    #define _XMAX (6.2831853)
    #define _YMIN (0)
    #define _YMAX (6.2831853)
#endif



//#define _GAM (2)
#define _GAM (1.666666666666)
#define _PI (3.14159265359)
#define _CH (5)

typedef float real;

// }}}

void conservatives(real Y[_M], real W[_M])
{
  int i;
  W[0] = Y[0];
  for(i = 1; i < 4; i++)
  {
    W[i] = Y[i]*Y[0];
  }
  W[4] = (_GAM-1)*(Y[4] + Y[0]*(Y[1]*Y[1] + Y[2]*Y[2] + Y[3]*Y[3])/2 + (Y[5]*Y[5] + Y[6]*Y[6] + Y[7]*Y[7])/2);
  for(i = 5; i < _M; i++)
  {
    W[i] = Y[i];
  }
}

void primitives(real Y[_M], real W[_M])
{
  int i;
  Y[0] = W[0];
  for(i = 1; i < 4; i++)
  {
    Y[i] = W[i]/W[0];
  }
  Y[4] = W[4]/(_GAM-1) + Y[0]*(Y[1]*Y[1] + Y[2]*Y[2] + Y[3]*Y[3])/2 + (W[5]*W[5] + W[6]*W[6] + W[7]*W[7])/2;
  for(i = 5; i < _M; i++)
  {
    Y[i] = W[i];
  }
}

real* flux(real W[_M], real n[3])
{
  real* F = (real*)malloc(_M*sizeof(real));
  real Y[_M];
  primitives(Y, W);
  int i;
  real BB = W[5]*W[5] + W[6]*W[6] + W[7]*W[7];
  real Bn = W[5]*n[0] + W[6]*n[1] + W[7]*n[2];
  real Bu = W[5]*Y[1] + W[6]*Y[2] + W[7]*Y[3];
  real un = Y[1]*n[0] + Y[2]*n[1] + Y[3]*n[2];

  F[0] = W[0]*un;
  for(i = 1; i<4; i++)
  {
    F[i] = un*W[i] + (Y[2] + BB/2)*n[i-1] - Bn*W[i+4];
  }
  F[4] = (W[4] + Y[4] + BB/2)*un - Bu*Bn;
  for(i = 5; i < 8; i++)
  {
    F[i] = un*W[i] - Bn*Y[i-4] + W[8] * n[i-5];
  }
  F[8] = _CH*_CH*Bn;
  return F;
}

// wexact {{{
void Wexact(real* x, real* y, real* W){

#ifdef _1D
    real YL[_M], YR[_M], WR[_M], WL[_M];

//    // Test de Choc fort
//    YL[0] = 3;
//    YL[1] = 1.3;
//    YL[4] = 0;
//    YL[2] = 0;
//    YL[3] = 3;
//    YL[6] = 1;
//    YL[7] = 1;
//    YL[5] = 1.5;
//    YL[8] = 0;
//
//    YR[0] = 1;
//    YR[1] = 1.3;
//    YR[4] = 0;
//    YR[2] = 0;
//    YR[3] = 1;
//    YR[6] = 0.0707372016677029;
//    YR[7] = 0.9974949866040544;
//    YR[5] = 1.5;
//    YR[8] = 0;

//    // Test de Brio et Wu
//    YL[0] = 1;
//    YL[1] = 0;
//    YL[4] = 0;
//    YL[2] = 0;
//    YL[3] = 1;
//    YL[6] = 1;
//    YL[7] = 0;
//    YL[5] = 0.75;
//    YL[8] = 0;
//
//    YR[0] = 0.125;
//    YR[1] = 0;
//    YR[4] = 0;
//    YR[2] = 0;
//    YR[3] = 0.1;
//    YR[6] = -1;
//    YR[7] = 0;
//    YR[5] = 0.75;
//    YR[8] = 0;

    //Test de Dai et Woodward
    YL[0] = 1.08;
    YL[1] = 1.2;
    YL[4] = 0.01;
    YL[2] = 0.5;
    YL[3] = 0.95;
    YL[6] = 1.0155412503859613165;
    YL[7] = 0.56418958354775628695;
    YL[5] = 1.1283791670955125739;
    YL[8] = 0;

    YR[0] = 1;
    YR[1] = 0;
    YR[4] = 0;
    YR[2] = 0;
    YR[3] = 1;
    YR[6] = 1.1283791670955125739;
    YR[7] = 0.56418958354775628695;
    YR[5] = 1.1283791670955125739;
    YR[8] = 0;


    conservatives(YL, WL);
    conservatives(YR, WR);


    if(*x < 0)
        for(int i=0; i<_M; i++){
            W[i] = WL[i];
        }
    else
        for(int i=0; i<_M; i++){
            W[i] = WR[i];
        }
#endif
#ifdef _2D
// Orzag-Tang
    real Y[_M];

    real gam = _GAM;

    Y[0] = gam*gam; //rho
    Y[1] = -sin(*y); //Ux
    Y[4] = gam; //p
    Y[2] = sin(*x); //Uy
    Y[3] = 0.0; //Uz
    Y[6] = sin(2*(*x)); //By
    Y[7] = 0.0; //Bz
    Y[5] = -sin(*y);  //Bx
    Y[8] = 0.0; //psi

    conservatives(Y, W);
#endif

}
// }}}

real* matFoisVect(real A[_M][_M], real V[_M])
{
  real* res = (real*)malloc(_M*sizeof(real));
  int i,j;
  for(i=0; i < _M; i++)
  {
    res[i]=0;
    for(j=0; j < _M; j++)
    {
      res[i] += A[i][j]*V[j];
      if(isnan(res[i]))
      {
        cout << A[i][j] << " fois " << V[j] << " = " << res[i] << endl;
        exit(1);
      }
    }
  }
  return res;
}

real* addVect(real V1[_M], real V2[_M], int soustraction)
{
  real* res = (real*)malloc(_M*sizeof(real));
  int i;
  for(i=0; i<_M; i++)
  {
    if(soustraction)
      res[i] = V1[i] - V2[i];
    else
      res[i] = V1[i] + V2[i];
  }
  return res;
}

void realFoisVect(real V[_M], real a)
{
  int i;
  for(i=0; i<_M; i++)
  {
    V[i] = V[i]*a;
  }
}

void matA(real A[_M][_M], real W[_M], real n[3])
{
  real Y[_M];
  int i,j;
  primitives(Y, W);
  real Bn = W[5]*n[0] + W[6]*n[1] + W[7]*n[2];
  real Bu = W[5]*Y[1] + W[6]*Y[2] + W[7]*Y[3];
  real un = Y[1]*n[0] + Y[2]*n[1] + Y[3]*n[2];
  for(i=0; i < _M; i++)
  {
    for(j=0; j<_M; j++)
    {
      A[i][j] = 0;
    }
  }

  for(i=0; i < 5; i++)
  {
    A[i][i] = un;
  }
  for(i=0; i < 3; i++)
  {
    A[0][i+1] = Y[0] * n[i];
    A[i+1][4] = n[i] / Y[0];
    A[4][i+1] = _GAM * Y[4]*n[i];
    A[i+1][i+5] = - Bn / Y[0];
    A[4][i+5] = (_GAM - 1)*Bu*n[i];
    A[8][i+5] = _CH*_CH*n[i];
    A[i+5][8] = n[i];
  }
  A[4][8] = (1 - _GAM)*Bn;
  A[5][2] = Y[5]*n[1];
  A[5][3] = Y[5]*n[2];
  A[6][1] = Y[6]*n[0];
  A[6][3] = Y[6]*n[2];
  A[7][1] = Y[7]*n[0];
  A[7][2] = Y[7]*n[1];
  A[5][5] = Y[2]*n[1] + Y[3]*n[2];
  A[6][6] = Y[1]*n[0] + Y[3]*n[2];
  A[7][7] = Y[1]*n[0] + Y[2]*n[1];
  A[5][6] = -Y[1]*n[1];
  A[5][7] = -Y[1]*n[2];
  A[6][5] = -Y[2]*n[0];
  A[6][7] = -Y[2]*n[2];
  A[7][5] = -Y[3]*n[0];
  A[7][6] = -Y[3]*n[1];
  A[5][1] = -Y[6]*n[1] - Y[7]*n[2];
  A[6][2] = -Y[5]*n[0] - Y[7]*n[2];
  A[7][3] = -Y[5]*n[0] - Y[6]*n[1];
  //moins sûr
  A[1][6] = (/*Y[6]*n[0]*/ - Y[5]*n[1])/Y[0];
  A[1][7] = (/*Y[7]*n[0]*/ - Y[5]*n[2])/Y[0];
  A[2][5] = (-Y[6]*n[0]/* + Y[5]*n[1]*/)/Y[0];
  A[2][7] = (-Y[6]*n[2]/* + Y[7]*n[1]*/)/Y[0];
  A[3][5] = (-Y[7]*n[0]/* + Y[5]*n[2]*/)/Y[0];
  A[2][7] = (/*Y[6]*n[2] */- Y[7]*n[1])/Y[0];
}

void Ref2PhysMap(real* xx, real* yy, real* x, real* y)
{
  *x = *xx * _LONGUEURX + _XMIN;
  *y = *yy * _LONGUEURY + _YMIN;
}

void InitData(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M])
{
  int i,j,k;
  real x,y;
  real xx, yy;
  real W[_M];
  for(j = 0; j < _NYTRANSBLOCK; j++)
  {
    for(i = 0; i < _NXTRANSBLOCK; i++)
    {
      xx = ((real)i)/_NXTRANSBLOCK;
      yy = ((real)j)/_NYTRANSBLOCK;
      Ref2PhysMap(&xx, &yy, &x, &y);
      Wexact(&x, &y, W);
      for(k = 0; k < _M; k++)
      {
        Wn1[j*_NXTRANSBLOCK + i + k*(_NXTRANSBLOCK*_NYTRANSBLOCK)] = W[k];
      }
    }
  }
}

real* fluxPolynomial(real WL[_M], real WR[_M], real n[3])
{
    int i;
    real* res = (real*)malloc(_M*sizeof(real));
    real* flux_WL = flux(WL,n);
    real* flux_WR = flux(WR,n);

    real* V = NULL;
    real* dW = addVect(WR, WL, 1);
    real A[_M][_M];
    real* Wmil = addVect(WL, WR, 0);
    realFoisVect(Wmil, 0.5);
    matA(A, Wmil, n);

    real* AdW = matFoisVect(A, dW);
    real* A2 = matFoisVect(A, AdW);
    realFoisVect(dW, -5./16);
    real* tmp = addVect(dW, A2, 0);
    realFoisVect(dW, -3);
    real* A3 = matFoisVect(A, tmp);
    real* A4 = matFoisVect(A, A3);
    free(tmp);
    tmp = addVect(dW, A4, 0);
    real* A5 = matFoisVect(A, tmp);
    real* A6 = matFoisVect(A, A5);
    realFoisVect(dW, 1./3);

    V = addVect(dW, A6, 0);

    for (i=0;i<_M;i++)
    {
      res[i] = 0.5*(flux_WL[i]+flux_WR[i]) - 0.5*V[i];
    }
    free(flux_WL);
    free(flux_WR);
    free(Wmil);
    free(dW);
    free(AdW);
    free(A2);
    free(A3);
    free(A4);
    free(A5);
    free(A6);
    free(V);
    free(tmp);
    return res;
}

real* fluxRusanov(real WL[_M],real WR[_M],real n[3]){
    real* res = (real*)malloc(_M*sizeof(real));
    real* flux_WL = flux(WL,n);
    real* flux_WR = flux(WR,n);
    int i;

    for (i=0;i<_M;i++)
      res[i] = 0.5*(flux_WL[i]+flux_WR[i]) - _CH*0.5*(WR[i]-WL[i]);
    free(flux_WL);
    free(flux_WR);
    return res;
}

void vectorWn1(real res[_M],real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M],int i, int j){
  int k;

  for (k=0;k<_M;k++)
    res[k]=Wn1[i+j*_NXTRANSBLOCK+k*_NXTRANSBLOCK*_NYTRANSBLOCK];
}


void TimeStepCPU1D(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M], real* dtt)
{
  real dxx=_LONGUEURX*1.0/_NXTRANSBLOCK;
  real Wn1bis[_NXTRANSBLOCK*_NYTRANSBLOCK*_M];
  real W_middle[_M],W_left[_M],W_right[_M];
  real* fluxMR = NULL;
  real* fluxLM = NULL;
  real nx[3]={1,0,0};
  int i,j,k;
  int caseijk;
  real zero = 0;
  real max = _NXTRANSBLOCK;

  /*Copie de Wn1 dans Wn1bis*/
  for (i=0;i<_NXTRANSBLOCK*_NYTRANSBLOCK*_M;i++)
    Wn1bis[i]=Wn1[i];

  for (j=0;j<_NYTRANSBLOCK;j++){
    for (i=0;i<_NXTRANSBLOCK;i++){
      vectorWn1(W_middle,Wn1bis,i,j);

      if (i==0) /*effet de bord à gauche*/
        Wexact(&zero,&zero,W_left);
      else
        vectorWn1(W_left,Wn1bis,i-1,j);

      if (i==_NXTRANSBLOCK-1) /*effet de bord à droite*/
        Wexact(&max,&zero,W_right);
      else
        vectorWn1(W_right,Wn1bis,i+1,j);

      fluxMR=fluxRusanov(W_middle,W_right,nx);
      fluxLM=fluxRusanov(W_left,W_middle,nx);
      for (k=0;k<_M;k++){
        caseijk=i+j*_NXTRANSBLOCK+k*_NXTRANSBLOCK*_NYTRANSBLOCK;
        Wn1[caseijk]=W_middle[k]-(*dtt/dxx)*(fluxMR[k]-fluxLM[k]);
      }
      free(fluxMR);
      free(fluxLM);
    }
  }
}

void TimeStepCPU2D(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M], real* dtt)
{
  real dxx=_LONGUEURX*1.0/_NXTRANSBLOCK;
  real dyy=_LONGUEURY*1.0/_NYTRANSBLOCK;
  real Wn1bis[_NXTRANSBLOCK*_NYTRANSBLOCK*_M];
  real W_middle[_M],W_left[_M],W_right[_M],W_up[_M],W_down[_M];
  real* fluxMR = NULL;
  real* fluxLM = NULL;
  real* fluxMU = NULL;
  real* fluxDM = NULL;
  real nx[3]={1,0,0};
  real ny[3]={0,1,0};
  int i,j,k;
  int caseijk;

  /*Copie de Wn1 dans Wn1bis*/
  for (i=0;i<_NXTRANSBLOCK*_NYTRANSBLOCK*_M;i++)
    Wn1bis[i]=Wn1[i];

  for (j=0;j<_NYTRANSBLOCK;j++){
    for (i=0;i<_NXTRANSBLOCK;i++){
      vectorWn1(W_middle,Wn1bis,i,j);

      if (i==0) /*effet de bord à gauche*/
        vectorWn1(W_left,Wn1bis,_NXTRANSBLOCK-1,j);
      else
        vectorWn1(W_left,Wn1bis,i-1,j);

      if (i==_NXTRANSBLOCK-1) /*effet de bord à droite*/
        vectorWn1(W_right,Wn1bis,0,j);
      else
        vectorWn1(W_right,Wn1bis,i+1,j);

      if (j == 0)  /*effet de bord en bas*/
        vectorWn1(W_down,Wn1bis,i,_NYTRANSBLOCK-1);
      else
        vectorWn1(W_down,Wn1bis,i,j-1);

      if (j == _NYTRANSBLOCK-1)  /*effet de bord en haut*/
        vectorWn1(W_up,Wn1bis,i,0);
      else
        vectorWn1(W_up,Wn1bis,i,j+1);

      fluxMR=fluxRusanov(W_middle,W_right,nx);
      fluxLM=fluxRusanov(W_left,W_middle,nx);
      fluxMU=fluxRusanov(W_middle,W_up,ny);
      fluxDM=fluxRusanov(W_down,W_middle,ny);

      /*fluxMR=fluxPolynomial(W_middle,W_right,nx);
      /*fluxLM=fluxPolynomial(W_left,W_middle,nx);
      fluxMU=fluxPolynomial(W_middle,W_up,ny);
      fluxDM=fluxPolynomial(W_down,W_middle,ny);*/

      for (k=0;k<_M;k++){
        caseijk=i + j*_NXTRANSBLOCK + k*_NXTRANSBLOCK*_NYTRANSBLOCK;
        Wn1[caseijk]=W_middle[k]-(*dtt/dxx)*(fluxMR[k]-fluxLM[k])-(*dtt/dyy)*(fluxMU[k]-fluxDM[k]);
      }
      free(fluxMR);
      free(fluxLM);
      free(fluxMU);
      free(fluxDM);
    }
  }
}

// gnuplot {{{
void GnuPlot(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M])
{

    system("if [ ! -d RESU ] ; then mkdir RESU; fi");

    // fichiers pour les tracés gnuplot
    ofstream ficr("RESU/r");
    ofstream ficux("RESU/ux");
    ofstream ficp("RESU/p");
    ofstream ficuy("RESU/uy");
    ofstream ficuz("RESU/uz");
    ofstream ficby("RESU/by");
    ofstream ficbz("RESU/bz");

    //  Wn1.GPU2CPU();

    // ici, ce sont les variables conservatives et non primitives pour rho ce n'est pas dérangeant
    // mais pour les autre, il faut y penser
    //cout << Wn1;

    for(int j=0;j<_NY;j++){
        for(int i = 0; i < _NX; i++) {
            real xx=1./_NX * (i+0.5);
            real yy=0;
            real x,y;

            real gam = _GAM;

            real r  = Wn1[i+0*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];
            real ux = Wn1[i+1*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];
            real p  = Wn1[i+2*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];
            real uy = Wn1[i+3*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];
            real uz = Wn1[i+4*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];
            real by = Wn1[i+5*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];
            real bz = Wn1[i+6*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];
            real bx = Wn1[i+7*_NXTRANSBLOCK*_NYTRANSBLOCK+j*_NXTRANSBLOCK];

            Ref2PhysMap(&xx,&yy,&x,&y);
            ficr << x<< " " << r << endl;
            ficux << x<< " " << ux/r << endl;
            ficp << x<< " " << (gam-1)*(p-0.5*r*((ux/r)*(ux/r)+(uy/r)*(uy/r)+(uz/r)*(uz/r))-0.5*(bx*bx+by*by+bz*bz)) << endl;
            ficuy << x<< " " << uy/r << endl;
            ficuz << x<< " " << uz/r << endl;
            ficby << x<< " " << by << endl;
            ficbz << x<< " " << bz << endl;
        }
        ficr << endl;
        ficux << endl;
        ficp << endl;
        ficuy << endl;
        ficuz << endl;
        ficby << endl;
        ficbz << endl;
    }
    ficr.close();
    ficux.close();
    ficp.close();
    ficuy.close();
    ficuz.close();
    ficby.close();
    ficbz.close();


    // creation du script de lancement de gnuplot
    ofstream gdat("gdat");

    gdat << "set terminal x11" << endl;
    gdat << "plot 'RESU/r' w l" << endl;
    gdat << "pause -1" << endl;
    gdat << "press return";

    gdat.close();

    system("gnuplot gdat");

}
// }}}

// Gmsh {{{
void PlotGmshBinary(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M]){


    ofstream fic("clmhd.msh",ios::binary);
    real x,y,xx,yy;

    int split=_SPLIT;//Pour le tracé gmsh, on ne trace qu'une maille sur _SPLIT car sinon on n'arrive pas à tracer quand le maillage est trop gros (1000*5000)


    fic << "$MeshFormat"<<endl;
    fic << 2.2 <<" ";
    fic << 1 <<" ";
    fic << sizeof(double) << endl;
    int one=1;
    fic.write((char*) &one,sizeof(int));
    fic << endl << "$EndMeshFormat"<<endl;
    fic << "$Nodes" << endl;

    // number of nodes
    fic << _NX/split * _NY/split * 4 <<endl;

    const int di[4]={0,1,1,0};
    const int dj[4]={0,0,1,1};


    for(int i=0;i<_NX/split;i++){
        for(int j=0;j<_NY/split;j++){
            for(int ii=0;ii<4;ii++){
                xx=1./_NX*(i+di[ii])*split;
                yy=1./_NY*(j+dj[ii])*split;
                Ref2PhysMap(&xx,&yy,&x,&y);
                int nnoe=4*(j*_NX/split+i)+ii+1;
                double xd=x,yd=y,zd=0;
                fic.write((char*) &nnoe,sizeof(int));
                fic.write((char*) &xd,sizeof(double));
                fic.write((char*) &yd,sizeof(double));
                fic.write((char*) &zd,sizeof(double));
            }
        }
    }


    fic << endl<<"$EndNodes"<<endl;
    fic << "$Elements"<<endl;
    fic << _NX/split * _NY/split <<endl;

    int elm_type=3;
    int num_elm_follow=_NX/split * _NY/split;
    int num_tags=0;

    fic.write((char*) &elm_type,sizeof(int));
    fic.write((char*) &num_elm_follow,sizeof(int));
    fic.write((char*) &num_tags,sizeof(int));

    for(int i=0;i<_NX/split;i++){
        for(int j=0;j<_NY/split;j++){

            int numelem=j*_NX/split+i+1;
            fic.write((char*) &numelem,sizeof(int));
            //fic << (j*_NX+i+1) <<" "<< 3 <<" 0 ";
            for(int ii=0;ii<4;ii++){
                int numnoe=4*(j*_NX/split+i) + ii +1;
                fic.write((char*) &numnoe,sizeof(int));
                //fic << 4*(j*_NX+i) + ii +1 <<" ";
            }
        }
    }

    fic << endl<<"$EndElements"<<endl;

    // data plots



    for(int typplot=0;typplot<_M;typplot++){
        //for(int typplot=1;typplot<2;typplot++){
        fic << "$NodeData"<<endl;
        fic << 1 <<endl;  // une valeur
        switch(typplot){
            case 0:
                fic << "\"RHO\""<<endl; // nom de la vue
                break;
            case 1:
                fic << "\"UX  \""<<endl; // nom de la vue
                break;
            case 2:
                fic << "\"UY  \""<<endl; // nom de la vue
                break;
            case 3:
                fic << "\"UZ  \""<<endl; // nom de la vue
                break;
            case 4:
                fic << "\"P\""<<endl; // nom de la vue
                break;
            case 5:
                fic << "\"BX  \""<<endl; // nom de la vue
                break;
            case 6:
                fic << "\"BY\""<<endl; // nom de la vue
                break;
            case 7:
                fic << "\"BZ  \""<<endl; // nom de la vue
                break;
            case 8:
                fic << "\"PSI \""<<endl; // nom de la vue
        }

        fic << 1 <<endl;  //une valeur
        fic << 0 << endl; // temps
        fic << 3 << endl; // 3 données
        fic << 0  << endl;  // numéro pas de temps
        fic << 1 << endl; //visu d'un scalaire

        fic << 4 * _NX/split * _NY/split <<endl;


        for(int j=0;j<_NY;j=j+split) {
            for(int i=0;i<_NX;i=i+split) {

                real rho,ux,uy,uz,bx,by,bz,p,gam,val,psi;
                //real phi, pinf;
                int k;

                gam = _GAM;

                k=(i+j*_NXTRANSBLOCK);
                rho=Wn1[k];

                k=(i+j*_NXTRANSBLOCK)+1*_NXTRANSBLOCK*_NYTRANSBLOCK;
                ux=(Wn1[k])/rho;

                k=(i+j*_NXTRANSBLOCK)+2*_NXTRANSBLOCK*_NYTRANSBLOCK;
                uy=(Wn1[k])/rho;

                k=(i+j*_NXTRANSBLOCK)+3*_NXTRANSBLOCK*_NYTRANSBLOCK;
                uz=(Wn1[k])/rho;

                k=(i+j*_NXTRANSBLOCK)+6*_NXTRANSBLOCK*_NYTRANSBLOCK;
                by=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+7*_NXTRANSBLOCK*_NYTRANSBLOCK;
                bz=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+5*_NXTRANSBLOCK*_NYTRANSBLOCK;
                bx=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+8*_NXTRANSBLOCK*_NYTRANSBLOCK;
                psi=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+4*_NXTRANSBLOCK*_NYTRANSBLOCK;
                p=(gam-1)*(Wn1[k] - 0.5*rho*(ux*ux+uy*uy+uz*uz) - 0.5*(bx*bx+by*by+bz*bz));


                switch(typplot){
                    case 0:
                        val=rho;
                        break;
                    case 1:
                        val=ux;
                        break;
                    case 2:
                        val=uy;
                        break;
                    case 3:
                        val=uz;
                        break;
                    case 4:
                        val=p;
                        break;
                    case 5:
                        val=bx;
                        break;
                    case 6:
                        val=by;
                        break;
                    case 7:
                        val=bz;
                        break;
                    case 8:
                        val=psi;
                }

                // same data at the four nodes
                for(int ii=0;ii<4;ii++){
                    int nodenumber=4*(j/split*_NX/split+i/split) + ii +1;
                    double value=val;
                    fic.write((char*) &nodenumber,sizeof(int));
                    fic.write((char*) &value,sizeof(double));
                }
            }
        }


        fic << endl<<"$EndNodeData"<<endl;
    }

    fic.close();

    }
    // }}}



int main(int argc, char const* argv[]){

    real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M];
    InitData(Wn1);

    int iter = 0;
    real dtt = _CFL*(_LONGUEURX*1.0/_NXTRANSBLOCK)/_CH;
    for(real t=0;t<_TMAX; t=t+dtt){

        cout << "Iter="<<iter++<< endl;;
        #ifdef _1D
            TimeStepCPU1D(Wn1,&dtt);
        #endif
        #ifdef _2D
            TimeStepCPU2D(Wn1,&dtt);
        #endif
        cout << t << endl;
    }

#ifdef _1D
    GnuPlot(Wn1);
#endif
    PlotGmshBinary(Wn1);
    return 0;
}
