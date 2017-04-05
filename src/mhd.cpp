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

// wexact {{{
void Wexact(real* x, real* y, real* W){

#ifdef _1D
    real YL[_M], YR[_M], WR[_M], WL[_M];

//    // Test de Choc fort
//    YL[0] = 3;
//    YL[1] = 1.3;
//    YL[3] = 0;
//    YL[4] = 0;
//    YL[2] = 3;
//    YL[5] = 1;
//    YL[6] = 1;
//    YL[7] = 1.5;
//    YL[8] = 0;
//
//    YR[0] = 1;
//    YR[1] = 1.3;
//    YR[3] = 0;
//    YR[4] = 0;
//    YR[2] = 1;
//    YR[5] = 0.0707372016677029;
//    YR[6] = 0.9974949866040544;
//    YR[7] = 1.5;
//    YR[8] = 0;

//    // Test de Brio et Wu
//    YL[0] = 1;
//    YL[1] = 0;
//    YL[3] = 0;
//    YL[4] = 0;
//    YL[2] = 1;
//    YL[5] = 1;
//    YL[6] = 0;
//    YL[7] = 0.75;
//    YL[8] = 0;
//
//    YR[0] = 0.125;
//    YR[1] = 0;
//    YR[3] = 0;
//    YR[4] = 0;
//    YR[2] = 0.1;
//    YR[5] = -1;
//    YR[6] = 0;
//    YR[7] = 0.75;
//    YR[8] = 0;

    //Test de Dai et Woodward
    YL[0] = 1.08;
    YL[1] = 1.2;
    YL[3] = 0.01;
    YL[4] = 0.5;
    YL[2] = 0.95;
    YL[5] = 1.0155412503859613165;
    YL[6] = 0.56418958354775628695;
    YL[7] = 1.1283791670955125739;
    YL[8] = 0;

    YR[0] = 1;
    YR[1] = 0;
    YR[3] = 0;
    YR[4] = 0;
    YR[2] = 1;
    YR[5] = 1.1283791670955125739;
    YR[6] = 0.56418958354775628695;
    YR[7] = 1.1283791670955125739;
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

    Y[0] = gam*gam;
    Y[1] = -sin(*y);
    Y[2] = gam;
    Y[3] = sin(*x);
    Y[4] = 0.0;
    Y[5] = sin(2*(*x));
    Y[6] = 0.0;
    Y[7] = -sin(*y);
    Y[8] = 0.0;

    conservatives(Y, W);
#endif

}
// }}}






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
                fic << "\"P  \""<<endl; // nom de la vue
                break;
            case 3:
                fic << "\"UY  \""<<endl; // nom de la vue
                break;
            case 4:
                fic << "\"UZ\""<<endl; // nom de la vue
                break;
            case 5:
                fic << "\"BY  \""<<endl; // nom de la vue
                break;
            case 6:
                fic << "\"BZ\""<<endl; // nom de la vue
                break;
            case 7:
                fic << "\"BX  \""<<endl; // nom de la vue
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

                k=(i+j*_NXTRANSBLOCK)+3*_NXTRANSBLOCK*_NYTRANSBLOCK;
                uy=(Wn1[k])/rho;

                k=(i+j*_NXTRANSBLOCK)+4*_NXTRANSBLOCK*_NYTRANSBLOCK;
                uz=(Wn1[k])/rho;

                k=(i+j*_NXTRANSBLOCK)+5*_NXTRANSBLOCK*_NYTRANSBLOCK;
                by=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+6*_NXTRANSBLOCK*_NYTRANSBLOCK;
                bz=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+7*_NXTRANSBLOCK*_NYTRANSBLOCK;
                bx=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+8*_NXTRANSBLOCK*_NYTRANSBLOCK;
                psi=(Wn1[k]);

                k=(i+j*_NXTRANSBLOCK)+2*_NXTRANSBLOCK*_NYTRANSBLOCK;
                p=(gam-1)*(Wn1[k] - 0.5*rho*(ux*ux+uy*uy+uz*uz) - 0.5*(bx*bx+by*by+bz*bz));


                switch(typplot){
                    case 0:
                        val=rho;
                        break;
                    case 1:
                        val=ux;
                        break;
                    case 2:
                        val=p;
                        break;
                    case 3:
                        val=uy;
                        break;
                    case 4:
                        val=uz;
                        break;
                    case 5:
                        val=by;
                        break;
                    case 6:
                        val=bz;
                        break;
                    case 7:
                        val=bx;
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
    real dtt = 0;
    for(real t=0;t<_TMAX; t=t+dtt){

        cout << "Iter="<<iter++<< endl;;
        TimeStepCPU(Wn1,&dtt);
        cout << t << endl;
    }

#ifdef _1D
    GnuPlot(Wn1);
#endif
    PlotGmshBinary(Wn1);
    return 0;
}
