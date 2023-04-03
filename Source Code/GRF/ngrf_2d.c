//Generate diffuse image at different freq. in FITS fromat
// uses input angular power spectrum of brightness temperature fluctuation 

#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<fitsio.h>
#include<unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double k_B=1.38e3 /* in Jy*/, c=3.0e8;



//return input angular power spectrum
double  P_T(double UU)// return input brightness temperature PS in K^2
{
  double ell,pu;  

  double Av=1e4, betav=-1.; // amplitude and multipole exponent of APS

  ell=2.0*M_PI*UU;// calculate ell --. multipole


  //pu=Av*pow(ell/1000.0,betav);
  
  /////////////////////////////////////////////////////////////

  pu=1.0; // unit angular power spectrum in K^2

  /////////////////////////////////////////////////////////////

  return(pu);
}


int main(int argc, char *argv[])
{

  double *out, *img;
  int NN, length;
  double LL, pixel, ff,nu0,nu; 
  unsigned long int seed;
  char  OPT, task[128];
  int Nchan;
  double chan0, deltanu, spindex;
  

  FILE *fp;

  void Fill_Modes(int NN,double LL,unsigned long int seed,double (*func)(double),double *out);
  void Make_Image(int NN,double LL,double *out);
  void Write_Fits_Image(int NN,double LL,double *out,char *outfile);
  void Non_Gauss(int NN,double LL,double *out,double ff);
  void Get_Modes(int NN,double LL,double *out);
  void power_spec(int NN,double LL,double *out,int NB,double *usum,double *pksum,int *no);
  //void PB_Flux(int NN,double LL, double *out, double nu); 
  void  PB_Flux_MF(int NN,double LL, double *out, double *img, int Nchan, double nu0,double chan0, double spindex, double deltanu);
  double Beam(double theta,double freq);
  void Write_Fits_MF_Image(int NN,double LL, int Nchan, double nu0,double chan0, double deltanu, double *img,char *outfile);



  
  // reads name of parameter file and output FITS file at runtime
  if(argc!=12)
   {
      printf("Usage: %s <N no. of pixs.><L gridsize (arcmin)><seed int <0><non-Gassian parameter ff><spectral index (Intensity)><nu0><channel width><chan0><Number of channels><ouput FITS file SF><ouput FITS file MF>\n", argv[0]);
      return 1;
    }
    
  //system(task);
  sscanf(argv[1],"%d",&NN); // 2048
  sscanf(argv[2],"%lf",&LL); //0.515= 1.453e-3 radian
  sscanf(argv[3],"%ld",&seed); //7999
  sscanf(argv[4],"%lf",&ff); //0.5
  sscanf(argv[5],"%lf",&spindex); //0.8
  sscanf(argv[6],"%lf",&nu0); //150.00e06 Hz
  sscanf(argv[7],"%lf",&deltanu); //62.50e3 in Hz// channel width
  sscanf(argv[8],"%lf",&chan0); //1
  sscanf(argv[9],"%d",&Nchan); //can vary  8/20/64
  printf("%d %lf %ld %lf %lf %lf %lf %lf %d \n",NN,LL,seed,ff,spindex,nu0,deltanu,chan0,Nchan);
   
 
  pixel=LL;// pixel resolution in arcmin, we have used this in HEADER
 
  LL=M_PI*LL/(180.*60.);// in rad, added
 
   
  
  //memory allocation for out and in  NN*(NN+2) array 

  out=(double*)calloc ((NN*(NN+2)),sizeof(double));
  img=(double*)calloc (NN*NN*Nchan,sizeof(double));
  // call the functions
  
  Fill_Modes(NN, LL, seed,(*P_T),out);
  
  Make_Image(NN,LL,out);

  Non_Gauss(NN,LL,out,ff);
   
  Write_Fits_Image(NN,LL,out,argv[10]); //In temerature unit at central freq
 
  PB_Flux_MF(NN, LL, out, img, Nchan, nu0, chan0, spindex, deltanu);// multi-frequency
 printf("Fine up to this."); 
  Write_Fits_MF_Image(NN, LL, Nchan, nu0,  chan0, deltanu, img, argv[11]); //In Intensity


 printf("Writing Fits properlys."); 

  
  //----------- Following is for validation -------------//
  
  //Get_Modes(NN,LL,out); 
  
  
  //length=LL*NN; // angular dimension of one edge of image in radian
  //############### for binned power spectrum calculation #################### //
  int i,NB=15,*num; 
  double *binu,*binps;
 
  num=(int*)calloc(NB,sizeof(int)); // stores number of U modes in each bin
  binu=(double*)calloc(NB,sizeof(double));// stores sum of |U| in each bin
  binps=(double*)calloc(NB,sizeof(double));// PS in each bin
    
  //power_spec(NN,LL,out,NB,binu,binps,num);// Calls the power spectrum function
  power_spec(NN,LL,img,NB,binu,binps,num);// Calls the power spectrum function

  //Write power spectrum in data file
  fp=fopen("psout.dat","w");

  for(i=0;i<NB;++i)
    if(num[i]>0)
      {
	fprintf(fp,"%e %e %e %d \n",binu[i],binps[i],P_T(binu[i]),num[i]); // These are in temperature unit
      }
  fclose(fp);
  //############################################################
  
  free(binu);
  free(num);
  free(binps);
}				  

