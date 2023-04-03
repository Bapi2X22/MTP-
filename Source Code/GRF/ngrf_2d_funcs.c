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

extern double k_B, c;



//Function to calculate the Angular power spectrum from Fourier modes

void power_spec(int NN,double LL,double *out,int NB,double *usum,double *pksum,int *no)
// NN -- image dimension, LL -- grid spacing, *out -- array of size NN*(NN+2) stores image or Fourier components, NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin 
   
{
  
  fftw_complex *in;  // stores Fourier component of temp fluct
  in=(fftw_complex*)&out[0];
 
  // local variables 

  double Umax,Umin,binsiz,uu,Length,Areainv;
  int i,j,index,ia,tmp;

  Length=(NN*LL); //in radian
  Areainv=pow(NN*LL,-2.0); //1/solid angle

  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }

  Umax=(NN/2.)/Length;//3439.5
  Umin=1./(Length);
  binsiz=(log10(1.*NN/2.)/(1.*NB));
  
  for(i=0;i<NN;++i)
    for(j=0;j<NN/2+1;++j)
      {
  	index=i*(NN/2+1)+j; //got it, column major
  	ia=(i>NN/2) ? (NN-i) : i ; //symmetric indentity
  	uu=sqrt(1.*(ia*ia+j*j))/Length;

  	if((uu>=Umin)&&(uu<Umax))
  	  {
  	    tmp=(int) floor(log10(uu/Umin)/binsiz);
	    tmp=(tmp<NB) ? tmp: tmp-1;
  	    no[tmp]++;
	    pksum[tmp]+=((in[index][0]*in[index][0])+(in[index][1]*in[index][1])); // X*X + Y*Y: numarical angular power spectrum, see paper samir 2014
  	    usum[tmp]+=uu;
  	  }
      }
  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]=usum[i]/no[i];
	pksum[i]=Areainv*pksum[i]/no[i]; // numarical  angular power spectrum
      }
}



// Given input brightness temperature power spectrum, fill modes assuming Gaussian random field 
void Fill_Modes(int NN,double LL,unsigned long int seed,double (*func)(double),double *out)
// NN -- image dimension, LL -- grid spacing, seed - for random number generation, (*func)(double) -- function returning input power spectrum, *out -- array of size NN*(NN+2) stores image or Fourier components

{  
  fftw_complex *in;
  in=(fftw_complex*)&out[0];
  
  // local variables 
  double Length,  uu,amp, fact;
  int i,j,ia,index,index1,xdim,ydim;
  
  //for random number generator
  gsl_rng *r;
  double sigma=1.;

  r= gsl_rng_alloc(gsl_rng_cmrg);
  gsl_rng_set (r, seed);
  //done

  //initalize local  variables 

  Length=NN*LL;
  fact=Length/sqrt(2.0); // factor for amplitude
  ydim=(NN/2+1); 
  xdim = NN;
  //Filling Fourier Components  

  //along axis (j-0 and j=N/2)

  for(j=0;j<ydim;j=j+NN/2)  // choosing y = 0 and y= NN/2 i.e. fills the x-axis only
    for(i=1;i<xdim/2;++i) // avoid filling end points twice, that's why starting from 1 and ending before N/2
      {
	// along + x (positive x)

	uu=sqrt(1.*(i*i+j*j))/Length;
	amp=fact*sqrt(func(uu)); // (func(uu)=P_T(double UU)
	
	index=i*ydim+j;  //ydim is period// column major, periodiciy is row
	
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma); // Tilda T(U) = X + iY, X =  in[index][0]
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma); // iY= in[index][1]
	
	
       
	// along  -x (negetive x)
	index1=(NN-i)*ydim+j;
	
	in[index1][0]=in[index][0]; //X
	in[index1][1]=-1.0*in[index][1];// using (-), putting complex conjugate (-iY) on other half. 
      }
  // printf("Random numbers %e",in[0][0]);
  
  
  // upper half plane excluding x axis
  
  for(i=0;i<xdim;++i)
    for(j=1;j<NN/2;++j)
      {
	ia= (i>NN/2) ? (NN-i) : i ;
	uu=sqrt(1.*(ia*ia+j*j))/Length;
	amp=fact*sqrt(func(uu));                // amplitude factor for Fourier modes calculation
	
	index=i*ydim+j;// column major
	
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
      }
      
  //4 points remain// (0,0) & (N/2,N/2) (re + im) 

  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      {
	if(i+j==0)  //(0,0)
	  {
	    in[0][0]=0.0;
	    in[0][1]=0.0;
	  }
	else    // (NN/2,NN/2)
	  {
	    uu=(NN/2.)*sqrt(1.*(i*i+j*j))/Length; //baseline calculation using grid index
	    amp=fact*sqrt(func(uu));// amplitude factor for Fourier modes calculation in tilde T(U)
	    
	    index=i*(NN/2)*ydim+j*(NN/2);
	    
	    in[index][0]=pow(-1.,(i*NN/2+j*NN/2))*amp*gsl_ran_gaussian(r,sigma);
	    in[index][1]=0.0; 
	  }
      }
 
  gsl_rng_free(r);
  // finished filling Fourier components
}

// Complete 


// Fourier Modes -> Image

void Make_Image(int NN,double LL,double *out) // image without beam // 
// NN*NN -- image dimension, LL -- grid spacing, *out -- array of size NN*(NN+2) stores image or Fourier components
{
  // local variables
  int i,j,index;
  double Areainv;

 
  fftw_complex *in; // stores Fourier component of temp fluct
  in=(fftw_complex*)&out[0];
  
  
  // fft plans: in => \Delta T(Ux,Uy), out = \deltaT(thetax,thetay) 
  fftw_plan p;
  p= fftw_plan_dft_c2r_2d (NN, NN, in,out, FFTW_ESTIMATE);// complex to real
  fftw_execute(p); // perform FT to make image in----> out

  Areainv=pow(NN*LL,-2.0); // 

  for(i=0;i<NN;++i)
    for(j=0;j<NN;++j)
      {
	index=i*(NN+2)+j; //column major
	out[index]*=Areainv; 
	//	printf("\n(i = %d \t,j = %d,\t index=%d)",i, j, index);
      }
 
  fftw_destroy_plan(p);
}
// Finished


// Given Gaussian random image, this introduces non-Gaussianity 
// d=d + ff*(d*2-var)/sqrt(var) 
void Non_Gauss(int NN,double LL,double *out,double ff)
// NN -- image dimension, LL -- grid spacing, *out -- array of size NN*(NN+2) stores image or Fourier components,, dimensionless non-Gaussianity parameter
{
  double var=0; // ff non-Gaussianity parameter, var - variance of Gaussian field
  double mean=0;
  
  int i,j,index;
  // calculate variance 

  for(i=0;i<NN;++i)
    for(j=0;j<NN;++j)
      {
	index=i*(NN+2)+j;
	var+=pow(out[index],2.);

      }
  var/=(1.*NN*NN);

  printf("sqrt(var)=%e\n",sqrt(var));
   
  // non-Gaussianity introduced 
  for(i=0;i<NN;++i)
    for(j=0;j<NN;++j)
      {
	index=i*(NN+2)+j;
       	out[index]+=ff*(pow(out[index],2.)-var)*pow(var,-0.5);// non-Gaussianity introduced in the image   	
	mean+= out[index];
      
      }
  mean/=(1.*NN*NN);
  
  printf("mean of NG field=%e\n",mean);

  // done 
}

//================ Multi-frequency  Flux(nu,theta)*A(nu,theta)==========  Starts===//


void PB_Flux_MF(int NN,double LL, double *out, double *img, int Nchan, double nu0,double chan0, double spindex,double deltanu)
//nu--dummu Central frequency
{
  double thetax, thetay, theta, rjfac, LLsq, nu, be; 
  int i,j,k, index,index1;
  double Beam(double theta,double freq);
  long fimpixel[3],limpixel[3];  //3D   

  LLsq=pow(LL,2.0);  
    
  for(k=0;k<Nchan;k++)  // chanel loop starts
    {
      nu=nu0+(k+1.-chan0)*deltanu;
      //nu=nu0; // no freq. dependence.
      
      fimpixel[0]=1;fimpixel[1]=1;fimpixel[2]=(long) (1+k); //3D
      limpixel[0]=(long)NN;limpixel[1]=(long)NN;limpixel[2]=(long)(1+k); //3D
 
      rjfac = 2.*k_B*nu0*nu0/(c*c); 
  
      for(i=0;i<NN;++i)
	for(j=0;j<NN;++j)
	  {
	    index=i*(NN+2)+j;
	    index1=NN*NN*k+i*NN+j;//column major
	    thetax= (i-NN/2)*LL;
	    thetay= (j-NN/2)*LL;
	    theta=sqrt((thetax*thetax)+(thetay*thetay));
	    be=Beam(theta,nu0);
	    img[index1]=out[index]*(LLsq*rjfac*be)*pow((nu0/nu),spindex); // in function : Beam = beam*beam
	  }

    }
 
}

// =================Image --> Fourier Modes=========================//

void Get_Modes(int NN,double LL,double *out)
// NN -- image dimension, LL -- grid spacing, *out -- array of size NN*(NN+2) stores image or Fourier components: Visibilities generation
{

  // local variables
  int i,j,index;
  double LLsq;

  // for filling Fourier components 
  fftw_complex *in;
  in=(fftw_complex*)&out[0];
  // we have generated visibilities for given APS and they are stored in the "in" array: see below
  fftw_plan p1;
  p1= fftw_plan_dft_r2c_2d (NN, NN, out,in, FFTW_ESTIMATE);
  fftw_execute(p1); // perform FFT on image, real to complex, gives an array "in"

  LLsq=pow(LL,2.0);

  for(i=0;i<NN;++i)
    for(j=0;j<NN;++j)
      {
	index=i*(NN+2)+j;// row-major  
	out[index]*=LLsq; //  specific intensity to flux 
      }
  fftw_destroy_plan(p1);
}
// Finished
