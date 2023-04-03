// This program reads a multichannel FITS image (diffuse) and UVFITS template then, it  interpolates the gridded visibilities for each sample baseline (from .FITS image after FFTW) to the nearest baseline of the uv track genereted from  the GMRT simulation
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>
# include <fftw3.h>
# include "read_fits_func.h"
# include "fitsprog.h"

extern double CRVAL[],CDELT[],CRPIX[];
extern long naxes[],naxesim[];

extern float del_chan,nu_chan0,chan0,Umaxc;
extern long nstokes,nchan,ncmplx,gcount;

int main(int argc, char *argv[])
{ 
  char IMAGEFITS[128],UVFITS[128];
  int chan,n_ave,ii,jj,ii1,jj1;
  int nbasln,NN;
  long index,index1;
  double dU,fac;
  
 
  if(argc!=4)
    {
      printf("Usage: %s <input FITS image file> <input UVFITS file> <factor=0(overwrite, 1=add)>\n", argv[0]); // we use 1
      return 1;
    }
  
  sscanf(argv[1],"%s",IMAGEFITS);
  sscanf(argv[2],"%s",UVFITS);
  fac=atof(argv[3]);

  if(access(IMAGEFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",IMAGEFITS);
    exit(0);
  }

  if(access(UVFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",UVFITS);
    exit(0);
  }
  //===================Read Image Header and print them  ==========================//
  printf("\nReading Simulated Images... (fits)\n");
  
  printf("Image Header\n");
  
  SREAD_HDR(IMAGEFITS); // it is written in fitsprog.c which is used to read image  header. It acts through fitsprog.h.
  
  NN=naxesim[0];n_ave=naxesim[2]; //naxesim[2]=(long) Nchan, defined in fitsprog_2d.c
  
  dU=(double) ((180.*60.)/(M_PI*NN*CDELT[0])); // Umin, arcmin -radian, du= 1/(LL*NN)(acrmin)= =3.85 ardaian, In image CDELT[0] is in arcmin

  printf("pixel=LL=%lf arcmin, naxes image={%ld,%ld,%ld}\ndU=%e\n",CDELT[0],naxesim[0],naxesim[1],naxesim[2],dU);
  printf("NN=%d n_ave=%d\n",NN,n_ave);

  //===================Read Image Header and print Done  ==========================//
  


  //====================Following parts will read image and FFT and fill in a array ==================================//
 

  fitsfile *fptr;
  long fpixel[naxis],lpixel[naxis],inc[naxis];
  int status=0,anynull;
  double nulval=0;
  fftw_complex *in,*reim;
  double *out,*img;
  fftw_plan p1;

  reim=(fftw_complex*)calloc ((NN*(NN/2+1)*n_ave),sizeof(fftw_complex));
  out=(double*)calloc ((NN*(NN+2)),sizeof(double));
  in=(fftw_complex*)&out[0];
  img=(double*)calloc(NN*NN,sizeof(double));//image
 
  p1= fftw_plan_dft_r2c_2d (NN, NN, out,in, FFTW_ESTIMATE);// FFT on real -plan declare

  fits_open_file(&fptr,IMAGEFITS,READONLY,&status); //Image fits file is openning 
  
  for(chan=0;chan<n_ave;++chan) //chanel loop starts
    {
      fpixel[0]=1;fpixel[1]=1;fpixel[2]=(long)(chan+1); //first pixel co-ordinate for a channel (given by channel loop)
      lpixel[0]=NN;lpixel[1]=NN;lpixel[2]=(long)(chan+1); //last pixel co-ordinate
      
      inc[0]=1;inc[1]=1;inc[2]=1;// increment
      
      fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,&nulval,img,&anynull,&status); //read fits-image  channelwise one by one (subset): first channel
      
      for(jj=0;jj<NN;++jj)
  	for(ii=0;ii<NN;++ii)
  	  {
  	    index=jj*NN+ii; 
  	    index1=ii*(NN+2)+jj;
  	    out[index1]=img[index]; //1st channel image :transfer image[index] to out[index1]
  	  }
      // continuing channel loop
     
      fftw_execute(p1);//FFT on Real image-out[index1] ----> in[index1] [Gives visibilities (complex)]
      

      for(ii=0;ii<NN;++ii) //visibilities are copied from "in" to "reim"
  	for(jj=0;jj<=NN/2;++jj)
  	  {
  	    index1=chan*(NN/2+1)*NN+ii*(NN/2+1)+jj;// channel wise visibilityies are stacked and stored through first part of index1
  	    index=ii*(NN/2+1)+jj;
  	    reim[index1][0]=pow(-1.,ii+jj)*in[index][0];//Image centre at (N/2,N/2) //real part [0]
  	    reim[index1][1]=pow(-1.,ii+jj)*(-1.*in[index][1]);//r2c use -1, to make it +1 exponent multiply with -1, imaginary part[1], index1-reffers to channel
  	  }
    }//chanel loop stops
  fftw_destroy_plan(p1);
  fits_close_file(fptr,&status);
  

//======================================read image and FFT and fill in a multichannel array======Done==================//
  



//========Only to read UVFITS header(switch on to  read_fits_func.c)=================)//
 

  nbasln=readfits(UVFITS,0.,0.,1.,1.);// The following function is used to understand the arguments here, hence it will not work to find the statistics for the given argument
 
//int readfits(char* in_file,double Umax,double Umin,int chan1,int chan2
  

//=======================put grid value in UV track===== Starts =================//
  

  float randpar[3],*data;
  long group;
  long nel,el1=1;  //nel-- no. of element
  double uuc,vvc,signv;
  int stokes,anynul=0;
  float nulval1=0;
  double *corrfact;

  nel=ncmplx*nstokes*nchan;
  

  printf("ncmplx=%ld nstokes =%ld nchan=%ld nel=%ld gcount=%ld\n",ncmplx,nstokes,nchan,nel,gcount); //To check
  
  printf("n_ave=%d (image) nchan=%ld (vis)\n",n_ave,nchan);//to check
  if(n_ave!=nchan)
    {
      printf("n_ave should be equal to nchan\n");  //To check
      return 1;
    }
   printf("chan0=%lf \n",chan0);//to check
  
  corrfact = (double*) calloc(n_ave, sizeof(double));
  for(ii=0; ii<n_ave; ii++) // channel loop to see correction factor for fre. changing about 1
    {
      //corrfact[ii] = 1.+ (del_chan/nu_chan0)*(ii+1.+0.5-chan0);//ok. lamda[chan0]/lambda[ii]=nu[ii]/nu0
      corrfact[ii] = nu_chan0 +  del_chan * (ii-(chan0-1)); // frequency of ith channel.
    }
 
  data=(float*)calloc(nel,sizeof(float));  //================data memory allocation=================//

  
  fits_open_file(&fptr,UVFITS,READWRITE,&status); //====uvfits is opening for reading and writing u and v values==========//
  
  for(group=1;group<=gcount;group++)// baseline taken from UVFITS loop starts
    {
      fits_read_grppar_flt(fptr,group,el1,3,randpar,&status); //pcount=3
      
      //rand parameter read done  u,v and w   
      
      signv= (randpar[1]<0.) ? -1. : 1. ; // v-- negetive, v*(-1)--> becomes positive, folded to the upper half plane
            
      fits_read_img_flt(fptr,group,el1,nel,nulval1,data,&anynul,&status); //reads data from UVFITS 
      
      for(chan=0;chan<n_ave;chan++)  //===== chanel loop starts====//
	for(stokes=0;stokes<nstokes;++stokes) //polarization loop
	  {
	    index=chan*nstokes+stokes;//indexing on uv plane of initial data
	    data[3*index+0]*=fac;// fac=1, atof(argv[3]);intital data, O*1=0
	    data[3*index+1]*=fac;//
	  }//===== channel loop ends
      
      for(chan=0;chan<n_ave;chan++)// Again channel loop starts
	{
	  uuc=signv*randpar[0]*corrfact[chan]; //freq sacling of u, observed  at central fre.
	  vvc=signv*randpar[1]*corrfact[chan]; // freq scaling of v, v observed at central fre.
	  
	  //uuc=signv*randpar[0]*corrfact[chan]*nu_chan0;
	  //vvc=signv*randpar[1]*corrfact[chan]*nu_chan0;
	  
	  if(abs(uuc)<(NN*dU/2.) && abs(vvc)<(NN*dU/2.)) //  read  uv compo of UVFITS and compare with the~( 6* Umax) value of Image : check (u,v) from UVFITS < Umax (FITS image)
	    {
	      ii1 = (int)roundf(uuc/dU); // nearest u value on  uv-grid indices, set u index
	      ii=(ii1<0) ? NN+ii1 : ii1 ;
	      jj = (int)roundf(vvc/dU); // nearest v value on  uv-grid indices, set v index
	      
  	      index1=chan*(NN/2+1)*NN+ii*(NN/2+1)+jj; // This indexing (index1) contains channel wise  visibilities data after FFTW of MF images: row-major with channel: reffered on the top:index1=chan*(NN/2+1)*NN+ii*(NN/2+1)+jj
  	      
	      for(stokes=0;stokes<nstokes;++stokes)
  		{
  		  index=chan*nstokes+stokes;// sdata index in UVFITS
  		  data[3*index+0]+=reim[index1][0]; //real part of Visibilities[0] are transfered  to UVFITS 
  		  data[3*index+1]+=(signv*reim[index1][1]); //Imaginary part of visibilities[1] are transfered  to UVFIT
		} //stokes loop closed along with Multi-Fre-Image data conveted to multifrequency [index-channel] gridded visibilities in the same input UVFITS
  	    } //nearest u, v  value grid indices
  	} //channel loop
      fits_write_img_flt(fptr, group, 1, nel, data, &status ); //Write data into the  UVFITS
    } //baseline loop
  if(fits_close_file(fptr, &status))
    printerror(status);
}
