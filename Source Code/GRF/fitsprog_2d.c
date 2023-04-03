//This programme for FITS HEADER
#include <stdio.h>
#include <fitsio.h>
#include <math.h>

void printerror(int status){
  if (status){
    fits_report_error(stderr, status);
    exit( status ); 
  }
}

void Write_Fits_Image(int NN,double LL,double *out,char *outfile)
{
  fitsfile *fitsp;
  int ii,jj,index,index1,status=0,naxis=2;
  //double nu0;

  // for header writing 
 //Writing in FITS format
  //Header information 4 parameters
  double CRVAL[naxis],CDELT[naxis],CRPIX[naxis];
  long naxes[naxis];
  double pixel= (180.*60*LL)/M_PI;// in acrmin
 
  char *CTYPE[2] = {"thetax", "thetay"};
  char keynam[6];

 // This is for single channel. If you want to make multi-chanel , then // terms (blue colour) are to be added
  
  naxes[0]=(long)NN;naxes[1]=(long)NN; //naxes[2]=(long) Nchan; dimension of each axis
  CRVAL[0]=0.;CRVAL[1]=0.;//CRVAL[2]=nu0; // first pixel postion   (theta),central frequency  in CENTRAL IMAGE
  CDELT[0]=1.*pixel;CDELT[1]=1.*pixel; // CDELT[2]=deltanu;pixel size in arcmin and los direction in MHz i.e (CW)
  CRPIX[0]=(1.*NN)/2.+1;CRPIX[1]=(1.*NN)/2.+1; // CRPIX[2]=chan0;center of image (pixel no.) 



  fits_create_file(&fitsp, outfile, &status);
  fits_create_img(fitsp, DOUBLE_IMG, naxis, naxes, &status);

 
  for(ii=0; ii<naxis; ii++){
    
    if(fits_make_keyn("CTYPE", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_str(fitsp, keynam, CTYPE[ii], " axis type", &status))
      printerror(status);
    
    if(fits_make_keyn("CRVAL", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRVAL[ii], 10, " ", &status))
      printerror(status);
    
    if(fits_make_keyn("CDELT", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CDELT[ii],  9, " ", &status))
      printerror(status);
     
    if(fits_make_keyn("CRPIX", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRPIX[ii],  9, " ", &status))
      printerror(status);
  }

  // header complete 

  double *img; // for storing image 
  img=(double*)calloc(NN*NN,sizeof(double));

  long fimpixel[2],limpixel[2];
  fimpixel[0]=1;fimpixel[1]=1;//first pixel co-ordinate (lower left corner--read from the FITS image )
  limpixel[0]=(long)NN;limpixel[1]=(long)NN;//last pixel co-ordinate (upper right corner--read from the FITS image )


  for(ii=0;ii<NN;++ii)
    for(jj=0;jj<NN;++jj)
      {
	index=ii*(NN+2)+jj;
	index1=jj*NN+ii;
	img[index1]=out[index];// (NN*(NN+2)) complex array to image of NN*NN
      }

  fits_write_subset(fitsp,TDOUBLE,fimpixel,limpixel,img,&status);

  fits_close_file(fitsp, &status);
}//void Write_Fits_Image



void Write_Fits_MF_Image(int NN,double LL, int Nchan, double nu0,double chan0, double deltanu, double *img,char *outfile)
{
  fitsfile *fitsp;
  int ii,jj,kk,index,index1,status=0;
  double nu;//In Hz
 

  // for header writing 
 //Writing in FITS format
  //Header information
  int naxis=3;
  double CRVAL[naxis],CDELT[naxis],CRPIX[naxis];
  long naxesim[naxis], naxes[naxis];
  double pixel= (180.*60*LL)/M_PI;// inarcmin
  
  
 
 
  char *CTYPE[4] = {"thetax", "thetay","3IMAGE", "CHAN"};
  
  char keynam[6];

 // This is for multi channel. If you want to make multi-chanel , then // terms (blue colour) are to be added
  
  naxesim[0]=(long)NN;naxesim[1]=(long)NN, naxesim[2]=(long) Nchan;// dimension of each axis// image dimension : NXNXNchan
  CRVAL[0]=0.;CRVAL[1]=0.,CRVAL[2]=nu0; // CENTRAL image center position (theta),central frequency 
  CDELT[0]=1.*pixel;CDELT[1]=1.*pixel, CDELT[2]=deltanu;//pixel size in arcmin
  CRPIX[0]=(1.*NN)/2.+1;CRPIX[1]=(1.*NN)/2.+1, CRPIX[2]=chan0;//center of the multi-channel-image (pixel no.) 

 
  fits_create_file(&fitsp, outfile, &status);
 
  fits_create_img(fitsp, DOUBLE_IMG, naxis, naxesim, &status);

  
  long fimpixel[3],limpixel[3];

  for(ii=0; ii<naxis; ii++){
    
    if(fits_make_keyn("CTYPE", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_str(fitsp, keynam, CTYPE[ii], " axis type", &status))
      printerror(status);
    
    if(fits_make_keyn("CRVAL", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRVAL[ii], 10, " ", &status))
      printerror(status);
    
    if(fits_make_keyn("CDELT", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CDELT[ii],  9, " ", &status))
      printerror(status);
     
    if(fits_make_keyn("CRPIX", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRPIX[ii],  9, " ", &status))
      printerror(status);
  }

  // header complete 
  
 
  double *img1; // for storing image 
  img1=(double*)calloc(NN*NN,sizeof(double));
  
  //first pixel co-ordinate (lower left corner--read from the FITS image )
 	
  // fits_open_file(&fptrim,argv[11],READWRITE,&status);

  for(kk=0;kk<Nchan;kk++)  // chanel loop starts
    {
      nu=nu0+(kk+1.+0.5)*deltanu;
      fimpixel[0]=1;fimpixel[1]=1;fimpixel[2]=(long) (1+kk);
      limpixel[0]=(long)NN;limpixel[1]=(long)NN;limpixel[2]=(long)(1+kk);
      // printf("chan=%d nu=%e image=[%ld %ld %ld] [%ld %ld %ld]\n",kk,nu,fimpixel[0],fimpixel[1],fimpixel[2],limpixel[0],limpixel[1],limpixel[2]);
     
  for(ii=0;ii<NN;++ii)
    for(jj=0;jj<NN;++jj)
      {
	index=kk*NN*NN + ii*NN +jj;
	index1=ii*NN+jj;
	img1[index1]=img[index];
      }

  fits_write_subset(fitsp,TDOUBLE,fimpixel,limpixel,img1,&status);
    }
  fits_close_file(fitsp, &status);
}
