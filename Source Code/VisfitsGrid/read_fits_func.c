//UVFITS reading functions: This program is used to read UVFITS file: 
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <fitsio.h>
# include <fftw3.h>
//# include <nr.h>
//# include <nrutil.h>
# include "read_fits_func.h"

static fitsfile *fptr;
long el1,nel; //element number 1, total number of elements
int status,anynul;
char object[100],obs_date[100],instru[100],telescope[100];
float ra,dec,pscu,pscv; // RA and DEC of target source
  
// extern in main progr visfitsgrid.c
long gcount,pcount,nstokes,nchan,ncmplx;     // pcount-8 number of random group parameters 
                                             // gcount- number of random groups 
float *data;                                 // Nchan, number of chanels to be anylysed
float *randpar; 
float  chan0,nu_chan0,del_chan;
float variancer=0.,variancei=0.;
float meanr=0.,meani=0.;
float *uu,*vv,**RRre,**RRim,**LLre,**LLim;
int flagLL;
FILE *fp, *fp1;

int read_fits_header(char* in_file)
{
  char key_simple[FLEN_KEYWORD]="SIMPLE"; //T- file does hold to FITS standard 
  char key_naxis[FLEN_KEYWORD]="NAXIS"; //32/64 /bits per data pixel
  char key_naxis2[FLEN_KEYWORD]="NAXIS2"; //data axis 1
  char key_naxis3[FLEN_KEYWORD]="NAXIS3"; //data axis 2
  char key_naxis4[FLEN_KEYWORD]="NAXIS4"; //data axis 3
  char key_gcount[FLEN_KEYWORD]="GCOUNT";
  char key_pcount[FLEN_KEYWORD]="PCOUNT";
  char keyvalue[FLEN_VALUE];
  char comment[FLEN_COMMENT];
  
  int ii,jj,kk,st,ttest;

  status=0;
  anynul=0;
  printf("-----------Start to read  UV - Fits header-------------------\n");
  fits_open_file(&fptr, in_file, READONLY, &status); // open uvfits file,
  if(status==1){fprintf(stderr,"unable to open %s", in_file);	return 1;}// input file pointer is fptr

  fits_read_keyword(fptr,key_simple,keyvalue,comment,&status);
  if(*keyvalue!='T'){ fprintf(stderr,"Input File is NOT a SIMPLE FITS file\n");return 1;}

  if(fits_read_key_lng(fptr,key_naxis2,&ncmplx,comment,&status))
    printerror( status );
  printf("NAXIS2=%ld\n",ncmplx);

  if(fits_read_key_lng(fptr,key_naxis3,&nstokes,comment,&status))
    printerror( status );
  printf("NAXIS3=%ld\n",nstokes);
  
  if(fits_read_key_lng(fptr,key_naxis4,&nchan,comment,&status))
    printerror( status );
  printf("NAXIS4=%ld\n",nchan);
  
  if(fits_read_key_lng(fptr,key_gcount,&gcount,comment,&status))
    printerror( status );
  printf("GCOUNT=%ld\n",gcount);
  
  if(fits_read_key_lng(fptr,key_pcount,&pcount,comment,&status))
    printerror( status );
  printf("PCOUNT=%ld\n",pcount);
  
  if(fits_read_key_flt(fptr, "CRVAL4", &nu_chan0, comment,&status))
    printerror( status );
  printf("CRVAL4=%e\n",nu_chan0); //CRVAL4-- frequency of first channel- 150 MHz
  
  if(fits_read_key_flt(fptr, "CDELT4", &del_chan, comment,&status))
    printerror( status );
  printf("CDELT4=%e\n",del_chan);// CDELT4-- channel width 62.5 KHz
  
  if(fits_read_key_flt(fptr, "CRPIX4", &chan0, comment,&status))
    printerror( status );
  printf("CRPIX4=%f\n",chan0);// CRPIX4- Starting, eg. channel here 65, Total Channel 128

//ref. pixel; CRPIXn can be a floating point no. x for which the physical value is CRVALn & increment CDELTn . We should interpret CRPIXn as the location of a FORTRAN index. For any given channel i the freq. is nu[i]=CRVAL4 + CDELT4*(i+1-CRPIX4).
  
  if(fits_read_key(fptr, TSTRING,"OBJECT",&object, comment,&status))
    printerror( status );
  printf("OBJECT=%s\n",object);
  
  if(fits_read_key(fptr, TSTRING,"DATE-OBS",&obs_date, comment,&status))
    printerror( status );
  printf("DATE-OBS=%s\n",obs_date);
  
  if(fits_read_key(fptr, TSTRING,"TELESCOP",&telescope, comment,&status))
    printerror( status );
  printf("TELESCOP=%s\n",telescope);
  
  if(fits_read_key(fptr, TSTRING,"INSTRUME",&instru, comment,&status))
    printerror( status );
  printf("INSTRUME=%s\n",instru);
  
  if(fits_read_key(fptr, TFLOAT,"CRVAL6",&ra, comment,&status))
    printerror( status );
  printf("RA=%e\n",ra);
  
  if(fits_read_key(fptr, TFLOAT,"CRVAL7",&dec, comment,&status))
    printerror( status );
  printf("DEC=%e\n",dec);
  
  if(fits_read_key(fptr, TDOUBLE,"PSCAL1",&pscu, comment,&status))
    printerror( status );
  printf("PSCAL1=%e\n",pscu);
  if(fits_read_key(fptr, TDOUBLE,"PSCAL2",&pscv, comment,&status))
    printerror( status );
  printf("PSCAL2=%e\n",pscv);
 
  printf("chan0=%f\t nu_chan0=%f\n",chan0,nu_chan0);
  
  // for reading random parameters , memory allocation
  randpar=(float*) calloc((pcount), sizeof(float));

 
  
  // for reading visibility data 
  nel=ncmplx*nstokes*nchan; //Total number of elements

  if((data = (float*)calloc(nel, sizeof(float)))==NULL)// memory allocation for data
      { fprintf(stderr,"malloc failure\n"); return 1;}
  
  return 0;
  
}

//----------------------------------------------
//============reads the random group parameters from uv data file=======================//
//---------------------------------------------

int read_ranpar(long grp)
{
  status=0;
  el1=1;  
  if(fits_read_grppar_flt(fptr,grp,el1,pcount,randpar,&status))
  printerror( status );
  return 0;
 }

//============reads the random group parameters from uv data file==== Done ===================//


//----------------------------------------------
// reads the visibilities from uv data file
//---------------------------------------------
int read_data(long group)
{
  float nulval;

  status=0;
  anynul=0;
  nulval=0.;
  el1=1;

  if(fits_read_img_flt(fptr,group,el1,nel,nulval,data,&anynul,&status)) // Image fits is not allowed
  printerror( status );
  return 0;
}

int close_fits()//close the input fits file
{
  if(fits_close_file(fptr,&status))
    printerror( status );
  return(0);   
}


//===================================== Done===========================================================

// After reading UVFITS, it finds the statistics of visibility data set(for given UVFITS file)

int readfits(char* in_file,double Umax,double Umin,int chan1,int chan2) 
{
  FILE *fp2, *fp1;
  int n_ave,ii,nbasln=0,flag=0;
  long group;
  float Uval,delnubynu,signv;
  unsigned long total=0;

  typedef struct visibility_type {float re,im,wt;} visibility;  // to read the visibility data
  visibility *v,*RR,*LL;  // The outputs are re, im and wt of the visibility data
  
  read_fits_header(in_file); 
  printf("------------------End UV Fits Header---------------\n\n");
  
  //==================================== Header reading  Done==========================//

  

  n_ave=(chan2-chan1+1); //number of chanels: Echan (chan2)-Bchan(chan1)+1: done

  delnubynu=del_chan/nu_chan0; //done, reads from UVFITS header
  
  // Memory allocation for following parametes
 
  uu = (float*)calloc(gcount,sizeof(float));
  vv = (float*)calloc(gcount,sizeof(float));
  RRre = (float**)calloc(gcount,sizeof(float*));
  RRim = (float**)calloc(gcount,sizeof(float*));
  RRre[0]=(float*)calloc(gcount*n_ave,sizeof(float));
  RRim[0]=(float*)calloc(gcount*n_ave,sizeof(float));
  LLre = (float**)calloc(gcount,sizeof(float*));
  LLim = (float**)calloc(gcount,sizeof(float*));
  LLre[0]=(float*)calloc(gcount*n_ave,sizeof(float));
  LLim[0]=(float*)calloc(gcount*n_ave,sizeof(float));
  
  for(ii=1;ii<gcount;++ii)
    {
      RRre[ii]=RRre[0]+ii*n_ave; //Real part of RR ploarization of the signal/ visibility/data
      RRim[ii]=RRim[0]+ii*n_ave;// Imaginary part of RR polarization of the signal/ visibility/data
      LLre[ii]=LLre[0]+ii*n_ave; //Real part of LL ploarization of the signal/ visibility/data
      LLim[ii]=LLim[0]+ii*n_ave; //Imaginary part of LL polarization of the signal/ visibility/data
    }

  
  for(group=1;group<=gcount;group++) //baseline loop
   
    { 
      read_ranpar(group);//random group parameters u,v
      randpar[0]*=nu_chan0;
      randpar[1]*=nu_chan0;
      
      //Determine U value
//      Uval=sqrt(SQR(1.*randpar[0])+SQR(1.*randpar[1]));//multiply PSCAL1*CRVAL4 for GMRT data     // Changed by Sri
	Uval=sqrt((1.*randpar[0]*randpar[0])+(1.*randpar[1]*randpar[1]));//multiply PSCAL1*CRVAL4 for GMRT data                // Changed by Sri
      
      
      if(Uval>=Umin && Uval<=Umax) 
     	{
	  read_data(group);// Reads  data for  group: reads the visibilities from uv data file

	  signv= (randpar[1]<0.) ? -1. : 1. ;// change sign if negetive
	  
 
	  //  uu1 = signv*randpar[0]; //added
	  //  vv1 = signv*randpar[1]; //added
     
	  //  printf("%e\t%e\n", (153.*uu1)/157., (153.*vv1)/157. );// change 157 MHz to 153 MHz
     

	  
	   uu[nbasln]=signv*randpar[0];
	   vv[nbasln]=signv*randpar[1];
	  

	  
	  
	  v=(visibility *)data;
	  v+=(chan1)*nstokes;// goes to the starting channel
	  RR=v;
	  LL=v+(nstokes-1);
	  
	  flag=0;
	  
	  for(ii=0; ii<n_ave;ii++) // channel loop
	    {
	      if(RR->wt>0.) // for RR
		{
		  RRre[nbasln][ii]=RR->re;
		  RRim[nbasln][ii]=RR->im;
		 
		 
		  		  
		  meanr+=RRre[nbasln][ii];
		  variancer+=pow(RRre[nbasln][ii],2.);
		  meani+=RRim[nbasln][ii];
		  variancei+=pow(RRim[nbasln][ii],2.);
		 
		  flag=1;
		  total++;	      /* done mean and variance of  RR */
		  
		  /* complex conjugate if v <0 */
		  RRim[nbasln][ii]=signv*RRim[nbasln][ii];
		}
	      else // the data is flagged
		{
		  RRre[nbasln][ii]=-1.0e8;// why not =0.;
		  RRim[nbasln][ii]=-1.0e8; // why not =0.;
		}
		// for LL			
	      if(LL->wt>0. && flagLL==1)
		{
		  LLre[nbasln][ii]=LL->re;
		  LLim[nbasln][ii]=LL->im;
		
		  
		  
		  meanr+=LLre[nbasln][ii];
		  variancer+=pow(LLre[nbasln][ii],2.);
		  meani+=LLim[nbasln][ii];
		  variancei+=pow(LLim[nbasln][ii],2.);
		  
		  flag=1;
		  total++;	      /* done mean and variance combing RR and LL*/
		  
		  /* complex conjugate if vv <0 */
		  LLim[nbasln][ii]=signv*LLim[nbasln][ii];
		}
	      else // the data is flagged
		{
		  LLre[nbasln][ii]=-1.0e8;//fine with arguments in vcg.c
		  LLim[nbasln][ii]=-1.0e8; 
		}
	      
	      v+=nstokes;// point to the next channel
	      RR=v;
	      LL=v+(nstokes-1);
	     
	    } // for(ii=0; ii<n_ave;ii++) channel loop
	  
	  if(flag==1){ ++nbasln;}
	
	}  // if(Uval>=Umin && Uval<=Umax) 
      //fclose(fp1);  
    }// for(group=1;group<=gcount;group++)
  close_fits();
  
  /* ===========finished putting data into array ========Done=======*/
  
   // Real Mean of visibility = [RR(real) + LL(real)]/2 
   // Imaginary  Mean of visibility = [RR(ima) + LL(ima)]/2 
/* print mean and rms */
  meanr= meanr/(1.*total);
  variancer=variancer/(1.*total)-pow(meanr,2.);
  meani= meani/(1.*total);
  variancei=variancei/(1.*total)-pow(meani,2.);
  
  printf("total= %ld\n ",total);
  printf("total baselines= %d\n ",nbasln);
  
  printf("total= %ld \n mean (re) = %e Jy rms. (re)  = %e Jy\n",total,meanr/2.,sqrt(variancer/2.));
  
  printf("mean (im) = %e Jy rms. (im)  = %e Jy\n",meani/2.,sqrt(variancei/2.));
  
  return(nbasln);
}

