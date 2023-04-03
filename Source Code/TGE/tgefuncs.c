# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include <fitsio.h>
# include <time.h>
# include "read_fits_func.h"
# include "tgefuncs.h"

extern long gcount,nstokes,nchan;
extern float *randpar,*data;
extern float del_chan,nu_chan0,chan0; 
extern char INFITS[128],input[128],GV[128];

int nc,Nm,Ng,Nu,chan1=1;        
double Umax,f,Uval,FWHM,dU,theta_0,theta_w,theta_eff;

long group;
double *nua;                   
double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL;


void read_inputs(char *in_file)
{
  FILE *fp;
  fp=fopen(in_file,"r"); 
  if( fscanf(fp,"%d%lf%lf%lf",&nc,&Umax,&FWHM,&f) == 4){
    printf("\nInput parameters: nc=%d, Umax=%.1lf, FWHM=%.1lf arcmin, f=%.1lf",nc,Umax,FWHM,f);
  }
  fclose(fp);
}

double winf(double udif)
{
  double y;
  y=M_PI*pow(theta_w,2.)*exp(-1.*M_PI*M_PI*pow(theta_w,2.)*udif); // udif = (U_g -U_i)^2
  return(y);
}

int initialize()
{
  int ii,jj,nn,sf,Ngrid; 

  read_fits_header(INFITS);  

  theta_0=0.6*FWHM*M_PI/(180.*60.);         
  theta_w=f*theta_0; 
  theta_eff=f*theta_0/sqrt(1+f*f);

  sf = 2; //sampling factor
  dU=sqrt(log(2))/(M_PI*theta_eff*sf); // sampling the effective PB twice.
  Nm=3*sf; //sqrt(ln(500)/ln(2) = 3)

  Ng = (int) ceil(Umax/dU) + Nm; // padding for convolution. 
  Nu = 2*Ng +1;                        
  Ngrid = Nu*Nu;

  nua = (double*)calloc(nc,sizeof(double));
  for(nn=0;nn<nc;nn++)
    {
      nua[nn]= (nu_chan0 + del_chan*(nn - (chan0-1)));
    }
                    
  
  GVRRre = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRim = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLre = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLim = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      GVRRre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVRRim[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLim[ii]=(double*)calloc(nc,sizeof(double)); 
    }

  SCRR = (double**)calloc(Ngrid,sizeof(double*));  
  SCLL = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      SCRR[ii]=(double*)calloc(nc,sizeof(double)); 
      SCLL[ii]=(double*)calloc(nc,sizeof(double)); 
    }

  printf("\nTotal no. of channels (nc) = %d",nc);
  printf("\nFWHM = %.2f arcmin.",FWHM);
  printf("\nGrid spacing, dU = %.2f wavelengths",dU);
}



void group_loop()
{
  float u1,v1,signv;
  int nn,ii,jj,kk,ll,nx,ny,index,ch1,ch2,dnu;
  double diff,wt,gRRr,gRRi,gLLr,gLLi,sRRr,sLLr;
  float *RRre,*RRim,*LLre,*LLim; 
  double *V2RR, *V2LL;
  RRre = (float*)calloc(nc,sizeof(float));
  RRim = (float*)calloc(nc,sizeof(float));
  LLre = (float*)calloc(nc,sizeof(float));
  LLim = (float*)calloc(nc,sizeof(float));
  V2RR = (double*)calloc(nc,sizeof(double));
  V2LL = (double*)calloc(nc,sizeof(double));

  typedef struct visibility_type {float re,im,wt;} visibility; 
  visibility *v,*RR,*LL;
  
  int count=0;
  time_t time_s = 0,time_e,time_d = 0,time_d1=0;
  for(group=1;group<=gcount;group++)
    {    
      read_ranpar(group);
      signv = (randpar[1]<0.) ? -1. :1. ;       // flip uv
      u1 = nu_chan0*signv*randpar[0];
      v1 = nu_chan0*signv*randpar[1];

      Uval=sqrt(u1*u1+v1*v1);

      if(Uval<=Umax) 
	{	  
	  read_data(group);
	  v=(visibility *)data;
	  v+=(chan1-1)*nstokes;
	  RR=v;
	  LL=v+1;
	    
	  for(nn=0;nn<nc;nn++)
	    {	         
	      if(RR->wt>0.)
		{
		  RRre[nn]=RR->re; 
		  RRim[nn]=signv*RR->im;
		}
	      else 
		{
		  RRre[nn]=0.; RRim[nn]=0.;
		}

	      if(LL->wt>0.)
		{
		  LLre[nn]=LL->re; 
		  LLim[nn]=signv*LL->im;
		}
	      else 
		{
		  LLre[nn]=0.; LLim[nn]=0.;
		}
	      v+=nstokes; 
	      RR=v;
	      LL=v+1; 
	    }



	  time(&time_s);
	  for(nn=0;nn<nc;nn++)
	    { 
	      V2RR[nn]=0.; V2LL[nn]=0.;
	    }
	  for(ch1=0;ch1<nc;ch1++)	  
	    for(ch2=ch1;ch2<nc;ch2++)
	      {
		dnu = ch2-ch1;
		V2RR[dnu] += (RRre[ch1]*RRre[ch2] + RRim[ch1]*RRim[ch2]);		
		V2LL[dnu] += (LLre[ch1]*LLre[ch2] + LLim[ch1]*LLim[ch2]);
	      }
	  time(&time_e);
	  time_d += difftime(time_e,time_s);



	  time(&time_s);	  
	  ii = (int) round(u1/dU); 
	  jj = (int) round(v1/dU); 

	  for (kk=-Nm;kk<=Nm;kk++)           
	    for(ll=-Nm;ll<=Nm;ll++)
	      {
		nx = ii+kk;
		ny = jj+ll;
		diff = pow(nx*dU-u1,2.) +  pow(ny*dU-v1,2.);  
		wt = winf(diff);
		    
		for(nn=0;nn<nc;nn++) // gridding & self-correlation.
		  {
		    gRRr = wt*RRre[nn]; 
		    gRRi = wt*RRim[nn];
		    gLLr = wt*LLre[nn]; 
		    gLLi = wt*LLim[nn];
		    sRRr = wt*wt*V2RR[nn];
		    sLLr = wt*wt*V2LL[nn];

		    		  
		    index = (Ng+nx)*Nu + (Ng+ny);		  
		    GVRRre[index][nn] += gRRr; 
		    GVRRim[index][nn] += gRRi;
		    GVLLre[index][nn] += gLLr; 
		    GVLLim[index][nn] += gLLi;
		    SCRR[index][nn] += sRRr;		    	  
		    SCLL[index][nn] += sLLr;		    	  

		    index = (Ng-nx)*Nu + (Ng-ny); // inversion points.
		    GVRRre[index][nn] += gRRr; 
		    GVRRim[index][nn] -= gRRi;
		    GVLLre[index][nn] += gLLr; 
		    GVLLim[index][nn] -= gLLi;
		    SCRR[index][nn] += sRRr;		    	  
		    SCLL[index][nn] += sLLr;		    	  
		  }
	      } 
	  time(&time_e);
	  time_d1 += difftime(time_e,time_s);
	  
	  count += 1;
	}
    }
  printf("\nTime taken in self-correlation: %d hr %d min %d sec\n", (int)(time_d/3600), (int)((time_d/60)%60),(int)(time_d%60));
  printf("\nTime taken in gridding: %d hr %d min %d sec\n", (int)(time_d1/3600), (int)((time_d1/60)%60),(int)(time_d1%60));
  printf("\nThere are %d baselines within the uv range.",count);       
  close_fits();
}

void WriteGVnAC()
{
  int ii,jj,nn, index;
  FILE *fp;
  int MM = 2*Nm;
  printf("\nSize of the output GV: %d*%d*%d.",Nu-2*MM,Nu-2*MM,nc);  

  fp=fopen(GV,"wb"); 
  fwrite(&FWHM,sizeof(FWHM),1,fp);
  fwrite(&dU,sizeof(dU),1,fp);
  fwrite(&f,sizeof(f),1,fp);
  for(ii=0;ii<nc;ii++)
    fwrite(&nua[ii],sizeof(nua[ii]),1,fp);

  for(ii=MM; ii<Nu-MM;ii++)
    for(jj=MM; jj<Nu-MM;jj++)
      {
	index=ii*Nu+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVRRre[index][nn],sizeof(GVRRre[index][jj]),1,fp);
	    fwrite(&GVRRim[index][nn],sizeof(GVRRim[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",GVRRre[index][nn],GVRRim[index][nn]);
	  }
      }
    
  for(ii=MM; ii<Nu-MM;ii++)
    for(jj=MM; jj<Nu-MM;jj++)
      {
	index=ii*Nu+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVLLre[index][nn],sizeof(GVLLre[index][jj]),1,fp);
	    fwrite(&GVLLim[index][nn],sizeof(GVLLim[index][jj]),1,fp);
	  }
      }

  //Self-correlation. (real) 
  for(ii=MM; ii<Nu-MM;ii++)
    for(jj=MM; jj<Nu-MM;jj++)
      {
	index=ii*Nu+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&SCRR[index][nn],sizeof(SCRR[index][jj]),1,fp);
	    fwrite(&SCLL[index][nn],sizeof(SCLL[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",SCRR[index][nn],SCRR[index][nn]);
	  }
      }
  fclose(fp); 
}
