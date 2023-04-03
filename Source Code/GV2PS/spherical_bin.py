#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# # Fourier coefficient matrix (for a single ell-bin)

# In[2]:


# returns A for a cos transform 
def calc_A(Na,Nb):
    A=np.outer(np.arange(Na),np.arange(Nb))
    A=np.cos(np.pi*A/(Nb-1.))
    A[:,0]=0.5*A[:,0]
    A[:,Nb-1]=0.5*A[:,Nb-1]
    return(A)


# # Fourier coefficient matrix (for multiple  ell-bin)

# In[3]:


# returns A for a cos transform 
def calc_AFT(Nc,Nell):
    Nm=Nc*Nell
    AFT=np.zeros((Nm,Nm))
    A=calc_A(Nc,Nc)
    for ii in range(Nell):
        AFT[ii*Nc:(ii+1)*Nc,ii*Nc:(ii+1)*Nc]=A
    return(AFT)


# # Identify FG modes and bins for 21-cm signal

# In[4]:


def FG_MODES(kper,kpar,dkl,fac):
    
    nc=np.shape(kpar)[0]
    nell=np.shape(kper)[0]
    
    #dkl = np.array([dkl])
    if(len(dkl)==1): # if dkl is passed as a float, make it an array.
        a = dkl[0]
        dkl = np.zeros(nell)+a
    
    NFG=np.zeros(nell,dtype=np.int32)
    # nc=500
    for ii in range(nell):
        for jj in range(nc):
            if(kpar[jj] <=fac*kper[ii]+dkl[ii]):
                NFG[ii]+=1
    #print(NFG)
    return(NFG)

# def FG_MODES(kper,kpar,dkl,fac):
#     nc=np.shape(kpar)[0]
#     nell=np.shape(kper)[0]
#     NFG=np.zeros(nell,dtype=np.int32)
#     # nc=500
#     for ii in range(nell):
#         for jj in range(nc):
#             if(kpar[jj] <=fac*kper[ii]+dkl):
#                 NFG[ii]+=1

#     print(NFG)
#     return(NFG)

def mk_log_bins(kper,kpar,NB,NFG):
# kL used wedge boundary of lowest kperp
 
    # botton left corner of TW
    kL=np.sqrt(kper[0]*kper[0]+kpar[NFG[0]]*kpar[NFG[0]]) 

    nell=np.shape(kper)[0]
    nc=np.shape(kpar)[0]

    # right top corner of TW
    if(NFG[-1]==nc):
        ii=0
        while  ((ii<nell) and (NFG[ii]<nc)) :
            ii+=1
    else:
        ii=nell
        
    kU=np.sqrt(kper[ii-1]*kper[ii-1]+kpar[-1]*kpar[-1])*(1.+.1)
#     kU=np.sqrt(kper[ii-1]*kper[ii-1])*(1.+.1)*np.sqrt(3.2)
    dd=np.log(kU/kL)/NB
    kk=np.arange(NB)
    kk=kL*np.exp((kk+0.5)*dd)
    return(kL,kU,dd,kk)

def calc_bins(kper, kpara, kL, dd, NB, NFG): # returns the mode indices for each bin.
    nl = len(kper)
    nc = len(kpara)
    nm = nl*nc
    A = np.zeros(nm, dtype='int64')
    for ii in range(nl):
        for jj in range(nc):
            m = ii*nc + jj        
            if(jj>=NFG[ii]):
                kk = np.sqrt(kper[ii]**2 + kpara[jj]**2)
                nb = int(np.floor(np.log(kk/kL)/dd))
            else:
                nb = NB # Last bin is FG bin
            A[m] = nb # A matrix has a mode's index, and its corresponding bin number.
            
    bins = np.empty(NB+1, dtype='object')
    for ii in range(NB+1):
        bins[ii] = np.array([], dtype='int16')
    for ii in range(nm):
        nb = A[ii]   
        bins[nb] = np.append(bins[nb], ii) 
        
    return bins


# # Rearrangement  matrix (AR)
# ## takes modes arranged according to $\ell$ and rearranes according to bin

# In[5]:


def calc_AR(bins, nm):        
    AR = np.zeros((nm, nm), dtype='int16')
    row = 0
    for jj in range(len(bins)): 
        for ii in range(len(bins[jj])):
            AR[bins[jj][ii],row] = 1
            row +=1
            
    return AR


# # Orthogonal basis transformation matrix (for a single signal-bin)

# In[6]:


# returns C, orthogonal basis co-efficients for a single bin.
def calc_C(N):
    C = np.zeros((N, N))
    C[0, :] = 1./np.sqrt(N)
    arg = 2.*np.pi*(np.arange(N))/N
    
    if(N%2 == 0): # even N
        for m in range (1, N//2):
            k = 2*(m-1) + 1
            C[k, ] = np.sqrt(2./N)*np.cos(arg*m) 
            C[k+1, ] = np.sqrt(2./N)*np.sin(arg*m)
            
        C[N-1, ] = np.sqrt(1./N)*np.cos(arg*(N//2)) 
        
    else: # odd
        for m in range (1, N//2+1):
            k = 2*(m-1) + 1
            C[k, ] = np.sqrt(2./N)*np.cos(arg*m) 
            C[k+1, ] = np.sqrt(2./N)*np.sin(arg*m)
            
    return(C)


# # Change-of-basis or transition matrix
# ## reaaranges each bin to signal and oprthogonal modes

# In[7]:


def calc_AB(bins, nm):
    AB = np.zeros((nm, nm))
    row = 0
    for ii in range(len(bins)-1): # signal bins
        ne = len(bins[ii])
        AB[row:row+ne, row:row+ne] = calc_C(ne)/np.sqrt(1.*ne)
        row += ne
        
# FG bin
    ne = len(bins[-1])
    AB[row:row+ne, row:row+ne] =np.identity(ne)
    row += ne
        
    return AB
    


# # Signal Isolation Matrix (Asg)
# ## Signal modes from each bin at top followd by all the ortho modes

# In[8]:


def calc_AS(bins, nm):
    iso = np.zeros((nm, nm), dtype='int16')
    nb = len(bins)-1   
    row = nb
    col=0
      
    for jj in range(nb): # Signal only.
        iso[jj, col] = 1
        nn=len(bins[jj])-1
        
        iso[row:row+nn,col+1:col+nn+1]=np.identity(nn)
        row+=nn
        col+=nn+1
    
        
    nfg = len(bins[-1])
    iso[-nfg:, -nfg:] = np.identity(nfg)
        
    return iso


# # Returns matrix which transforms from (Signal, Ortho, FG) to $C_{\ell}(\Delta \nu)$

# In[9]:


def calc_AA(kper,kpara,dkl,fac,NB):

    nl = len(kper)
    nc = len(kpara)
    nm = nl*nc # number of modes 

    NFG = FG_MODES(kper,kpara,dkl,fac) 
    
    kL, kU, dd, kk = mk_log_bins(kper,kpara,NB,NFG)

    bins=calc_bins(kper, kpara, kL, dd, NB, NFG)

    AFT=calc_AFT(nc,nl)

    AR=calc_AR(bins, nm)

    AB=calc_AB(bins, nm)

    AS=calc_AS(bins, nm)

    AA=AFT@AR@np.linalg.inv(AS@AB)
    
#    AA=AFT@AR@np.linalg.inv(AB)@np.linalg.inv(AS)

    return AA,kk


# # Following does not include ortho modes

# # Trnsformation matrix for binned HI signal 

# In[10]:


def calc_MT(kper,kpara,kL,dd,NFG,NB):
#HI contribution to M


    nl = len(kper)
    nc = len(kpara)
    
    MT=np.zeros((nc*nl,NB))
    A=calc_A(nc,nc)
    #counts number of modes in each bin
    nk=np.zeros(NB)


    for ii in range(nl):
        for kk in range(NFG[ii],nc,1):   
            km=np.sqrt(kper[ii]*kper[ii]+kpara[kk]*kpara[kk])
            nb=int(np.floor(np.log(km/kL)/dd))
            nk[nb]+=1.
            for jj in range(nc):
                MT[ii*nc+jj,int(nb)]+=A[jj,kk]
    print(nk)           
    return(MT)


# # Trnsformation matrix for FG modes

# In[11]:


# foreground contribution to M

def calc_MFG(nell,nc,NFG):
    A=calc_A(nc,nc)
    NFGT=np.sum(NFG)
    M=np.zeros((nc*nell,NFGT))
# FG contribution 
    dy=0
    for ii in range(nell):
        M[ii*nc:(ii+1)*nc,int(dy):int(dy+NFG[ii])]=A[:,:NFG[ii]]      
        dy+=NFG[ii]

    return(M)


# # Transformation matrix for (signal,FG) to $C_{\ell}(\Delta \nu)$

# In[15]:


def calc_AM(kper,kpara,dkl,fac,NB):

    nl = len(kper)
    nc = len(kpara)
    nm = nl*nc # number of modes 

    NFG = FG_MODES(kper,kpara,dkl,fac)
    
    kL, kU, dd, kk = mk_log_bins(kper,kpara,NB,NFG)
    
    MT = calc_MT(kper,kpara,kL,dd,NFG,NB)
    
    MFG = calc_MFG(nl,nc,NFG)
    
    AM = np.append(MT,MFG,1)
    
    return AM, kk


# # Given transformation matrix, $C_{\ell}(\Delta \nu)$ and covi, 
# # Returns P  Binned signal and its error covariance COV

# In[13]:


def func_PE(AA,cl,covi,NB,vfac):
    AAt=np.transpose(AA)
    COV=np.linalg.inv(AAt@covi@AA)
    AAf=COV@AAt@covi
    M=AAf@cl    
#     return M[0:NB], COV[0:NB,0:NB]*vfac**2
    return M*vfac, COV*vfac**2



# # Returns $\Delta T$, its rms. and 2 $\sigma$ upper limit in mK units

# def func_dT(kk,p,COV): # SB
#     dp=np.sqrt(np.diagonal(COV))
#     dT=np.sqrt(1.e6*kk**3*p/(2.*np.pi*np.pi))
#     ddT=np.sqrt(1.e6*kk**3*dp/(2.*np.pi*np.pi))
#     ul=np.sqrt(dT**2+2.*ddT**2)
#     return dT, ddT, ul

def func_dT(kk,p,COV): # AE
    dp=np.sqrt(np.diagonal(COV))
    dT2=1.e6*kk**3*p/(2.*np.pi*np.pi)
    dT_Signal = np.copy(dT2)
    dT_Signal [dT_Signal < 0] = 0 # if dT2 is negative, the noise determines the limit. (Hera Collab. 2020)
    ddT2=1.e6*kk**3*dp/(2.*np.pi*np.pi)
    ul2=dT_Signal+2.*ddT2 #dT_Signal+2.*ddT2
    return dT2, ddT2, ul2

# In[ ]:
# Given cylindrical ps, this etimates log-binned pk in traditional way.
def binpk(pk,kper,kpara,dkl,fac,NB):
    NFG = FG_MODES(kper,kpara,dkl,fac)
    kL, kU, dd, kb = mk_log_bins(kper,kpara,NB,NFG)
#     print(kL, kU, dd, kb)
#     print(NFG)
    Nbin = len(kper)
    Nchan = len(kpara)
    
    pks = np.zeros(NB)
    kbin = np.zeros(NB)
    num = np.zeros(NB)
    #Nbin=20
    #Nchan = 300        
    for ii in range (Nbin):
        for jj in range(NFG[ii],Nchan,1): # CB with NFG.
            #if (kpara[jj]>=(fac*kper[ii]+dkl)):
            kk = np.sqrt(kper[ii]**2 + kpara[jj]**2)
            #if ((kk>=kL) and (kk < kU)):
            nb = int(np.floor(np.log(kk/kL)/dd))
            pks[nb] += pk[ii, jj] 
            kbin[nb] += kk 
            num[nb] += 1

    for ii in range (NB):
        if(num[ii]!=0):
            pks[ii] = pks[ii]/num[ii]
            kbin[ii] = kbin[ii]/num[ii]
    # print(num)
    return pks, kb, num




# for Omh1bh1 estimation using PT
def pk_intp(kv, karr, pkarr): # kv, k-array , pk-array
    return np.interp(kv, karr, pkarr)


def calc_MT_PT(kper,kpara,kL,dd,NFG, NB, karr, pkarr): #takes pTM
#HI contribution to M


    nl = len(kper)
    nc = len(kpara)
    
    MT=np.zeros((nc*nl,NB))
    A=calc_A(nc,nc)
    #counts number of modes in each bin
    nk=np.zeros(NB)


    for ii in range(nl):
        for kk in range(NFG[ii],nc,1):   
            km=np.sqrt(kper[ii]*kper[ii]+kpara[kk]*kpara[kk])
            nb=int(np.floor(np.log(km/kL)/dd))
            nk[nb]+=1.
            # interpolate PT at km
            pT = np.interp(km, karr, pkarr)
            
            for jj in range(nc):
                MT[ii*nc+jj,int(nb)]+=A[jj,kk]*pT
    print(nk) 
    # MT = a.flatten()[..., np.newaxis] # pass the model (full cl21)  
    return(MT)

def calc_AM_PT(kper,kpara,dkl,fac,NB, karr, pkarr,vfac):

    nl = len(kper)
    nc = len(kpara)
    nm = nl*nc # number of modes 

    NFG = FG_MODES(kper,kpara,dkl,fac)
    
    kL, kU, dd, kk = mk_log_bins(kper,kpara,NB,NFG)
    
    MT= calc_MT_PT(kper,kpara,kL,dd,NFG,NB, karr, pkarr)
    print(MT.shape)
    
    MFG= calc_MFG(nl,nc,NFG)
    print(MFG.shape)
    
    AM=np.append(MT,MFG,1) #if p21 is passed
    print(AM.shape)
    return AM, kk

def func_PE_PT(AA,cl,covi,NB,vfac):
    AAt=np.transpose(AA)
    COV=np.linalg.inv(AAt@covi@AA)
    AAf=COV@AAt@covi
    M=AAf@cl
    return M*vfac, COV*vfac**2


# pass the model (full cl21)
# pass the model (full cl21)

def calc_MT_A(kper,kpara,kL,dd,NFG,NB,a):
#HI contribution to M
    MT = a.flatten()[..., np.newaxis] # pass the model (full cl21)    
    return(MT)

def calc_AM_A(kper,kpara,dkl,fac,NB,a,vfac):

    nl = len(kper)
    nc = len(kpara)
    nm = nl*nc # number of modes 

    NFG = FG_MODES(kper,kpara,dkl,fac)
    
    kL, kU, dd, kk = mk_log_bins(kper,kpara,NB,NFG)
    
    MT= calc_MT_A(kper,kpara,kL,dd,NFG,NB,a)
    print(MT.shape)
    
    MFG= calc_MFG(nl,nc,NFG)*vfac
    print(MFG.shape)
    
    AM=np.append(MT,MFG,1) # if p21 is passed
    print(AM.shape)
    return AM, kk

def func_PE_A(AA,cl,covi,NB,vfac):
    AAt=np.transpose(AA)
    COV=np.linalg.inv(AAt@covi@AA)
    AAf=COV@AAt@covi
    M=AAf@cl
    return M, COV
