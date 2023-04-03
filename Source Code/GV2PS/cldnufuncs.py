import numpy as np
import multiprocessing as mp # multiprocessing

def collapse(gv, NC): 
    na = np.shape(gv)[2]//NC
    Nu = np.shape(gv)[0]
    cgv = np.zeros((Nu,Nu,na), dtype='complex128')
    bb=0
    for ii in range(na):
        aa = bb
        bb = aa+NC
        cgv[:,:,ii] = np.mean(gv[:,:,aa:bb], axis=-1)  
    return cgv

def read_header(filename, headersize, nc):    
    aa = np.fromfile(filename, dtype='float64', count=headersize, offset=0)      
    FWHM, dUg, f, nua = aa[0], aa[1], aa[2], aa[3:nc+3]
    return FWHM, dUg, f, nua

def read_vcg(filename, nc, NC, xn):
    headersize = (nc+3)*np.float64().itemsize
    aa = np.fromfile(filename, dtype='float64', count=headersize, offset=0)      
    FWHM, dUg, f, nua = aa[0], aa[1], aa[2], aa[3:nc+3] # header
    hinfo = (FWHM, dUg, f, nua)
    aa = np.fromfile(filename, dtype='complex128', count=-1, offset=headersize)
    Nu = int(np.sqrt(len(aa)/(3*nc)))
    NN = Nu*Nu*nc
    GR = collapse(np.reshape(aa[:NN], (Nu, Nu, nc))[xn:Nu-xn, xn:Nu-xn, :], NC)
    GL = collapse(np.reshape(aa[NN:2*NN], (Nu, Nu, nc))[xn:Nu-xn, xn:Nu-xn, :], NC)
    AA = collapse(np.reshape(aa[2*NN:3*NN], (Nu, Nu, nc))[xn:Nu-xn, xn:Nu-xn, :], NC)
    return hinfo, GR, GL, np.real(AA), np.imag(AA)

def gvcorr(vig, aig): # old TGE, self-term subtrated only from \delta\nu=0
    na = len(vig)
    kg = np.zeros(na)
    cov = np.real(np.outer(vig, vig.conjugate())) # gives maps nc*nc
    for dnu in range(na):
        kg[dnu] = np.trace(cov, offset=dnu)-aig[dnu] # noise bias subtraction from all dnu
        # kg[dnu] = np.trace(cov, offset=dnu) # no noise bias subtraction
    # kg[0] -= aig[0] # noise bias subtraction only from dnu=0   
    return kg

def func_cldnu_multi(fvcg, nc, NC, xn, Ncores):       
    hinfo, vR, vL, aR, aL = read_vcg(fvcg, nc, NC, xn)
    vcg = vR+vL 
    ag = aR+aL
    
    Nu = np.shape(vcg)[0]
    na = np.shape(vcg)[2] 
    
    pool = mp.Pool(Ncores)
    eg1 = np.array(pool.starmap(gvcorr, [(vcg[ii, jj], ag[ii, jj]) for ii in range (Nu) for jj in range (Nu)]))
    pool.close()    
    # eg1 = np.array([gvcorr(vcg[ii, jj], ag[ii, jj]) for ii in range(Nu) for jj in range(Nu)]) # serial 
    
    eg1 = eg1.reshape(Nu,Nu,na)         
    return eg1

def gvcorr_CCP(viR, viL): # Cross-correlation
    na = len(viR)
    kg = np.zeros(na)
    
    cov = np.real(np.outer(viR, viL.conjugate()) + np.outer(viL, viR.conjugate())) # gives maps nc*nc
    for dnu in range(na):
        kg[dnu] = np.trace(cov, offset=dnu) # no noise bias subtraction
        # kg[dnu] = np.trace(cov, offset=dnu)-aig[dnu] # noise bias subtraction from all dnu
    # kg[0] -= aig[0] # noise bias subtraction only from dnu=0   
    return kg

def func_cldnu_multi_CCP(fvcg, nc, NC, xn, Ncores): # Cross-correlation      
    hinfo, vR, vL, aR, aL = read_vcg(fvcg, nc, NC, xn)
    #vcg = vR+vL 
    #ag = aR+aL
    
    Nu = np.shape(vR)[0]
    na = np.shape(vR)[2] 
    
    pool = mp.Pool(Ncores)
    eg1 = np.array(pool.starmap(gvcorr_CCP, [(vR[ii, jj], vL[ii, jj]) for ii in range (Nu) for jj in range (Nu)]))
    pool.close()    
    # eg1 = np.array([gvcorr(vcg[ii, jj], ag[ii, jj]) for ii in range(Nu) for jj in range(Nu)]) # serial 
    
    eg1 = eg1.reshape(Nu,Nu,na)         
    return eg1

def mk_lin_bin(mg, dUg, U1, U2, NB): # linear binning
    dd = (U2-U1)/NB    
    Ng = np.shape(mg)[0]    
    
    el = np.zeros(NB)
    wl = np.zeros(NB)
    
    Nm = (Ng-1)//2    
    for ii in range(Ng):
        for jj in range(Nm,Ng):
            UU=np.sqrt((ii-Nm)*(ii-Nm)+(jj-Nm)*(jj-Nm))*dUg
            if ((UU>=U1) and (UU < U2)):
                nb=int(np.floor((UU-U1)/dd)) # nb=int(np.floor(np.log(UU/U1)/dd)) # log-binning 
                wg = mg[ii,jj,0]
                el[nb] += wg*UU
                wl[nb] += wg
                
    ub = el/wl            
    return ub, dd

def func_bincl(eg, dUg, U1, U2, NB): # linear binning
    dd = (U2-U1)/NB    
    Ng = np.shape(eg)[0]    
    na = np.shape(eg)[2]    
    el = np.zeros((NB,na))    
    Nm = (Ng-1)//2
    
    for ii in range(Ng):
        for jj in range(Nm,Ng):
            UU=np.sqrt((ii-Nm)*(ii-Nm)+(jj-Nm)*(jj-Nm))*dUg
            if ((UU>=U1) and (UU < U2)):
                nb=int(np.floor((UU-U1)/dd)) # nb=int(np.floor(np.log(UU/U1)/dd)) # log-binning 
                el[nb] += eg[ii,jj]                
    return el    
