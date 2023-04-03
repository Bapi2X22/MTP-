import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# returns A for a cos transform 
def calc_A(Na, Nb):
    A=np.outer(np.arange(Na),np.arange(Nb))
    A=np.cos(np.pi*A/(Nb-1.))
    A[:,0]=0.5*A[:,0]
    A[:,Nb-1]=0.5*A[:,Nb-1]
    return A

# PS using MLE
def func_pk(cl, w, covi, vfac): 
    nell = cl.shape[0]      
    nc = cl.shape[1]      
    A = calc_A(nc, nc)    
    At = np.transpose(A) # A^dagger     
    pk = np.zeros((nell, nc))
    for ii in range(nell):
        vari = np.diag(covi[ii]) # use noise-variance.
        X = np.linalg.inv(At@vari@A)@At@vari
        pk[ii,:] = (X@(w*cl[ii,:]))*vfac
        
    return pk

# PS w/o using MLE
def func_pki(cl, w, vfac): 
    nell = cl.shape[0]      
    nc = cl.shape[1]      
    A = calc_A(nc, nc)    
    X = np.linalg.inv(A) 
    pk = np.zeros((nell, nc))
    for ii in range(nell):
        pk[ii,:] = (X@(w*cl[ii,:]))*vfac
        
    return pk

def func_cl(pk, vfac): # returns cl(dnu) from given p(kper, kpara)
    nell = pk.shape[0]      
    nc = pk.shape[1]      
    A = calc_A(nc, nc)    
    cl = np.zeros((nell, nc))
    for ii in range(nell):
        cl[ii,:] = A@pk[ii,:]/vfac
        
    return cl

def plot_pk(pk, kper, kpara, m, n, fac, dkl):# m,n = No of kper, kpara to exclude from plotting
    y, x = np.meshgrid(kpara, kper)
    plt.pcolormesh(x[:-m,n:], y[:-m,n:], abs(pk[:-m,n:])*1e6, norm = LogNorm(), shading='auto', cmap='coolwarm')
    plt.colorbar()

    plt.plot(x, fac*x, 'k-')
    plt.plot(x, fac*x + dkl, 'b-')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim(y[:-m,n:].min(), y[:-m,n:].max())
    plt.xlim(x[:-m,n:].min(), x[:-m,n:].max())    
    return 

def linbin(array, nb):
    nc = len(array)
    ab = np.zeros(nb)
    bins = np.linspace(0, nc, nb+1, dtype='int64')
    for jj in range(nb):
        ab[jj] = np.average(array[bins[jj]: bins[jj+1]])        
    return ab

def linbin2(array, nb, nx):
    nc = len(array)
    ab = np.zeros(nb)
    nk = nb+1
    bins = np.linspace(nx, nc, nk, dtype='int64')
    if (nx==0):        
        for jj in range(nb):
            ab[jj] = np.average(array[bins[jj]: bins[jj+1]])
        
    else:
        ab[0] = array[0]
        for jj in range(1,nb):
            ab[jj] = np.average(array[bins[jj]: bins[jj+1]])     
            
    return ab

def dct(cl, w, vfac):      
    nc = cl.shape[0]      
    A = calc_A(nc, nc)    
    X = np.linalg.inv(A)    
    pk = (X@(w*cl))*vfac        
    return pk

def idct(pk, vfac): # returns cl(dnu) from given p(kper, kpara)
    nc = pk.shape[0]      
    A = calc_A(nc, nc)    
    cl = A@pk/vfac        
    return cl