import numpy as np

# Returns widnow of same size as X - center at X[0]
def window(typ,x):
    if(typ=='BN'):
        A=[3.6358e-1, 4.891775e-1, 1.3659e-1, 1.06411e-2, 0., 0., 0.]
    if(typ=='BH4'):
        A=[3.5875e-1, 4.8829e-1, 1.4128e-1, 1.1680e-2, 0., 0., 0.]
    if(typ=='MS5'):
        A=[3.2321e-1, 4.7149e-1, 1.7553e-1, 2.8496e-2, 1.2613e-3, 0., 0.]
    if(typ=='MS6'):
        A=[2.9355e-1, 4.5193e-1, 2.0141e-1, 4.7926e-2, 5.0261e-3, 1.3755e-4, 0.]
    if(typ=='MS7'):
        A=[2.7122e-1, 4.3344e-1, 2.1800e-1, 6.5785e-2, 1.07618e-2, 7.7001e-4, 1.3680e-5]
    NC=len(x)
    w=np.zeros(NC)    
    
    for ii in range (NC):
        for jj in range(len(A)):
            ph=np.cos(2.*(ii+NC-1)*jj*np.pi/(2.*NC-2.))
            w[ii]+=(-1.)**jj*A[jj]*ph        

    if(typ=='none'):
        w[:]=1.
    
    return(w)