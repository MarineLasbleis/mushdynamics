import numpy as np
from numba import jit, f8

from scipy.sparse import diags
import scipy.sparse as sparse
import scipy.sparse.linalg

import time

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
@jit(f8[:] (f8[:],f8[:],f8[:],f8[:] ))
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)


    WARNING: in this function, b and d have same size, but a and c are shorter!
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc


def TDMA_2(a,b,c,d):
    n = len(d)
    w= np.zeros(n-1,float)
    g= np.zeros(n, float)
    p = np.zeros(n,float)

    w[0] = c[0]/b[0]
    g[0] = d[0]/b[0]

    for i in range(1,n-1):
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
    p[n-1] = g[n-1]
    for i in range(n-1,0,-1):
        p[i-1] = g[i-1] - w[i-1]*p[i]
    return p




## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolverold(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    print(ac, bc, cc, dc)

    for it in range(1, nf):

        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
                
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc


def inversion_matrice(a, b, c, d):

    N = len(d)
    #print(a, b, c, d)
    # check sizes of a, b, c, d
    if len(a) == len(b): 
        _a = a[1:]
    else:
        _a = a
    if len(c) == len(b): 
        _c = c[:-1]
    else:
        _c = c


    _b, _d = b, d #copy of arrays
    diagonals = [-1,0,1]
    M = diags([[_a], [_b], [_c]], diagonals, [N,N])
    try: 
        output = sparse.linalg.spsolve( M,  _d )
    except Warning as warni:
        print(warni)
        output = np.ones_like(_d)


    return output








if __name__ == "__main__":

    print("======1=======")
    a = np.array([0., 0.])
    b = np.array([1., 1., 1.])
    c = np.array([0., 0.])
    d = np.array([1., 1., 1.])
    print(TDMAsolver(a, b, c, d))
    print(inversion_matrice(a, b, c, d))
    print("=============")

    print("======2=======") 
    a = np.array([2, 2])
    b = np.array([1, 1, 1])
    c = np.array([ 2, 2])
    d = np.array([ 1, 1, 1])
    print(TDMAsolver(a, b, c, d))
    print(TDMAsolverold(a, b, c, d))
    print(TDMA_2(a, b, c, d))
    A = np.array([[1, 2, 0], [2, 1, 2], [0, 2, 1]], dtype=float)
    print(np.linalg.solve(A, d))
    print(inversion_matrice(a, b, c, d))
    print("=============")

    print("======3=======") 
    a = np.array([1, 2])
    b = np.array([1, 1, 1])
    c = np.array([ 2, 0])
    d = np.array([ 1, 2, 3])
    print(TDMAsolver(a, b, c, d))
    print(TDMAsolverold(a, b, c, d))
    print(TDMA_2(a, b, c, d))
    A = np.array([[1, 2, 0], [1, 1, 0], [0, 2, 1]], dtype=float)
    print(np.linalg.solve(A, d))
    print(inversion_matrice(a, b, c, d))
    print("=============")

    print("======4=======") 
    a = np.array([1, 1])
    b = np.array([1, 1, 1])
    c = np.array([ 1, 1])
    d = np.array([ 1, 1, 1])
    A = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]], dtype=float)
    print(np.linalg.solve(A, d))
    print(inversion_matrice(a, b, c, d))
    
    print("=============")

    print("======5=======") 

    a = np.ones(9)
    b = np.ones(10)
    c = np.ones(9)
    d = np.ones(10)
    print(inversion_matrice(a, b, c, d))
    print("=============")


    print("======6=======") 

    A = np.array([[10,2,0,0],[3,10,4,0],[0,1,7,5],[0,0,3,4]],dtype=float)

    a = np.array([3.,1,3])
    b = np.array([10.,10.,7.,4.])
    c = np.array([2.,4.,5.])
    d = np.array([3,4,5,6.])

    print("Test results:")
    print(TDMAsolverold(a, b, c, d))
    print(TDMAsolver(a, b, c, d))
    print(TDMA_2(a, b, c, d))
    print(np.linalg.solve(A, d))
    print(inversion_matrice(a, b, c, d))

    print("Speed results:")
    t0 = time.time()
    for i in range(100000):
        TDMAsolver(a, b, c, d)
    t1 = time.time()
    print("jit_new {}".format(t1-t0))
    t0 = time.time()
    for i in range(100000):
        TDMAsolverold(a, b, c, d)
    t1 = time.time()
    print("old_without {}".format(t1-t0))
    t0 = time.time()
    for i in range(100000):
        np.linalg.solve(A, d)
    t1 = time.time()
    print("control {}".format(t1-t0))

    


    






