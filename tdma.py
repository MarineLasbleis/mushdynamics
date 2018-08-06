import numpy as np

from scipy.sparse import diags
import scipy.sparse as sparse
import scipy.sparse.linalg



def inversion_matrice(a, b, c, d):
    """ Solve the system Ax=d, with A a tridiagonal matrix construct with a,b,c. """
    N = len(d)
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

    






