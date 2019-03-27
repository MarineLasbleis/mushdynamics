import numpy as np
from scipy.optimize import leastsq
import numpy.ma as ma


def average(variable, R, options):
    """ Average of a variable, weigthed by volume of elements. """
    dr = R[1]-R[0] # constant steps in radius/heigth
    if options["coordinates"] == "cartesian":
        dV = dr*np.ones_like(R)
    elif options["coordinates"] == "spherical":
        dV = 4*np.pi*dr*R**2
    return np.sum(variable*dV)/np.sum(dV) #in cartesian

def flux_top(phi, velocity):
    """ Flux of solid from the top boundary """
    return (1-phi[-1])*velocity[-1]

def thickness_boundary_layer_old(phi, R):
    """ Thickness of the mushy zone at the top """
    # find the first inflexion point (starting from top)
    # of the porosity field
    dr = R[1]-R[0]
    dphi = phi[2:] - phi[:-2]
    d2phi = phi[2:] + phi[:-2] -2*phi[1:-1]
    find_it = False
    it = 0
    init_sign = np.sign(d2phi[-1])
    while find_it == False and it <= len(phi)-2:
        if it == len(phi)-2:
            delta = 0.
            it += 1
        elif np.sign(d2phi[-it]) == init_sign:
            it = it+1
        else:
            delta = - dr/dphi[-it]    # -it???
            find_it = True
    return delta

def thickness_boundary_layer(phi, R):
    """ Thickness of the mushy zone at the top """
    def logistic4(x, A, B, C, D):
        """4PL logistic equation."""
        return ((A-D)/(1.0+((x/C)**B))) + D
    def residuals(p, y, x):
        """Deviations of data from fitted 4PL curve"""
        A,B,C,D = p
        err = y-logistic4(x, A, B, C, D)
        return err
    # find 1st minimum of function
    min_phi = np.argmin(phi)
    # fit with sigmoid function
    dr = R[1] -R[0]
    y_values = phi[min_phi:]
    radius = R[min_phi:-1] +dr/2
    # initial guess
    p0 = [phi[min_phi], np.abs((phi[-1]-phi[min_phi])/(radius[-1]-radius[0])), (radius[-1]+radius[0])/2 , phi[-1]]
    plsq = leastsq(residuals, p0, args=(y_values, radius))
    # get the thickness of the sigmoid
    A, B, C, D = plsq[0]
    phi_c = logistic4(C, A, B, C, D)
    return radius[-1]-C+phi_c*4*C/B/(D-A)

def thickness(phi, R):
    phi_inverse = phi[::-1]
    dr = R[1]-R[0]
    min_phi = np.argmin(phi_inverse)
    phi_inverse = phi_inverse[0:min_phi+1]
    max_value = phi[-1]
    mean_phi = (phi_inverse[-1] + phi_inverse[0])/2
    if min_phi <3:
        delta = 0.
    else:
        i_middle = np.argmin(np.abs(phi_inverse-mean_phi))
        dphi = phi_inverse[1:] - phi_inverse[:-1]
        slope = dphi[i_middle]/dr
        delta = -mean_phi/slope
    return delta

def porosity_compacted_region(phi, R, delta, options):
    """ Porosity in the compacted region

    delta: thickness of mushy zone
    """
    if 1.5*delta > R[-1]:
        print("mushy zone larger than core. No compacted region.")
        return 0.
    R_max = R[-1]-1.5*delta
    mask = R<R_max
    return average(ma.masked_array(phi_comp, mask=mask), ma.masked_array(R_comp, mask=mask), options)

def porosity_given_depth(phi, depth, R):
    """ extract the porosity value at the given depth

    if depth is larger than max radius, return 0.
    """
    index = np.argmin(np.abs(R-depth))
    return phi[index]