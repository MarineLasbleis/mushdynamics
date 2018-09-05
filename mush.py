""" Solver for 2-phase flow equations in 1D

To be used to study compaction in the Earth's inner core.
Boundary conditions are specified to either study
- close boundary at top
- open boundary with fluid expelled from the system.
"""


import numpy as np
import os
import pandas as pd
# import scipy as sc
# import matplotlib.pyplot as plt

import tdma


def schema():
    """ Description of the discretization for numerical scheme"""
    print(" <-1 volume-> \n\
||-----.-----|-----.-----|       ...      .-----|-----.-----|-----.        ...       -----.-----|-----.-----||  \n\
       0     0     1     1                                                                N-1  N-1    N          \n\
 0   1/2dr  dr   3/2dr              |           |           |                       |           |                \n\
             ->          ->               |     ->    |     ->    |                       |     ->    |            \n\
 V0          V1          V2               |    V_i-1  |    V_i    |                       |    V_N-1  |     V_N   \n\
      phi0        phi1                  phi_i-1     phi_i      phi_i+1                phi_N-1      phi_N         \n\
            DP0         DP1                  DP_i-1                                           DP_N-1           \n")


# Nr = 10 #number of points in space
# Nt = 10 #number of points in time


def fluxlimiterscheme(velocity, variable, dr, options):
    """ Advection scheme using the flux limiter scheme

    Coefficients for the scheme are lambda+, lambda-, v+ and v-
    (p and m are used instead of + and m)

    The code return the coefficients for the tridiagonal matrix a, b, c

    The scheme is used to solve the equation
    D/Dt variable = D/Dr (variable * V)
    where D/Dr means partial derivative

    Equation 3.54 in Sramek thesis gives:
    DF/Dx (at the point i) ~ 1/dx *(a_1 variable_i-1 + b_i variable_i +c_i variable_i+1)

    """

    # detect size of the system and initialize variables:
    lambdap, lambdam, vp, vm = np.zeros_like(velocity), np.zeros_like(
        velocity), np.zeros_like(velocity), np.zeros_like(velocity)
    _a, _b, _c, _d = np.zeros_like(variable), np.zeros_like(
        variable), np.zeros_like(variable), np.zeros_like(variable)

    try:
        option = options['advection']
    except KeyError:
        option = 'upwind'

    if option == 'upwind':
        pass  # lambdap and lambdam == zeros.
    elif option == "centered":
        lambdam[:] = np.ones(len(velocity))
        lambdap[:] = np.ones(len(velocity))
    elif option == "FLS":
        Rp = (variable[1:-2] - variable[0:-3]) / \
            (variable[2:-1] - variable[1:-2])
        Rm = (variable[3:] - variable[2:-1]) / \
            (variable[2:-1] - variable[1:-2])
        # minmod
        lambdap[1:-1] = np.fmax(0., np.minimum(1., Rp))
        lambdam[1:-1] = np.fmax(0., np.minimum(1., Rm))
        # at the points [0], [-1]: lambda =0.
    else:
        print("Problem with the choosen scheme for advection. Default is upwind.")

    if len(velocity) == 1:
        if velocity > 0: vp[:] = velocity
        else: vm[:] = - velocity
    else:
        vp[:] = 0.5 * (velocity[:] + np.abs(velocity[:]))
        vm[:] = 0.5 * (velocity[:] - np.abs(velocity[:]))

    _a[1:-1] = -vp[:-1] * (1 - lambdap[:-1] / 2.) - vm[:-1] * lambdam[:-1] / 2.
    _b[1:-1] = vp[1:] * (1 - lambdap[1:] / 2.) + vm[1:] * lambdam[1:] / 2. \
        - vm[:-1] * (1 - lambdam[:-1] / 2.) - vp[:-1] * lambdap[:-1] / 2.
    _c[1:-1] = vm[1:] * (1 - lambdam[1:] / 2.) + vp[1:] * lambdap[1:] / 2.
    _d[1:-1] = (_a[1:-1] * variable[:-2] + _b[1:-1] *
                variable[1:-1] + _c[1:-1] * variable[2:])

    # boundary conditions:
    # velocity fixed at U0 and UN
    # porosity fixed at phi0 and phiN
    # (phi0 and phiN correspond to phi_-1 and phi_N+1,
    # the 2 that are not in the array phi)
    # default is all 0.

    try:
        if options["BC"] == "V==0":
            U0, UN = 0., 0. # velocity[0], velocity[-1]
        elif options["BC"] == "dVdz==0":
            U0, UN = velocity[0], velocity[-1]
    except KeyError:
        U0, UN = 0., 0.
        print("BC not well defined. U0, UN = 0., 0.")

    try:
        psi0, psiN = options["psi0"], options["psiN"]
    except KeyError:
        try:
            psi0, psiN = 1-options["phi0"], 1-options["phiN"]
        except KeyError:
            psi0, psiN = 0.5, 0.5
            print("Be careful, boundary conditions not well specified. \
                    You need to specify psi0, psiN or phi0, phiN")

        # minmod
    # lambdap_up = 0.
    # lambdam_up = 0.
    lambdam_down = 0.
    lambdap_down = 0.

    U0p, U0m = 0.5 * (U0 + np.abs(U0)), 0.5 * (U0 - np.abs(U0))
    UNp, UNm = 0.5 * (UN + np.abs(UN)), 0.5 * (UN - np.abs(UN))

    _a[0] = -U0p * (1 - lambdap_down / 2.) - U0m * lambdam_down / 2.
    _b[0] = vp[0] * (1 - lambdap[0] / 2.) + vm[0] * lambdam[0] / 2. \
        - U0m * (1 - lambdam_down / 2.) - U0p * lambdap_down / 2.
    _c[0] = vm[0] * (1 - lambdam[0] / 2.) + vp[0] * lambdap[0] / 2.
    _d[0] = (_a[0]*psi0+ _b[0] * variable[0] + _c[0] * variable[1]) +_a[0]*psi0

    _a[-1] = -vp[-1]
    _b[-1] = UNp - vm[-1]
    _c[-1] = UNm
    _d[-1] = _a[-1] * variable[-2] + _b[-1] * variable[-1] + _c[-1] * psiN
    _d[-1] = _d[-1] + psiN * _c[-1]

    return _a / (2 * dr), _b / (2 * dr), _c / (2 * dr), _d / (2 * dr)


def CrankNicholson(variable, dr, options):
    """ Calculate the coefficients for the diffusion

    Variable is only here to provide the good shape for the arrays.
    """

    # TODO : find another name for Ra (not the good name, as it is actually
    # 1/Ra)
    try:
        Ra = options["Ra"]
    except KeyError:
        Ra = 0.  # Default value is 0.

    _a, _b, _c, _d = np.zeros_like(variable), np.zeros_like(
        variable), np.zeros_like(variable), np.zeros_like(variable)
    diff = Ra / (dr**2)
    _a[:-1] = - diff / 2.
    _b[:-1] = diff
    _c[:-1] = - diff / 2.
    _d[1:-1] = (variable[:-2] - 2 * variable[1:-1] + variable[2:]) / 2.

    return _a, _b, _c, _d


def velocity_Sramek(variable, radius, options, verbose=False):
    r""" Sramek thesis p46, equation 3.22
    $$ \frac{V}{\delta**2\phi**2} = \frac{d}{dz} [ \frac{(K0+4/3\phi)(1-\phi)}{\phi} \frac{d}{dz}V]-s(1-\phi) $$

    Variable: length N
    Output: velocity is length N-1s
    a, b, c, d are length N-1, N-1, N-1 and N-1
    and a, c are injected as length N-2 for calculating the tridiagonal matrix.
    """
    dr = radius[1] - radius[0]
    try:
        s = options['sign']
    except KeyError:
        s = -1.
        if verbose:
            print("s (sign of density difference) was not defined, please consider defining it for later. Default value is {}".format(s))

    try:
        K0 = options['K0']
    except KeyError:
        K0 = 1.
        if verbose:
            print(
                "K0 was not defined, please consider defining it for later. Default value is {}".format(K0))

    try:
        delta = options["delta"]
    except KeyError:
        delta = 1.
        if verbose:
            print("Delta (compaction length) was not defined, please consider defining it for later. Default value is {}".format(delta))

    _inter = (K0 + 4. / 3. * variable) * (1. - variable) / variable

    if options["coordinates"] == "cartesian":
        _a = _inter[:-1] / dr**2 * variable[:-1] * variable[1:]
        _b = -1. / (delta**2) \
            - _inter[:-1] / dr**2 * variable[:-1] * variable[1:]\
            - _inter[1:] / dr**2 * variable[:-1] * variable[1:]
        _c = _inter[1:] / dr**2 * variable[:-1] * variable[1:]
        _d = s * \
            (1 - np.sqrt(variable[:-1] * variable[1:])) * \
            variable[:-1] * variable[1:]
    elif options["coordinates"] == "spherical":
        _a = _inter[:-1] / dr**2 * variable[:-1] * variable[1:] * \
            radius[:-2]**2 / (radius[0:-2] + dr / 2)**2
        _b = -1. / (delta**2) \
            - _inter[:-1] / dr**2 * variable[:-1] * variable[1:] * radius[1:-1]**2 / (radius[0:-2] + dr / 2)**2\
            - _inter[1:] / dr**2 * variable[:-1] * variable[1:] * \
            radius[1:-1]**2 / (radius[1:-1] + dr / 2)**2
        _c = _inter[1:] / dr**2 * variable[:-1] * variable[1:] * \
            radius[2:]**2 / (radius[1:-1] + dr / 2)**2
        # TODO add a version if no time max (so no growth)
        _d = s / options["Ric_adim"] \
            * (1 - np.sqrt(variable[:-1] * variable[1:])) * \
            variable[:-1] * variable[1:] * radius[1:-1]
        # correction with the term in 4\mu dphi/dr /r
        correction = + 4.*(variable[1:]-variable[:-1])*variable[1:]*variable[:-1]/dr/radius[1:-1]
        _b = _b + correction

    # boundary conditions: V is solved between 0 and N-1,
    # and boundary conditions are forced for V_-1=0 and V_N=0
    # for line 0: V_-1 a_0+ V_0b_0 + V_1 c_0 = d_0,
    # for line N-1: V_N-2 a_N-1+ V_N-1b_N-1 + V_N c_N-1 = d_N-1,
    # so it ends up... Doing nothing should be OK.
    #_a[-1], _b[-1], _c[-1], _d[-1] = 0,1,0,0  #for V=0 at the boundaries
    #_a[0], _b[0], _c[0], _d[0] = 0,1,0,0
    # (would be OK if we wanted to force boundary conditions on 0 and N-1, but we want to calculate the values there)
    # if we wanted V_N = U for example, then d[-1] = d[-1] - U*c[-1]

    # boundary conditions:
    if options["BC"] == "dVdz==0":
        _b[-1] = _b[-1] + _c[-1]
    elif options["BC"] == "V==U":
        _d[-1] = _d[-1] - _c[-1]
    elif options["BC"] == "V==0":
        pass
    new_velocity = tdma.inversion_matrice(_a[1:], _b, _c[:-1], _d)
    return new_velocity


def velocity_Sumita(variable, radius, options, verbose=False):
    """ Solve the velocity of the solid

    based on equation 9 in Sumita et al. 1996
    Cartesian and Spherical coordinates (change in options["coordinates"])
    """
    dr = radius[1] - radius[0]  # assuming no variations of dr

    try:
        K0 = options['K0']
    except KeyError:
        K0 = 1.
        if verbose:
            print(
                "K0 was not defined, please consider defining it for later. Default value is {}".format(K0))

    try:
        eta = options["eta"]
    except KeyError:
        eta = 1.
        if verbose:
            print(
                "eta was not defined, please consider defining it for later. Default value is {}".format(eta))

    try:
        eta0 = options["eta0"]
    except KeyError:
        eta0 = 1.
        if verbose:
            print(
                "eta0 was not defined, please consider defining it for later. Default value is {}".format(eta0))

    try:
        K = options["K"]
    except KeyError:
        K = 1.
        if verbose:
            print(
                "K was not defined, please consider defining it for later. Default value is {}".format(K))

    try:
        grain = options["grain"]
    except KeyError:
        grain = 1
        if verbose:
            print(
                "grain was not defined, please consider defining it for later. Default value is {}".format(grain))

    try:
        sign = options["sign"]
    except KeyError:
        sign = -1.
        if verbose:
            print(
                "sign was not defined, please consider defining it for later. Default value is {}".format(sign))
    if options["coordinates"] == "cartesian":
        _a = - ((1. / (dr**2.)) * ((1. - variable[0:-1])**2.)
                * (4. / (3. * variable[0:-1])) * (eta / eta0))
        _b = ((1. - np.sqrt(variable[1:] * variable[0:-1]))**2 / (variable[0:-1] * variable[1:])**(3. / 2.)) * ((K * K0) / grain**2.) \
                + (1. / dr**2.) * (((1. - variable[0:-1])**2.) * (4. / (3. * variable[0:-1])) * (
                    eta / eta0) + ((1. - variable[1:])**2.) * (4. / (3. * variable[1:])) * (eta / eta0))
        _c = - ((1. / (dr**2.)) *
                ((1. - variable[1:])**2.) * (4. / (3. * variable[1:])) * (eta / eta0))
        _d = - sign * ((1. - np.sqrt(variable[1:] * variable[0:-1])))
    elif options["coordinates"] == "spherical":
        _a = -((4 * (1 - variable[0:-1])**2) / (3 * variable[0:-1])) * (eta / eta0) \
            * (1 / ((radius[0:-2]) + (dr / 2))**2) * ((radius[0:-2])**2 / dr**2)
        _b = ((K * K0) / grain**2) * ((1. - np.sqrt(variable[1:] * variable[0:-1]))**2 / (variable[0:-1] * variable[1:])**(3. / 2.)) \
            + (((4 * (1 - variable[1:])**2) / (3 * variable[1:])) * (eta / eta0)
                * (1 / ((radius[1:-1]) + (dr / 2))**2) * (((radius[1:-1])**2 / dr**2))) \
            + (((4 * (1 - variable[0:-1])**2) / (3 * variable[0:-1])) * (eta / eta0)
                * (1 / ((radius[0:-2]) + (dr / 2))**2) * (((radius[1:-1])**2 / dr**2)))
        _c = -((4 * (1 - variable[1:])**2) / (3 * variable[1:])) * (eta / eta0) \
            * (1 / ((radius[1:-1]) + (dr / 2))**2) * (((radius[2:])**2 / dr**2))
         # TODO add a version if no time max (so no growth)
        _d = - sign / options["Ric_adim"] * \
            ((1. - np.sqrt(variable[1:] * variable[0:-1])) * (radius[1:-1]))

    # boundary conditions:
    if options["BC"] == "dVdz==0":
        _b[-1] = _b[-1] + _c[-1]
    elif options["BC"] == "V==U":
        _d[-1] = _d[-1] - _c[-1]
    elif options["BC"] == "V==0":
        pass

    # if phi is too close to 1 for the system to converge to a velocity
    too_large = (variable[:-1] > 1. - 1e-6)
    _a = np.where(too_large, 0., _a)
    _b = np.where(too_large, 1., _b)
    _c = np.where(too_large, 0., _c)
    _d = np.where(too_large, 0., _d)

    new_velocity = tdma.inversion_matrice(_a[1:], _b, _c[:-1], _d)
    return new_velocity


def source_spherical_advection(psi, velocity, radius, options):
    """ Source term for advection in spherical coordiantes

    Equation 3.60 in Sramek PhD thesis
    psi: (1-porosity) np.array length: N
    velocity: np.array length: N-1
    radius: np.array length: N+1
    output: d (np.array ) length N
    """
    dr = radius[1]-radius[0]
    d = np.zeros_like(psi)
    d[1:-1] = psi[1:-1] * (velocity[1:]+velocity[:-1])/(radius[2:-1]+dr/2)
    # V==0 at z=0
    d[0] = psi[0]*velocity[0]/(dr/2)
    if options['BC'] == 'V==0':
        d[-1] = psi[-1]*velocity[-1]/(radius[-1]-dr/2)
    elif options['BC'] == 'dVdz==0':
        d[-1] = 2.*psi[-1]*velocity[-1]/(radius[-1]-dr/2)
    return d


def update(V, psi, dt, radius, options):
    """ Solve the porosity evolution equation

    advection of the variable psi (1-porosity) by the velocity V
    options is a dictionnary, with a bunch of parameters for everything.
    """
    dr = radius[1]-radius[0]
        # a_adv, b_adv, c_adv, d_adv = fluxlimiterscheme(V, phi, dr, options)
        # a_diff, b_diff, c_diff, d_diff = CrankNicholson(phi, dr, options)
        # _a, _b, _c, _d = a_adv+a_diff, b_adv+b_diff, c_adv+c_diff, d_adv+d_diff
    _a, _b, _c, _d = fluxlimiterscheme(V, psi, dr, options)

    if options["coordinates"] == "spherical":
        d_spherical = source_spherical_advection(psi, V, radius, options)
        _d = _d + d_spherical

    _a = _a * dt
    _b = 1. + _b * dt
    _c = _c * dt
    _d = psi - _d * dt
    psi2 = np.zeros_like(psi)
    _psi = tdma.inversion_matrice(_a[1:], _b, _c[:-1], _d)
    psi2[:] = _psi
    return psi2


def output(time, psi, velocity, R, fig=False, file=False, output_folder="", ax=[]):
    """ function to write output

    fig: True for output a figure. Require to provide a figure ax variable
    file: True for output a file
    output_folder: name of the folder for the output
    """
    # verifiy if folder exist
    if not os.path.isdir(output_folder):
         os.makedirs(output_folder)

    # write the outputs
    if fig:
        dr = R[1]-R[0]
        ax[0].plot(1 - psi, R[:-1] + dr / 2.)
        ax[1].plot(velocity, R[1:-1])
    if file:
        file = output_folder+"/output_{:5.2f}.csv".format(time)
        _data = {"radius": pd.Series(R), 'porosity': pd.Series(1-psi), 'velocity': pd.Series(velocity)}
        data = pd.DataFrame(_data)
        data.to_csv(file)


if __name__ == '__main__':

    pass
    # plt.savefig("sumita_phi03.pdf")
