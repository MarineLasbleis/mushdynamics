""" Compaction of a growing layer """

import numpy as np
import matplotlib.pyplot as plt
import yaml
import os
import pandas as pd

import mush

def compaction_column_growth(calcul_velocity, **options):
    """ Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    output_folder = options["output"] + "/"
    if not os.path.isdir(output_folder):
         os.makedirs(output_folder)
    options["Ric_adim"] = radius(options["time_max"], options)
    param_file = output_folder + options["filename"]+'_param.yaml'
    with open(param_file, 'w') as f:
        yaml.dump(options, f) # write parameter file with all input parameters

    psi0 = 1 - options["phi_init"]
    N = 20
    R_init = radius(options["t_init"], options)
    R = np.linspace(0, R_init, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)
    time = options["t_init"]
    dt_print = options["dt_print"]
    time_p = time
    time_max = options["time_max"]
    it = 0
    iter_max = 10000000

    velocity = calcul_velocity(1 - psi, R, options)
    v_m = np.amax(np.abs(velocity))
    dt = min(0.5 * dr / (v_m), 0.5)
    dt = min(dt, dr/growth_rate(time, options))

    stat_file = output_folder + options["filename"]+'_statistics.txt'
    with open(stat_file, 'w') as f:
        f.write("iteration_number time radius radius_size sum_phi r_dot velocity_top max_velocity RMS_velocity thickness_boundary\n")
        f.write('{:d} {:.4e} {:.4e} {:d} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}\n'.format(it, time, R[-1], len(R), average(1-psi, R[1:], options), growth_rate(time, options), velocity[-1], np.max(velocity), average(velocity, R[1:-1], options), thickness_boundary_layer(1-psi, R)))

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        if R[-1]+dr < radius(time, options):
            psi, R = append_radius(psi, R, options)
        velocity = calcul_velocity(1 - psi, R, options)
        psi = mush.update(velocity, psi, dt, R, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.001 * dr / (v_m))
        dt = min(dt, 0.5*dr/growth_rate(time, options))

        with open(stat_file, 'a') as f:
            f.write('{:d} {:.4e} {:.4e} {:d} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}\n'.format(it, time, R[-1], len(R), average(1-psi, R[1:], options), growth_rate(time, options), velocity[-1], np.max(velocity), average(velocity, R[1:-1], options), thickness_boundary_layer(1-psi, R)))

        if time_p > dt_print:
        # if it % 100 == 0:
            data = {"radius": pd.Series(R), 'porosity': pd.Series(1-psi), 'velocity': pd.Series(velocity)}
            data = pd.DataFrame(data)
            mush.output(time, data, fig=False, file=True, output_folder=output_folder, ax=[])
            time_p = time_p - dt_print


def radius(time, options):
    """ Radius of the IC, as function of time. """
    return options["coeff_velocity"]*(time)**options["growth_rate_exponent"]

def growth_rate(time, options):
    """ Growth of the IC, as function of time.

    Correspond to d(radius)/dt
    """
    return options["coeff_velocity"]*time**(1-options["growth_rate_exponent"])

def append_radius(psi, R, options):
    """ Add one element in radius """
    psi = np.append(psi, [options["psiN"]])
    dr = R[1] - R[0]
    R = np.append(R, [R[-1] + dr])
    return psi, R

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

def thickness_boundary_layer(phi, R):
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
            delta = - dr/dphi[it]
            find_it = True
    return delta

def porosity_given_depth(phi, depth, R):
    """ extract the porosity value at the given depth

    if depth is larger than max radius, return 0.
    """
    index = np.argmin(np.abs(R-depth))
    return phi[index]


if __name__ == "__main__":

    r_max = 10.
    t_max = (10/2.)**2
    dt = t_max/20

    options = {'advection': "FLS",
               'delta': 1.,
               'eta': 1.,
               'psi0': 1.,
               'psiN': 0.6,
               'phi_init': 0.4,
               'sign': 1,
               'BC': "dVdz==0",
               'coordinates': "spherical",
               "t_init": 0.01,
               "growth_rate_exponent": 0.5,
               'filename': 'IC_ref',
               'time_max': t_max,
               'dt_print': dt,
               'coeff_velocity': 2., 
               'output': "compaction/"}
    print("Time to be computer: {}, dt for print: {}".format(t_max, dt))
    compaction_column_growth(mush.velocity_Sumita, **options)
    # mush.fig_stat("compaction/IC_ref_statistics.txt")