""" Compaction of a growing layer """

import numpy as np
import matplotlib.pyplot as plt
import yaml
import os
import pandas as pd
from scipy.optimize import leastsq
import numpy.ma as ma

import mush, data_analysis



def verify_parameters(options):
    """ Verify if the parameters given in options are compatible, then write param file. """

    # calculate the missing Ric_adim, time_max, coeff velocity
    if "Ric_adim" in options and "time_max" in options and "coeff_velocity" in options:
        print("Ric_adim, time_max and coeff_velocity should not be given together as options. Coeff_velocity overwritten by the system.")
        options["coeff_velocity"] = options["Ric_adim"]*options["time_max"]**(-options["growth_rate_exponent"])
    elif not "Ric_adim" in options:
        options["Ric_adim"] = options["coeff_velocity"]*options["time_max"]**options["growth_rate_exponent"]
    elif not "time_max" in options:
        options["time_max"] = (options["Ric_adim"]/options["coeff_velocity"])**(1./options["growth_rate_exponent"])
    elif not "coeff_velocity" in options:
        options["coeff_velocity"] = options["Ric_adim"]/options["time_max"]**options["growth_rate_exponent"]

    options["Ric_adim"] = radius(options["time_max"], options)

    # calculate the R_init / N_init
    if "t_init" in options and "R_init" in options:
        if not radius(options["t_init"], options) == options["R_init"]:
            options["R_init"] = radius(options["t_init"], options)
            print("t_init and R_init should not be both given in options. R_init overwritten to value {}".format(options["R_init"]))
    elif "t_init" in options:
        options["R_init"] = radius(options["t_init"], options)
    elif "R_init" in options:
        options["t_init"] = (R_init/options["coeff_velocity"])**(1./options["growth_rate_exponent"])
    else:
        print("Please provide either t_init or R_init. R_init set to 0.1*R_ic_adim")
        options["R_init"] = 0.1*options["Ric_adim"]
        options["t_init"] = (R_init/options["coeff_velocity"])**(1./options["growth_rate_exponent"])

    try:
        N = options["N_init"]
    except Exception:
        N=20

    output_folder = options["output"] + "/"
    if not os.path.isdir(output_folder):
         os.makedirs(output_folder)
    param_file = output_folder + options["filename"]+'_param.yaml'
    with open(param_file, 'w') as f:
        yaml.dump(options, f) # write parameter file with all input parameters
    return options



def compaction_column_growth(calcul_velocity, **options):
    """ Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    options = verify_parameters(options)

    time = options["t_init"]
    R_init =  options["R_init"]
    N = options["N_init"]
    psi0 = 1 - options["phi_init"]
    
    R = np.linspace(0, R_init, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)
    dt_print = options["dt_print"]
    time_p = time
    time_max = options["time_max"]
    it = 0
    iter_max = 100000000

    velocity = calcul_velocity(1 - psi, R, options)
    v_m = np.amax(np.abs(velocity))
    dt = min(0.5 * dr / (v_m), 0.5)
    dt = min(dt, dr/growth_rate(time, options))

    stat_file = output_folder + options["filename"]+'_statistics.txt'
    with open(stat_file, 'w') as f:
        f.write("iteration_number time radius radius_size sum_phi r_dot velocity_top max_velocity RMS_velocity thickness_boundary\n")
        f.write('{:d} {:.4e} {:.4e} {:d} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}\n'.format(it, time, R[-1], len(R), data_analysis.average(1-psi, R[1:], options), growth_rate(time, options), velocity[-1], np.max(velocity), data_analysis.average(velocity, R[1:-1], options), data_analysis.thickness_boundary_layer(1-psi, R)))

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        # print(it)
        time = time + dt
        time_p = time_p + dt
        if R[-1]+dr < radius(time, options):
            psi, R = append_radius(psi, R, options)
        velocity = calcul_velocity(1 - psi, R, options)
        psi = mush.update(velocity, psi, dt, R, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.001 * dr / (v_m))
        dt = min(dt, 0.5*dr/growth_rate(time, options))

        stat = False
        if it > 1e3:
            if it%100==0:
                stat = True
        else:
            stat = True
        if stat:
            with open(stat_file, 'a') as f:
                f.write('{:d} {:.4e} {:.4e} {:d} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}\n'.format(it, time, R[-1], len(R), data_analysis.average(1-psi, R[1:], options), growth_rate(time, options), velocity[-1], np.max(velocity), data_analysis.average(velocity, R[1:-1], options), data_analysis.thickness_boundary_layer(1-psi, R)))

        if time_p > dt_print:
        # if it % 100 == 0:
            data = {"radius": pd.Series(R), 'porosity': pd.Series(1-psi), 'velocity': pd.Series(velocity)}
            data = pd.DataFrame(data)
            mush.output(time, data, fig=False, file=True, output_folder=output_folder, ax=[])
            time_p = time_p - dt_print


def radius(time, options):
    """ Radius of the IC, as function of time. """
    try:
        if options["supercooling"] == True:
           pass
    except Exception as e:
        pass

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
