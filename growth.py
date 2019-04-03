""" Compaction of a growing layer """

import numpy as np
import matplotlib.pyplot as plt
import yaml
import os
import pandas as pd
from scipy.optimize import leastsq
import numpy.ma as ma

import mush, data_analysis


class Compaction():

    def __init__(self, calcul_velocity, **options):
        self.options = options
        self.verify_parameters()
        print(self.options)
        self.print_param()
        self.iter_max = 100000000
        self.calcul_velocity = calcul_velocity
        self.output_folder = options["output"] + "/"
        self.initialisation()

    def verify_parameters(self):
        self.options = verify_parameters(self.options)

    def run(self):
        self.initialisation()
        while self.time < self.time_max and self.it < self.iter_max:
            self.it += 1
            self.time += self.dt
            self.time_p += self.dt
            self.one_step()
            self.write_stat()
            self.write_profile()

    def initialisation(self):
        self.it = 0
        self.time = self.options["t_init"]
        self.R_init =  self.options["R_init"]
        self.N = self.options["N_init"]
        self.psi0 = 1 - self.options["phi_init"]
        self.R = np.linspace(0, self.R_init, self.N + 1)
        self.dr = self.R[1] - self.R[0]
        self.psi = self.psi0 * np.ones(self.N)
        self.dt_print = self.options["dt_print"]
        self.time_p = self.time
        self.time_max = self.options["time_max"]
        # 1st run
        self.velocity = self.calcul_velocity(1 - self.psi, self.R, self.options)
        v_m = np.amax(np.abs(self.velocity))
        dt = min(0.5 * self.dr / (v_m), 0.5)
        self.dt = min(dt, self.dr/self.growth_rate(self.time))
        # init stat file
        self.stat_file = self.output_folder + self.options["filename"]+'_statistics.txt'
        with open(self.stat_file, 'w') as f:
            f.write("iteration_number time radius radius_size sum_phi r_dot velocity_top max_velocity RMS_velocity thickness_boundary\n")
            f.write('{:d} {:.4e} {:.4e} {:d} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}\n'.format(self.it, self.time,\
                                                self.R[-1], len(self.R), data_analysis.average(1-self.psi, self.R[1:], \
                                                self.options), self.growth_rate(self.time), self.velocity[-1], np.max(self.velocity), \
                                                data_analysis.average(self.velocity, self.R[1:-1], self.options), \
                                                data_analysis.thickness(1-self.psi, self.R)))

    def one_step(self):
        if self.R[-1]+self.dr < self.radius(self.time):
            self.psi, self.R = append_radius(self.psi, self.R, self.options)
        self.velocity = self.calcul_velocity(1 - self.psi, self.R, self.options)
        self.psi = mush.update(self.velocity, self.psi, self.dt, self.R, self.options)
        v_m = np.amax(np.abs(self.velocity))
        dt = min(0.5, 0.05 * self.dr / (v_m))
        self.dt = min(dt, 0.05*self.dr/self.growth_rate(self.time))

    def write_stat(self):
        stat = False
        if self.it > 1e3:
            if self.it%100==0:
                stat = True
        else:
            stat = True
        if stat:
            with open(self.stat_file, 'a') as f:
                f.write('{:d} {:.4e} {:.4e} {:d} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}\n'.format(self.it, self.time,\
                                                self.R[-1], len(self.R), data_analysis.average(1-self.psi, self.R[1:], \
                                                self.options), self.growth_rate(self.time), self.velocity[-1], np.max(self.velocity), \
                                                data_analysis.average(self.velocity, self.R[1:-1], self.options), \
                                                data_analysis.thickness(1-self.psi, self.R)))

    def write_profile(self):
        if self.time_p > self.dt_print:
            # if it % 100 == 0:
            data = {"radius": pd.Series(self.R), 'porosity': pd.Series(1-self.psi), 'velocity': pd.Series(self.velocity)}
            data = pd.DataFrame(data)
            mush.output(self.time, data, fig=False, file=True, output_folder=self.output_folder, ax=[])
            self.time_p += -self.dt_print

    def radius(self, time):
        return radius(time, self.options)

    def growth_rate(self, time):
        return growth_rate(time, self.options)

    def print_param(self):
        output_folder = self.options["output"] + "/"
        if not os.path.isdir(output_folder):
            os.makedirs(output_folder)
        param_file = output_folder + self.options["filename"]+'_param.yaml'
        with open(param_file, 'w') as f:
            yaml.dump(self.options, f, default_flow_style=False) # write parameter file with all input parameters


class Compaction_Supercooling(Compaction):

    def verify_parameters(self):
        self.options = verify_parameters(self.options)
        #self.options["r0_supercooling"] = self.options["Ric_adim"]/self.options["time_max"]**self.options["growth_rate_exponent"]\
        #                                    *(self.options["t0_supercooling"]+self.options["Dt_supercooling"])**self.options["growth_rate_exponent"]
        self.options["Dt_supercooling"] = (self.options["r0_supercooling"]/self.options["Ric_adim"])**(1./self.options["growth_rate_exponent"])\
                                            *self.options["time_max"] - self.options["t0_supercooling"]
        print("Dt supercooling {}".format(self.options["Dt_supercooling"] ))
        self.options["tic"] = self.options["time_max"]
        self.options["time_max"] += -self.options["Dt_supercooling"]

    def radius(self, time):
        if time<self.options["t0_supercooling"]:
            radius = self.options["r0_supercooling"]/self.options["t0_supercooling"]*time
        else:
            radius = self.options["Ric_adim"]/self.options["tic"]**self.options["growth_rate_exponent"]\
                                            *(time+self.options["Dt_supercooling"])**self.options["growth_rate_exponent"]
        return radius

    def growth_rate(self, time):
        if time<self.options["t0_supercooling"]:
            growth_rate = self.options["r0_supercooling"]/self.options["t0_supercooling"]
        else:
            growth_rate = self.options["Ric_adim"]/self.options["tic"]**self.options["growth_rate_exponent"]\
                                            *(time+self.options["Dt_supercooling"])**(self.options["growth_rate_exponent"]-1)* self.options["growth_rate_exponent"]
            #growth_rate = self.options["Ric_adim"]/self.options["tic"]**self.options["growth_rate_exponent"]\
            #                                *self.options["growth_rate_exponent"]*(time)**(self.options["growth_rate_exponent"]-1)
        return growth_rate


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
        delta = data_analysis.thickness_boundary_layer(1-psi, R)
        f.write("iteration_number time radius radius_size sum_phi r_dot velocity_top max_velocity RMS_velocity thickness_boundary\n")
        f.write('{:d} {:.4e} {:.4e} {:d} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}\n'.format(it, time, R[-1], len(R), data_analysis.average(1-psi, R[1:], options), growth_rate(time, options), velocity[-1], np.max(velocity), data_analysis.average(velocity, R[1:-1], options), delta , data_analysis.porosity_compacted_region(1-psi, R, delta, options)))

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
        dt = min(0.5, 0.05 * dr / (v_m))
        dt = min(dt, 0.05*dr/growth_rate(time, options))

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
    return options["coeff_velocity"]*(time)**options["growth_rate_exponent"]

def growth_rate(time, options):
    """ Growth of the IC, as function of time.

    Correspond to d(radius)/dt
    """
    #return options["Ric_adim"]/options["tic"]**options["growth_rate_exponent"]\
    #                                        *options["growth_rate_exponent"]*(time)**(options["growth_rate_exponent"]-1)
    if options["coeff_velocity"] == 0.:
        return 0.
    else:
        return options["coeff_velocity"]*time**(options["growth_rate_exponent"]-1)*options["growth_rate_exponent"]

def append_radius(psi, R, options):
    """ Add one element in radius """
    psi = np.append(psi, [options["psiN"]])
    dr = R[1] - R[0]
    R = np.append(R, [R[-1] + dr])
    return psi, R

def verify_parameters(options):
    """ Verify if the parameters given in options are compatible, then write param file. """
    # calculate the missing Ric_adim, time_max, coeff velocity

    # ce serait plus elegant de faire autrement...
    if options["coeff_velocity"] ==0.:
        return options

    if "Ric_adim" in options and "time_max" in options and "coeff_velocity" in options:
        print("Ric_adim, time_max and coeff_velocity should not be given together as options. Coeff_velocity overwritten by the system.")
        options["coeff_velocity"] = options["Ric_adim"]*options["time_max"]**(-options["growth_rate_exponent"])
    elif not "Ric_adim" in options:
        options["Ric_adim"] = options["coeff_velocity"]*options["time_max"]**options["growth_rate_exponent"]
    elif not "time_max" in options:
        options["time_max"] = (options["Ric_adim"]/options["coeff_velocity"])**(1./options["growth_rate_exponent"])
    elif not "coeff_velocity" in options:
        options["coeff_velocity"] = options["Ric_adim"]/options["time_max"]**options["growth_rate_exponent"]

    # calculate the R_init / N_init
    if "t_init" in options and "R_init" in options:
        if not radius(options["t_init"], options) == options["R_init"]:
            options["R_init"] = radius(options["t_init"], options)
            print("t_init and R_init should not be both given in options. R_init overwritten to value {}".format(options["R_init"]))
    elif "t_init" in options:
        options["R_init"] = radius(options["t_init"], options)
    elif "R_init" in options:
        options["t_init"] = (options["R_init"]/options["coeff_velocity"])**(1./options["growth_rate_exponent"])
    else:
        print("Please provide either t_init or R_init. R_init set to 0.1*R_ic_adim")
        options["R_init"] = 0.1*options["Ric_adim"]
        options["t_init"] = (options["R_init"]/options["coeff_velocity"])**(1./options["growth_rate_exponent"])

    try:
        N = options["N_init"]
    except Exception:
        N=20
    return options

def print_param(options):
    output_folder = options["output"] + "/"
    if not os.path.isdir(output_folder):
         os.makedirs(output_folder)
    param_file = output_folder + options["filename"]+'_param.yaml'
    with open(param_file, 'w') as f:
        yaml.dump(options, f) # write parameter file with all input parameters

def plot_growth():
    # we need to provide some options, but the only ones we will use are the growth history's ones!
    def options(r, t_max, exp):
        coeff = r/(t_max**exp)
        print(coeff)
        #t_max **exp = r/coeff
        #t_max = (r/coeff)**(1/exp)
        opt = {     'output': " ",
                    'filename': "",
                    'coordinates': "spherical",
                    "growth_rate_exponent": exp,
                    'time_max': t_max,
                    'coeff_velocity': coeff, 
                    'N_init':4,
                    "R_init": 1., 
                    "phi_init": 0., 
                    "phi_0": 0., 
                    "dt_print": 0.,
                    "BC": ""
                    }
        return opt


    Model_t1 = Compaction(mush.velocity_Sramek, **options(10., 1., 1.))
    Model_t05 = Compaction(mush.velocity_Sramek, **options(10., 1., 0.5))
    Model_t03 = Compaction(mush.velocity_Sramek, **options(10., 1., 1./3.))

    opt = options(10., 1., 1./3)
    opt["t0_supercooling"] = 1e-3
    opt["r0_supercooling"] = 7.
    Model_supercooling = Compaction_Supercooling(mush.velocity_Sramek, **opt)

    models = [Model_t1, Model_t05, Model_t03] # , Model_supercooling]
    labels = ["Linear growth", "$r \sim t^{1/2}$", "$r \sim t^{1/3}$"]
    lines = ["-", ":", "--"]

    fig, ax = plt.subplots(2, 1, figsize=[6, 5], sharex=True)
    for j, mod in enumerate(models): 
        tmax = mod.time_max
        time = np.linspace(0, tmax, 100)
        radius = np.zeros_like(time)
        growth = np.zeros_like(time)
        for i, t in enumerate(time):
            radius[i] = mod.radius(t)
            growth[i] = mod.growth_rate(t)
        ax[0].plot(time, radius, lines[j], label=labels[j])
        ax[1].plot(time, growth, lines[j], label=labels[j])

    mod = Model_supercooling
    tmax = mod.time_max
    time = np.linspace(0, tmax, 100)
    radius = np.zeros_like(time)
    growth = np.zeros_like(time)
    for i, t in enumerate(time):
            radius[i] = mod.radius(t)
            growth[i] = mod.growth_rate(t)
    ax[0].plot(time+mod.options["Dt_supercooling"], radius, "-.", label="Delayed nucleation")
    ax[1].plot(time+mod.options["Dt_supercooling"], growth, "-.", label="Delayed nucleation")
    #self.options["Dt_supercooling"] 

    ax[1].set_ylim([0,15])
    ax[0].set_ylim([0,10])
    ax[0].set_xlim([0,1])
    ax[0].legend()
    ax[1].legend()

    ax[1].set_xlabel("Time")
    ax[0].set_ylabel("Radius")
    ax[1].set_ylabel("Growth rate")
    plt.savefig("growth.pdf")
    plt.show()

if __name__ == "__main__":

    plot_growth()
