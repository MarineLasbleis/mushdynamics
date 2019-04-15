""" Series of runs using yaml file as input """


import sys
import argparse
import yaml
import numpy as np
import numpy.random as random

import mush
import growth





def run(options):
    """ Run a model with the provided options"""
    print(options["output"])
    print(options)
    #Model = growth.Compaction_Supercooling(mush.velocity_Sramek, **options)
    model = growth.Compaction(mush.velocity_Sramek, **options)
    model.run()

def run_supercooling(options):
    print(options["output"])
    print(options)
    Model = growth.Compaction_Supercooling(mush.velocity_Sramek, **options)
    #Model = growth.Compaction(mush.velocity_Sramek, **options)
    Model.run()

def read_options(filename):
    with open(filename, 'w') as infile:
        data = yaml.read(infile)
    return data

def write_options(options, filename):
    print(options)
    with open(filename, 'w') as outfile:
        yaml.dump(options, outfile, default_flow_style=False)

def default_yaml():
    r_max = 10.
    t_max = (10/2.)**2
    dt = t_max/20
    options = {'advection': "FLS",
                'delta': 1.,
                'eta': 1.,
                'psi0': 1.,
                'psiN': 0.6,
                'phi_init': 0.4,
                'K0': 1.,
                'sign': 1.,
                'BC': "dVdz==0",
                'coordinates': "spherical",
                "growth_rate_exponent": 0.5,
                'filename': 'IC_Sramek',
                'time_max': t_max,
                'dt_print': dt,
                'coeff_velocity': 2.,
                'output': "compaction/",
                "R_init": 0.001,
                "N_init": 5,
                "t_init": 0.} # min(5, int(4e3/r_max))}
    return options

def modify_param(options, new_options):
    return {**options, **new_options}

def param_growth(r, exp, coeff, n=2, N_fig=20, basefolder="", R_init=1e-3, N_max=5000):
    t_max = (r/coeff)**(1/exp)
    dt = t_max/N_fig
    folder_name = basefolder+"/exp_{:.2e}_coeff_{:.2e}_radius_{:.2e}".format(exp, coeff, r)
    options = {'advection': "FLS",
                'n': n,
                'delta': 1.,
                'eta': 1.,
                'psi0': 1.,
                'psiN': 0.6,
                'phi_init': 0.4,
                'K0': 1.,
                'sign': 1.,
                'BC': "dVdz==0",
                'coordinates': "spherical",
                "growth_rate_exponent": exp,
                'filename': 'IC',
                'time_max': t_max,
                'dt_print': dt,
                'coeff_velocity': coeff,
                'output': folder_name,
                "R_init": R_init*r,
                "N_init": max(5, int(N_max*R_init))}
    return options

def param_no_growth(R, t_max, N_time, n=2, N=2000, output="output/"):
    coeff = 0. # no growth
    options = {'advection': "FLS",
                'n': 2,
                'delta': 1.,
                'eta': 1.,
                'psi0': 1.,
                'psiN': 0.6,
                'phi_init': 0.4,
                'K0': 1.,
                'delta': 1.,
                'sign': 1.,
                'BC': "dVdz==0",
                'coordinates': "spherical",
                "growth_rate_exponent": 0.5,
                'filename': 'IC_Sramek',
                'time_max': t_max,
                'dt_print': t_max/N_time,
                'coeff_velocity': coeff,
                'output': output,
                "R_init": R,
                "N_init": N, 
                "t_init": 0., 
                "Ric_adim": R}
    return options

def param_supercooling():
    return {}


def run_no_growth():
    Ric = 10**(np.linspace(0, 1, 3))
    for radius in Ric:
        output = "output_{:.2e}".format(radius)
        options = param_no_growth(radius, 2e3, 50, N=10000, output=output)
        run(options)


def run_growth():
    radius = 10**np.linspace(-1, 2, 1)# [100., 200., 300.]
    exponents = [1., 1.]
    coefficients = [1., 1.] #10**np.linspace(-2, 2, 5)#[1.]

    for r in radius:
        for exp in exponents:
            for coeff in coefficients:
                options = param_growth(r.item(), exp, coeff, basefolder="./test_param/", R_init=5e-3, N_max=2000)
                run(options)


def run_growth_random(Nr=20, Nc=20):
    random.seed()
    
    logradius = np.linspace(-3, 3, Nr)# [100., 200., 300.]
    dr = np.abs(np.diff(logradius)[0])
    exp= 1.
    logcoefficients = np.linspace(3, -4, Nc)#[1.]
    dc = np.abs(np.diff(logcoefficients)[0])

    n = 2

    for r in logradius:
        for coeff in logcoefficients:
            rand = random.normal([0.,0.], [dr/2., dc/2.])
            r = 10**(r+rand[0])
            coeff = 10**(coeff+rand[1])
            N_max = 2000
            if coeff < 1.:
                if r>900: N_max = 15000
                elif r> 200.: N_max = 5000
            options = param_growth(r.item(), exp, coeff, n=n, basefolder="./random_n2/", R_init=5e-3, N_max=N_max)
            run(options)





if __name__ == "__main__":

    #name_script = sys.argv[0]
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-v", "--verbose", help="increase output verbosity",
    #                action="store_true")
    #parser.add_argument("-i", "--input", help="choose input file (yaml file)", 
    #                action="store_true")
    #parser.add_argument("-o", "--output", help="specify output folder", action="store_true")
    #args = parser.parse_args()

    #if args.verbose:
    #    print("verbosity turned on")

    run_growth_random()