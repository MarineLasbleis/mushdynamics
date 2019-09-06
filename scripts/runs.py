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
                'output': "compaction/",
                "R_init": 0.001,
                "N_init": 5,
                "t_init": 0.} # min(5, int(4e3/r_max))}
    return options

def modify_param(options, new_options):
    return {**options, **new_options}

def param_growth(r, exp, t_max, n=2, N_fig=20, basefolder="", R_init=1e-3, N_max=5000):
    dt = t_max/N_fig
    folder_name = basefolder+"/exp_{:.2e}_t_max_{:.2e}_radius_{:.2e}".format(exp, t_max, r)
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
                'output': folder_name,
                "R_init": R_init*r,
                "N_init": max(5, int(N_max*R_init)),
                "Ric_adim": r}
    return options

def param_no_growth(R, t_max, N_time, n=2, N=2000, output="output/"):
    coeff = 0. # no growth
    options = {'advection': "FLS",
                'n': n,
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

def param_supercooling(r, exp, t_max, r0_supercooling, n=2, N_fig=20, basefolder="", R_init=1e-3, N_max=5000):
    #t_max = (r/coeff)**(1/exp)
    dt = t_max/N_fig
    # r0 is the percentage of radius, not the actual radius!
    folder_name = basefolder+"/exp_{:.2e}_t_max_{:.2e}_radius_{:.2e}_r0_{}".format(exp, t_max, r, r0_supercooling)
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
                'output': folder_name,
                "R_init": R_init*r,
                "N_init": max(5, int(N_max*R_init)),
                "Ric_adim": r}
    options["t0_supercooling"] = 1e-3*t_max
    options["r0_supercooling"] = r0_supercooling/100.*r     # r0 is the percentage of radius, not the actual radius!
    return options


def run_no_growth(n=3):
    #Ric = 10**(np.linspace(-1, 2, 3))
    Ric = np.array([100.])
    for radius in Ric:
        output = "no_growth_n3/output_{:.2e}".format(radius)
        options = param_no_growth(radius.item(), 1e3, 10, N=5000, output=output, n=n)
        run(options)


def run_growth():
    radius = 10**np.array([-0.5]) #10**np.linspace(1., 3, 5 )# [100., 200., 300.]
    exponents = [0.5]
    times = 10**np.linspace(2, 0, 5)#[1.]

    for r in radius:
        for exp in exponents:
            for time in times:
                options = param_growth(r.item(), exp, time.item(), basefolder="./diag_n3_exp05/", R_init=5e-3, N_max=8000)
                run(options)


def run_growth_random(n=3, Nr=20, Nc=9):
    random.seed()
    logradius = np.array([2.75])# np.linspace(-3, 3, Nr)# [100., 200., 300.]
    dr = (0.5) #np.abs(np.diff(logradius)[0])
    print(logradius)
    exp= 0.5
    logcoefficients = np.linspace(2, -2, Nc)#[1.]
    dc = np.abs(np.diff(logcoefficients)[0])

    for logr in logradius:
        for logcoeff in logcoefficients:
            rand = random.normal([0.,0.], [dr/8., dc/8.])
            print(logr, logcoeff)
            r = 10**(logr+min(rand[0], dr/2. ))
            coeff = 10**(logcoeff+min(rand[1], dc/2.))
            print(r, coeff, rand)
            N_max = 2000
            if coeff < 1.:
                if r>900: N_max = 15000
                elif r> 200.: N_max = 5000
            options = param_growth(r.item(), exp, coeff.item(), n=n, basefolder="./diag_random_n{}_exp05/".format(n), R_init=5e-3, N_max=N_max)
            run(options)



def run_all_supercooling(Nr=1, Nc=9, N_r0=45):

    logradius = np.linspace(1., 1., Nr)
    exp= 0.5
    logcoefficients = np.array([]) #np.linspace(0, -2, Nc)#[1.]
    list_r0 = np.linspace(48, 10, N_r0)
    n = 3

    for logr in logradius:
        for logcoeff in logcoefficients:
            for r0 in list_r0:
                print(logr, logcoeff, r0, type(logr.item()))
                r = 10**(logr.item())
                coeff = 10**(logcoeff.item())
                N_max = 2000
                if coeff < 1.:
                    if r>900: N_max = 15000
                    elif r> 200.: N_max = 5000
                options = param_supercooling(r, exp, coeff, r0.item(), n=n, basefolder="./supercooling_4/", R_init=5e-3, N_max=N_max)
                run_supercooling(options)




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
