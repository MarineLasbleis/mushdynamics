""" Compaction in a layer as in Sramek phd thesis   """


import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd

from mush import *


def compaction_column(calcul_velocity, output_fig=True, **options):
    """ Compaction of a column of sediments.

    Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) 
    if output_fig False: output data files and no figures (OK for cluster)
    """

    psi0 = 1 - options["phi_init"]
    N = 1000
    R = np.linspace(0, 1, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)

    velocity = calcul_velocity(1 - psi, R, options)
    v_m = np.amax(np.abs(velocity))
    dt = min(0.5 * dr / (v_m), 0.5)

    if output_fig:
        fig, ax = plt.subplots(1, 4)
        ax[0].plot(1 - psi, R[:-1] + dr / 2.)
        ax[1].plot(velocity, R[1:-1])

    time = 0.
    dt_print = 20.
    time_p = time
    time_max = 400.
    it = 0
    iter_max = 2000

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        psi = update(velocity, psi, dt, R, options)
        velocity = calcul_velocity(1 - psi, R, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.1 * dr / (v_m))
        if time_p > dt_print:
            time_p = time_p - dt_print
        # if it % 1000 == 0:
            if output_fig:
                print(it, dt, time)
                # reinitinalize the mark to know if we need to print/plot
                # something.
                ax[0].plot(1 - psi, R[:-1] + dr / 2.)
                ax[1].plot(velocity, R[1:-1])
                ax[2].plot(sum_phi(1-psi, R[1:], options), time, 'x')
                ax[3].plot(velocity[-1], time, '+')
            else:
                file = "output/output_{:5.2f}.csv".format(time)
                _data = {"radius": pd.Series(R), 'porosity': pd.Series(1-psi), 'velocity': pd.Series(velocity)}
                data = pd.DataFrame(_data)
                data.to_csv(file)


    if output_fig:
        ax[0].set_xlim([0, 1])
        ax[0].set_ylim([0, 1])
        ax[1].set_ylim([0, 1])
        ax[0].set_xlabel("Porosity")
        ax[0].set_ylabel("Height (non-dim)")
        ax[1].set_xlabel("Solid velocity (non-dim)")
        ax[2].set_ylabel("time (un dim)")
        ax[2].set_xlabel("Total porosity")
        ax[2].set_xlim([0, 1])
        ax[3].set_xlabel("Velocity at top")

def sum_phi(phi, R, options):
    dr = R[1]-R[0] # constant steps in radius/heigth
    if options["coordinates"] == "cartesian":
        dV = dr*np.ones_like(R)
    elif options["coordinates"] == "spherical":
        dV = 4*np.pi*dr*R**2
    return np.sum(phi*dV)/np.sum(dV) #in cartesian

def flux_top(phi, velocity):
    return (1-phi[-1])*velocity[-1]


def analytic_Sumita_cart(phi0, R, options):
    """ Solution analytique pour resolution Sumita in cartesian coordinates. """
    s = options["sign"]
    x1 = np.sqrt(1 / phi0**2) * np.sqrt(3. / 4.)
    x2 = -x1
    c3 = -s*(phi0**3 / ((1 - phi0)))
    c2 = (c3 * (np.exp(x1) - 1)) / (np.exp(x2) - np.exp(x1))
    c1 = -c2 - c3
    return c1 * np.exp(x1 * R) + c2 * np.exp(x2 * R) + c3


def analytic_Sramek_cart(phi0, R, options):
    """ Solution analytique pour resolution Sramek in cartesian coordinates. """
    h = np.sqrt(options["delta"]**2 * phi0 *
                (1. + 4. / 3. * phi0) * (1 - phi0))
    C3 = -options["sign"] * phi0**2 * (1 - phi0) * options["delta"]**2
    if options["BC"] == "V==0":
        C1 = C3 * (np.exp(-1 / h) - 1) / (np.exp(1 / h) - np.exp(-1 / h))
        C2 = -C1 - C3
    elif options["BC"] == "dVdz==0":
        C1 = -C3 / (1 + np.exp(2 / h))
        C2 = C1 * np.exp(2 / h)
    analytical_solution = C1 * np.exp(R / h) + C2 * np.exp(-R / h) + C3
    return analytical_solution


def analytic_Sramek_spher(phi0, R, options):
    h = np.sqrt(options["delta"]**2 * phi0 *
                (1. + 4. / 3. * phi0) * (1 - phi0))
    C3 = -options["sign"] * phi0**2 * (1 - phi0) * options["delta"]**2
    div = np.cosh(1 / h) - h * np.sinh(1 / h)
    op = (np.cosh(R / h) - h / R * np.sinh(R / h)) / R
    return (R - op / div) * C3





def figures_compaction_only(output="/fig"):

    options = {'advection': "FLS",
               'eta': 1.,
               'phi0': 0.,
               'phiN': 1.,
               'phi_init': 0.3,
               'U0': 0.,
               'UN': 0.,
               'sign': 1,
               'BC': "V==0",
               'coordinates': "cartesian"}

    # top boundary impermeable
    compaction_column(velocity_Sramek, delta=1., **options)
    plt.savefig("fig/phi03_Sramek_delta_1.pdf")
    compaction_column(velocity_Sramek, delta=.5, **options)
    plt.savefig("fig/phi03_Sramek_delta_05.pdf")
    compaction_column(velocity_Sramek, delta=.1, **options)
    plt.savefig("fig/phi03_Sramek_delta_01.pdf")

    compaction_column(velocity_Sumita, K=1., **options)
    plt.savefig("fig/phi03_Sumita_K_1.pdf")
    compaction_column(velocity_Sumita, K=.5, **options)
    plt.savefig("fig/phi03_Sumita_K_05.pdf")


    options["coordinates"] = "spherical"
    compaction_column(velocity_Sramek, delta=1., **options)
    plt.savefig("fig/phi03_Sramek_delta_1_sph.pdf")
    compaction_column(velocity_Sramek, delta=.5, **options)
    plt.savefig("fig/phi03_Sramek_delta_05_sph.pdf")
    compaction_column(velocity_Sramek, delta=.1, **options)
    plt.savefig("fig/phi03_Sramek_delta_01_sph.pdf")


    compaction_column(velocity_Sumita, K=1., **options)
    plt.savefig("fig/phi03_Sumita_K_1_sph.pdf")
    compaction_column(velocity_Sumita, K=.5, **options)
    plt.savefig("fig/phi03_Sumita_K_05_sph.pdf")

    ## top boundary free
    options["BC"] = "dVdz==0"
    options["coordinates"] = "cartesian"
    options["phiN"] = options["phi_init"]

    compaction_column(velocity_Sramek, delta=1., **options)
    plt.savefig("fig/top_phi03_Sramek_delta_1.pdf")
    compaction_column(velocity_Sramek, delta=.5, **options)
    plt.savefig("fig/top_phi03_Sramek_delta_05.pdf")
    compaction_column(velocity_Sramek, delta=.1, **options)
    plt.savefig("fig/top_phi03_Sramek_delta_01.pdf")

    compaction_column(velocity_Sumita, K=1., **options)
    plt.savefig("fig/top_phi03_Sumita_K_1.pdf")
    compaction_column(velocity_Sumita, K=.5, **options)
    plt.savefig("fig/top_phi03_Sumita_K_05.pdf")

    options["coordinates"] = "spherical"
    compaction_column(velocity_Sramek, delta=1., **options)
    plt.savefig("fig/top_phi03_Sramek_delta_1_sph.pdf")
    compaction_column(velocity_Sramek, delta=.5, **options)
    plt.savefig("fig/top_phi03_Sramek_delta_05_sph.pdf")
    compaction_column(velocity_Sramek, delta=.1, **options)
    plt.savefig("fig/top_phi03_Sramek_delta_01_sph.pdf")

    compaction_column(velocity_Sumita, K=1., **options)
    plt.savefig("fig/top_phi03_Sumita_K_1_sph.pdf")
    compaction_column(velocity_Sumita, K=.5, **options)
    plt.savefig("fig/top_phi03_Sumita_K_05_sph.pdf")





if __name__ == "__main__":

    options = {'advection': "FLS",
               'Ra': 0.,
               'eta': 1.,
               'phi0': 1.,
               'phiN': 0.,
               'phi_init': 0.3,
               'sign': 1,
               'BC': "V==0",
               'coordinates': "spherical", 
               "Ric_adim": 1.}

    compaction_column(velocity_Sramek, output_fig=False, delta=.1, **options)


    #compaction_column_dVdz()
    #figures_compaction_only()

    plt.show()
