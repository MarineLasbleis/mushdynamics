""" Compaction in a layer as in Sramek phd thesis   """


import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from mush import *


def compaction_column(calcul_velocity, **options):
    """ Compaction of a column of sediments. No velocity up or down.
    
    Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    psi0 = 1 - options["phi_init"]
    N = 1000
    R = np.linspace(0, 1, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)

    velocity = calcul_velocity(1 - psi, R, options)
    v_m = np.amax(np.abs(velocity))
    dt = min(0.5 * dr / (v_m), 0.5)

    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].plot(1 - psi, R[:-1] + dr / 2.)
    ax[1].plot(velocity, R[1:-1])

    time = 0.
    dt_print = 2.
    time_p = time
    time_max = 4000.
    it = 0
    iter_max = 5000

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        psi = update(velocity, psi, dt, R, options)
        velocity = calcul_velocity(1 - psi, R, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.1 * dr / (v_m))
        # if time_p > dt_print:
        if it % 1000 == 0:
            print(it, dt, time)
            # reinitinalize the mark to know if we need to print/plot
            # something.
            time_p = time_p - dt_print
            ax[0].plot(1 - psi, R[:-1] + dr / 2.)
            ax[1].plot(velocity, R[1:-1])
    print(it)
    ax[0].set_xlim([0, 1])
    ax[0].set_ylim([0, 1])
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")



def compaction_column_dVdz():
    """ Compaction of a column of sediments.

    Liquids can go through top (dVdz =0). Top boundary set at phi0.
    Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    calcul_velocity = velocity_Sumita
    options = {'advection': "upwind",
               'Ra': 0.,
               'eta': 1.,
               'bc': '',
               'phi0': 0.,
               'phiN': 0.6,
               'sign': 1,
               'BC': "dVdz==0",
               'coordinates': "cartesian"}
    psi0 = 0.6 #1 - options["phi0"]
    N = 1000
    R = np.linspace(0, 1, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)

    velocity = calcul_velocity(1 - psi, R, options)
    v_m = np.amax(np.abs(velocity))
    dt = min(0.5 * dr / (v_m), 0.5)
    print(dt)

    fig, ax = plt.subplots(1, 2, sharey=True)
    fig2, ax2 = plt.subplots(2)
    ax[0].plot(1 - psi, R[:-1] + dr / 2.)
    ax[1].plot(velocity, R[1:-1])

    time = 0.
    dt_print = 2.
    time_p = time
    time_max = 4000.
    it = 0
    iter_max = 10000
    ax2[0].plot(time, sum_phi(1-psi), '+')
    ax2[1].plot(time, flux_top(1-psi, velocity), '+')

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        psi = update(velocity, psi, dt, R, options)
        velocity = calcul_velocity(1 - psi, R, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.1 * dr / (v_m))
        ax2[0].plot(time, sum_phi(1-psi), '+')
        ax2[1].plot(time, flux_top(1-psi, velocity), '+')
        # if time_p > dt_print:
        if it % 1000 == 0:
            print(it, dt, time)
            # reinitinalize the mark to know if we need to print/plot
            # something.
            time_p = time_p - dt_print
            ax[0].plot(1 - psi, R[:-1] + dr / 2.)
            ax[1].plot(velocity, R[1:-1])
    print(it)
    ax[0].set_xlim([0, 1])
    ax[0].set_ylim([0, 1])
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")

def sum_phi(phi):
    return np.sum(phi)/len(phi) #in cartesian
def flux_top(phi, velocity):
    return (1-phi[-1])*velocity[-1]

def analytic_Sumita_cart(phi0, R):
    """ Solution analytique pour resolution Sumita in cartesian coordinates. """
    x1 = np.sqrt(1 / phi0**2) * np.sqrt(3. / 4.)
    x2 = -x1
    c3 = -(phi0**3 / ((1 - phi0)))
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





def figures_compaction_only():
    options = {'advection': "upwind",
               'Ra': 0.,
               'eta': 1.,
               'bc': '',
               'phi0': 0.,
               'phiN': 1.,
               'U0': 0.,
               'UN': 0.,
               'sign': -1,
               'BC': "V==0",
               'coordinates': "cartesian"}

    compaction_column(velocity_Sramek, delta=1., **options)
    plt.savefig("fig/phi03_Sramek_delta_1.pdf")
    compaction_column(velocity_Sramek, delta=.5, **options)
    plt.savefig("fig/phi03_Sramek_delta_05.pdf")
    compaction_column(velocity_Sramek, delta=.2, **options)
    plt.savefig("fig/phi03_Sramek_delta_02.pdf")
    compaction_column(velocity_Sramek, delta=.1, **options)
    plt.savefig("fig/phi03_Sramek_delta_01.pdf")
    compaction_column(velocity_Sramek, delta=.05, **options)
    plt.savefig("fig/phi03_Sramek_delta_005.pdf")

    compaction_column(velocity_Sumita, K=1., **options)
    plt.savefig("fig/phi03_Sumita_K_1.pdf")
    compaction_column(velocity_Sumita, K=.5, **options)
    plt.savefig("fig/phi03_Sumita_K_05.pdf")
    compaction_column(velocity_Sumita, K=.2, **options)
    plt.savefig("fig/phi03_Sumita_K_02.pdf")
    compaction_column(velocity_Sumita, K=.1, **options)
    plt.savefig("fig/phi03_Sumita_K_01.pdf")
    compaction_column(velocity_Sumita, K=.05, **options)
    plt.savefig("fig/phi03_Sumita_K_005.pdf")

    options["coordinates"] = "spherical"
    compaction_column(velocity_Sramek, delta=1., **options)
    plt.savefig("fig/phi03_Sramek_delta_1_sph.pdf")
    compaction_column(velocity_Sramek, delta=.5, **options)
    plt.savefig("fig/phi03_Sramek_delta_05_sph.pdf")
    compaction_column(velocity_Sramek, delta=.2, **options)
    plt.savefig("fig/phi03_Sramek_delta_02_sph.pdf")
    compaction_column(velocity_Sramek, delta=.1, **options)
    plt.savefig("fig/phi03_Sramek_delta_01_sph.pdf")
    compaction_column(velocity_Sramek, delta=.05, **options)
    plt.savefig("fig/phi03_Sramek_delta_005_sph.pdf")

    compaction_column(velocity_Sumita_spher, K=1., **options)
    plt.savefig("fig/phi03_Sumita_K_1_sph.pdf")
    compaction_column(velocity_Sumita_spher, K=.5, **options)
    plt.savefig("fig/phi03_Sumita_K_05_sph.pdf")
    compaction_column(velocity_Sumita_spher, K=.2, **options)
    plt.savefig("fig/phi03_Sumita_K_02_sph.pdf")
    compaction_column(velocity_Sumita_spher, K=.1, **options)
    plt.savefig("fig/phi03_Sumita_K_01_sph.pdf")
    compaction_column(velocity_Sumita_spher, K=.05, **options)
    plt.savefig("fig/phi03_Sumita_K_005_sph.pdf")




def compaction_column_growth(calcul_velocity, **options):
    """ Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    psi0 = 1 - options["phi_init"]
    N = 10
    R_init = 0.001
    R = np.linspace(0, R_init, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)

    velocity = calcul_velocity(1 - psi, R, options)
    v_m = np.amax(np.abs(velocity))
    dt = min(0.5 * dr / (v_m), 0.5)
    print(dt)

    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].plot(1 - psi, R[:-1] + dr / 2.)
    ax[1].plot(velocity, R[1:-1])

    time = 0.
    dt_print = 2.
    time_p = time
    time_max = 50.
    it = 0
    iter_max = 10000

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        if R[-1]+dr < radius(time, R_init):
            psi, R = append_radius(psi, R, options)
        #psi = np.append(psi, [psi0, psi0])
        #R = np.append(R, [R[-1]+dr, R[-1]+2*dr])
        velocity = calcul_velocity(1 - psi, R, options)
        psi = update(velocity, psi, dt, R, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.001 * dr / (v_m))
        # if time_p > dt_print:
        if it % 1000 == 0:
            print(it, dt, time, R[-1])
            # reinitinalize the mark to know if we need to print/plot
            # something.
            time_p = time_p - dt_print
            ax[0].plot(1 - psi, R[:-1] + dr / 2.)
            ax[1].plot(velocity, R[1:-1])
    print(it)

    #ax[0].set_xlim([0.3, 0.7])
    # ax[0].set_ylim([0,1])
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")


def radius(time, R_init):
    return R_init+ 1.* time


def append_radius(psi, R, options):
    psi = np.append(psi, [1-options["phiN"]])
    dr = R[1] - R[0]
    R = np.append(R, [R[-1] + dr])
    return psi, R


if __name__ == "__main__":

    options = {'advection': "FLS",
               'Ra': 0.,
               'eta': 1.,
               'phi0': 1.,
               'phiN': 0.,
               'phi_init': 0.3,
               'sign': -1,
               'BC': "V==0",
               'coordinates': "spherical"}

    # compaction_column(velocity_Sramek, delta=1., **options)

    options = {'advection': "FLS",
               'Ra': 0.,
               'delta': 1.,
               'eta': 1.,
               'phi0': 1.,
               'phiN': 0.5,
               'phi_init': 0.5,
               'sign': -1,
               'BC': "dVdz==0",
               'coordinates': "cartesian"}
    compaction_column_growth(velocity_Sramek, **options)
    #compaction_column_dVdz()
    #figures_compaction_only()
    

    
    
    plt.show()
