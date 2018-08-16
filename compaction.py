""" Compaction in a layer as in Sramek phd thesis   """


import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from mush import *


def compaction_column(calcul_velocity, **options):
    """ Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    psi0 = 1 - options["phi0"]
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
    iter_max = 100000

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        psi = update(velocity, psi, dt, dr, options)
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
    C3 = -options["s"] * phi0**2 * (1 - phi0) * options["delta"]**2
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
    C3 = -options["s"] * phi0**2 * (1 - phi0) * options["delta"]**2
    div = np.cosh(1 / h) - h * np.sinh(1 / h)
    op = (np.cosh(R / h) - h / R * np.sinh(R / h)) / R
    return (R - op / div) * C3


def compaction_column_growth(calcul_velocity, **options):
    """ Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    psi0 = 1 - options["phi0"]
    N = 1000
    R = np.linspace(0, 0.01, N + 1)
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
    iter_max = 500

    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        # if R[-1]+dr < radius(time):
        #    psi, R = append_radius(psi, R)
        #psi = np.append(psi, [0.5, 0.5])
        #R = np.append(R, [R[-1]+dr, R[-1]+2*dr])
        velocity = calcul_velocity(1 - psi, R, options)
        psi = update(velocity, psi, dt, dr, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.1 * dr / (v_m))
        # if time_p > dt_print:
        if it % 10 == 0:
            print(it, dt, time)
            # reinitinalize the mark to know if we need to print/plot
            # something.
            time_p = time_p - dt_print
            ax[0].plot(1 - psi, R[:-1] + dr / 2.)
            ax[1].plot(velocity, R[1:-1])
    print(it)

    ax[0].set_xlim([0.3, 0.7])
    # ax[0].set_ylim([0,1])
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")


def radius(time):
    return 1. * time


def append_radius(psi, R):
    psi = np.append(psi, [0.5])
    dr = R[1] - R[0]
    R = np.append(R, [R[-1] + dr])
    return psi, R


if __name__ == "__main__":

    options = {'advection': "upwind",
               'Ra': 0.,
               'eta': 1.,
               'bc': '',
               'phi0': .3,
               'phiN': 0.,
               'U0': 0.,
               'UN': 0.,
               'sign': 1,
               'BC': "dVdz==0",
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

    # compaction_column(velocity_Sumita, K0=1, **options)
    # plt.savefig("fig/Sumita_K_1.pdf")
    # compaction_column(velocity_Sumita, K0=10, **options)
    # plt.savefig("fig/Sumita_K_10.pdf")
    # compaction_column(velocity_Sumita, K0=0.1, **options)
    # plt.savefig("fig/Sumita_K_01.pdf")

    compaction_column_growth(velocity_Sumita, delta=1., **options)
    plt.show()
