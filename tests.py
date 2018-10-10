
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mush import *
import compaction
import tdma


def advection_point():

    N = 200
    V0 = 1.
    V = V0 * np.ones([N - 1])
    R = np.linspace(-2, 5, N)
    phi = np.zeros_like(R)
    phi_sin =  np.where(np.abs(R) > 1, 0, 1 + np.cos(R * np.pi))
    phi_rec = np.where(np.abs(R) > 1, 0, 1.)

    #phi[30:60] = 2.

    dr = R[1] - R[0]
    dt = np.abs(0.5 * dr / V0)

    fig, ax = plt.subplots(3, 4)
#	fig, ax = plt.subplots(2, 1, sharex = True)

    def run(init, axis, options, correction_V=False):
        phi = init
        time = 0.
        axis.plot(R, init, 'k', linewidth=2)
        for it in range(1, 100):
            phi_0 = phi
            phi = update(V, phi, dt, R, options)
            time = time + dt
            if it % 20 == 0:
                if correction_V:
                    correction = time * V[0]
                else:
                    correction = 0.
                axis.plot(R - correction, phi)
                axis.set_title(options["advection"])

    options = {'advection': "FLS",
               'psi0': 0.,
               'psiN': 0.,
               'delta': 1.,
               'sign': 1, 
               'coordinates': "cartesian",
               'BC': "dVdz==0"}
    options["advection"] = "upwind"
    run(phi_sin, ax[0, 0], options)
    options["advection"] = "centered"
    run(phi_sin, ax[1, 0], options)
    options["advection"] = "FLS"
    run(phi_sin, ax[2, 0], options)
    options["advection"] = "upwind"
    run(phi_sin, ax[0, 1], options, True)
    options["advection"] = "centered"
    run(phi_sin, ax[1, 1], options, True)
    options["advection"] = "FLS"
    run(phi_sin, ax[2, 1], options, True)
    options["advection"] = "upwind"
    run(phi_rec, ax[0, 2], options)
    options["advection"] = "centered"
    run(phi_rec, ax[1, 2], options)
    options["advection"] = "FLS"
    run(phi_rec, ax[2, 2], options)
    options["advection"] = "upwind"
    run(phi_rec, ax[0, 3], options, True)
    options["advection"] = "centered"
    run(phi_rec, ax[1, 3], options, True)
    options["advection"] = "FLS"
    run(phi_rec, ax[2, 3], options, True)


def analytical_solutions():
    """ Test the analytical solutions for the various boundary conditions. """
    options = {'advection': "FLS",
               'phi0': 0.3,
               'delta': 1.,
               'sign': -1, 
               'Ric_adim': 1.}
    options["delta"] = 1. / np.sqrt(4 / 3 * 0.3 / 0.7**2)

    psi0 = 1 - options["phi0"]
    N = 100
    R = np.linspace(0, 1, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)
    phi0 = 1 - psi0

    fig, ax = plt.subplots(1, 2)  # cartesian and spherical solutions
    # function from Sumita
    options["BC"] = "V==0"
    options["coordinates"] = "cartesian"
    velocity = velocity_Sumita(1 - psi, R, options)
    ax[0].plot(velocity, R[1:-1], label="Sumita_cart_V=0")
    velocity = compaction.analytic_Sumita_cart(phi0, R, options)
    ax[0].plot(velocity, R[:], '--', label="analytic_Sumita_cart_V=0")
    options["coordinates"] = "spherical"
    velocity = velocity_Sumita(1 - psi, R, options)
    ax[1].plot(velocity, R[1:-1], label="Sumita_spher_V=0")

    # function from Sramek
    options["coordinates"] = "cartesian"
    velocity = velocity_Sramek(1 - psi, R, options)
    ax[0].plot(velocity, R[1:-1], label="Sramek_cart_V=0")
    velocity = compaction.analytic_Sramek_cart(phi0, R, options)
    ax[0].plot(velocity, R[:], '--', label="analytic_Sramek_cart_V=0")

    options["BC"] = "dVdz==0"
    velocity = velocity_Sramek(1 - psi, R, options)
    ax[0].plot(velocity, R[1:-1], label="Sramek_cart_dVdz=0")
    velocity = compaction.analytic_Sramek_cart(phi0, R, options)
    ax[0].plot(velocity, R[:], '--', label="analytic_Sramek_cart_V=0")
    options["coordinates"] = "spherical"
    options["BC"] = "V==0"
    velocity = velocity_Sramek(1 - psi, R, options)
    ax[1].plot(velocity, R[1:-1], label="Sramek_spher_V=0")
    velocity = compaction.analytic_Sramek_spher(phi0, R[1:], options)
    ax[1].plot(velocity, R[1:], '--', label="analytic_Sramek_spher_V=0")

    ax[0].legend()
    ax[1].legend()


def advection_Vcst():

    N = 200
    V0 = 1.
    V = V0 * np.ones([N - 1])
    R = np.linspace(-5, 5, N)
    phi = np.zeros_like(R)
    phi_sin = 0.5+ np.where(np.abs(R) > 1, 0, 1 + np.cos(R * np.pi))
    phi_rec = np.where(np.abs(R) > 1, 0, 1.)

    options = {'advection': "FLS",
               'psi0': 0.3,
               'psiN': 0.8,
               'delta': 1.,
               'sign': 1,
               'coordinates': "cartesian", 
               'BC': "dVdz==0"}

    dr = R[1] - R[0]
    dt = np.abs(0.5 * dr / V0)

    fig, ax = plt.subplots()

    def run(init, axis, options, correction_V=False):
        phi = init
        time = 0.
        axis.plot(R, init, 'k', linewidth=2)
        for it in range(1, 100):
            phi_0 = phi
            phi = update(V, phi, dt, R, options)
            time = time + dt
            if it % 20 == 0:
                if correction_V:
                    correction = time * V[0]
                else:
                    correction = 0.
                axis.plot(R - correction, phi)

    run(phi_sin, ax, options)
    ax.plot(0., 0.3, 'o')

def test_output():
    """ Test function output in mush """
    options = {'advection': "FLS",
               'phi0': 0.3,
               'delta': 1.,
               'sign': -1}
    options["delta"] = 1. / np.sqrt(4 / 3 * 0.3 / 0.7**2)

    psi0 = 1 - options["phi0"]
    N = 100
    R = np.linspace(0, 1, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)
    phi0 = 1 - psi0

    fig, ax = plt.subplots(1, 2)  
    # function from Sumita
    options["BC"] = "V==0"
    options["coordinates"] = "cartesian"
    velocity = velocity_Sumita(1 - psi, R, options)
    data = {"radius": pd.Series(R), 'porosity': pd.Series(1-psi), 'velocity': pd.Series(velocity)}
    data = pd.DataFrame(data)
    output(0, data, True, True, "test", ax)
    velocity = compaction.analytic_Sumita_cart(phi0, R, options)
    data = {"radius": pd.Series(R), 'porosity': pd.Series(1-psi), 'velocity': pd.Series(velocity)}
    data = pd.DataFrame(data)
    output(1, data, True, True, "test", ax)
    # output(time, psi, velocity, R, fig=False, file=False, output_folder="", ax=[]):


if __name__ == "__main__":

    # test_TDMA()
    schema()
    advection_point()
    analytical_solutions()
    advection_Vcst()
    test_output()
    plt.show()
