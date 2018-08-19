import mush
import numpy as np
import matplotlib.pyplot as plt



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
               'delta': 1.,
               'eta': 1.,
               'phi0': 1.,
               'phiN': 0.5,
               'phi_init': 0.5,
               'sign': -1,
               'BC': "dVdz==0",
               'coordinates': "cartesian"}
    compaction_column_growth(velocity_Sramek, **options)

    plt.show()