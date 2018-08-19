import mush
import numpy as np
import matplotlib.pyplot as plt



def compaction_column_growth(calcul_velocity, **options):
    """ Calcul_Velocity is a function (velocity_Sramek or velocity_Sumita) """

    psi0 = 1 - options["phi_init"]
    N = 10
    R_init = radius(options["t_init"], options)
    R = np.linspace(0, R_init, N + 1)
    dr = R[1] - R[0]
    psi = psi0 * np.ones(N)
    time = options["t_init"]
    dt_print = .1
    time_p = time
    time_max = 1.
    it = 0
    iter_max = 100000

    velocity = calcul_velocity(1 - psi, R, options)
    v_m = np.amax(np.abs(velocity))
    dt = min(0.5 * dr / (v_m), 0.5)
    dt = min(dt, dr/growth_rate(time, options))
    print(dt)

    fig, ax = plt.subplots(1, 4, sharey=True)
    ax[0].plot(1 - psi, R[:-1] + dr / 2.)
    ax[1].plot(velocity, R[1:-1])



    while time < time_max and it < iter_max:
        # for it in range(0,10000):
        it = it + 1
        time = time + dt
        time_p = time_p + dt
        if R[-1]+dr < radius(time, options):
            psi, R = append_radius(psi, R, options)
        #psi = np.append(psi, [psi0, psi0])
        #R = np.append(R, [R[-1]+dr, R[-1]+2*dr])
        velocity = calcul_velocity(1 - psi, R, options)
        psi = mush.update(velocity, psi, dt, R, options)
        v_m = np.amax(np.abs(velocity))
        dt = min(0.5, 0.001 * dr / (v_m))
        dt = min(dt, 0.5*dr/growth_rate(time, options))

        if time_p > dt_print:
        # if it % 100 == 0:
            print(it, dt, time, R[-1], len(R))
            # reinitinalize the mark to know if we need to print/plot
            # something.
            time_p = time_p - dt_print
            ax[0].plot(1 - psi, R[:-1] + dr / 2.)
            ax[1].plot(velocity, R[1:-1])
            ax[2].plot(sum_phi(1-psi, R[1:], options), R[-1], 'x')
            ax[3].plot((options["psiN"])*velocity[-1], R[-1], 'x')
    print(it)

    #ax[0].set_xlim([0.3, 0.7])
    # ax[0].set_ylim([0,1])
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")



def radius(time, options):
    return (time)**options["growth rate exponent"]

def growth_rate(time, options):
    return time**(1-options["growth rate exponent"])

def append_radius(psi, R, options):
    psi = np.append(psi, [options["psiN"]])
    dr = R[1] - R[0]
    R = np.append(R, [R[-1] + dr])
    return psi, R


def sum_phi(phi, R, options):
    dr = R[1]-R[0] # constant steps in radius/heigth
    if options["coordinates"] == "cartesian":
        dV = dr*np.ones_like(R)
    elif options["coordinates"] == "spherical":
        dV = 4*np.pi*dr*R**2
    return np.sum(phi*dV)/np.sum(dV) #in cartesian

def flux_top(phi, velocity):
    return (1-phi[-1])*velocity[-1]


if __name__ == "__main__":

    options = {'advection': "FLS",
               'Ra': 0.,
               'delta': 1.,
               'eta': 1.,
               'psi0': 1.,
               'psiN': 0.5,
               'phi_init': 0.5,
               'sign': 1,
               'BC': "dVdz==0",
               'coordinates': "cartesian",
               "t_init": 0.01,
               "growth rate exponent": 1}
    compaction_column_growth(mush.velocity_Sramek, **options)

    plt.show()