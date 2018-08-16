
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mush import *
import compaction
import tdma

#TDMAsolver = tdma.TDMAsolver

def test_TDMA():
        
    a = np.zeros([10])
    b = np.ones([10])
    c = np.zeros([10])
    d = np.ones([10])
    x = inversion_matrice(a, b, c, d)
    #print('Test TDMA is OK if value is 0: {}'.format(np.sum(x-1)))	
    assert(np.sum(x-1)==0.)

    x = inversion_matrice(a, 2*b, c, d)
    assert(np.sum(x-0.5)==0.)
    print("10:0, 2, 0, 1: {}".format(x))

    x = inversion_matrice(a, b, c, 2*d)
    assert(np.sum(x-2)==0.)
    print("10:0, 1, 0, 2: {}".format(x))

    x = inversion_matrice(a, a, 2*a, d)
    print("10:1, 1, 2, 1: {}".format(x))

    x = inversion_matrice(2*a, a, c, d)
    print("10:2, 1, 0, 1: {}".format(x))

    a = np.zeros([3])
    b = np.ones([3])
    c = np.zeros([3])
    d = np.ones([3])

    x = inversion_matrice(a, a, a, d)
    print("3:1, 1, 1, 1: {}".format(x))

    x = inversion_matrice(2*a, a, c, d)
    print("3:2, 1, 0, 1: {}".format(x))

def test_fluxlimiterscheme():

    #input: velocity, variable (set the size of the output), dr, options={}
    #options['advection']: 'upwind' (default), "centered", "FLS"
    #output:_a/(2*dr), _b/(2*dr), _c/(2*dr), _d/(2*dr)

    #a, b, c, d = fluxlimiterscheme()
    pass


def advection_point():

    N = 200
    V0 = 1.
    V = V0 * np.ones([N-1])
    R = np.linspace(-2, 5, N)
    phi = np.zeros_like(R)
    phi_sin = np.where(np.abs(R)>1, 0, 1+np.cos(R*np.pi))
    phi_rec = np.where(np.abs(R)>1, 0, 1.)

    #phi[30:60] = 2.

    dr = R[1]-R[0]
    dt = np.abs(0.5*dr/V0)

    fig, ax = plt.subplots(3, 4)
#	fig, ax = plt.subplots(2, 1, sharex = True)

    def run(init, axis, scheme, correction_V=False):
        phi = init
        time = 0.
        axis.plot(R, init, 'k', linewidth=2)
        for it in range(1,100):
            phi_0 = phi
            phi= update(V, phi, dt, dr, scheme)
            time = time + dt
            if it%20 ==0:
                if correction_V:
                    correction = time * V[0]
                else: correction = 0.
                axis.plot(R-correction,phi)
                axis.set_title(scheme)

    run(phi_sin, ax[0,0], {'advection':"upwind"})
    run(phi_sin, ax[1,0], {'advection':"centered"})
    run(phi_sin, ax[2,0], {'advection':"FLS"})	
    run(phi_sin, ax[0,1], {'advection':"upwind"}, True)
    run(phi_sin, ax[1,1], {'advection':"centered"}, True)
    run(phi_sin, ax[2,1], {'advection':"FLS"}, True)
    run(phi_rec, ax[0,2], {'advection':"upwind"})
    run(phi_rec, ax[1,2], {'advection':"centered"})
    run(phi_rec, ax[2,2], {'advection':"FLS"})
    run(phi_rec, ax[0,3], {'advection':"upwind"}, True)
    run(phi_rec, ax[1,3], {'advection':"centered"}, True)
    run(phi_rec, ax[2,3], {'advection':"FLS"}, True)


def advection_gradient_velocity():
        
    N = 200
    V0 = 1.
    V = V0 * np.linspace(0,1,N-1)
    R = np.linspace(-2, 5, N)
    phi = np.zeros_like(R)
    phi_sin = np.where(np.abs(R)>1, 0, 1+np.cos(R*np.pi))
    #phi_rec = np.where(np.abs(R)>1, 0, 1.)

    #phi[30:60] = 2.

    dr = R[1]-R[0]
    dt = 0.01*dr/V0

    fig, ax = plt.subplots(3, 1)
#	fig, ax = plt.subplots(2, 1, sharex = True)

    def run(init, axis, scheme, correction_V=False):
        phi = init
        time = 0.
        axis.plot(R, init, 'k', linewidth=2)
        for it in range(1,2000):
            phi_0 = phi
            phi= update(V, phi, dt, dr, scheme)
            time = time + dt
            if it%1000 ==0:
                if correction_V:
                    correction = time * V[0]
                else: correction = 0.
                axis.plot(R-correction,phi)
                axis.set_title(scheme)
    run(phi_sin, ax[0], {'advection':"upwind"})
    run(phi_sin, ax[1], {'advection':"centered"})
    run(phi_sin, ax[2], {'advection':"FLS"})	


def diffusion():  ##maybe should change boundary conditions?

    def run(init, axis, scheme, correction_V=False):
        
        phi = init[:]
        time = 0.
        axis[0].plot(R, init, 'k', linewidth=2)
        for it in range(1,500):
            phi_0 = phi[:]
            phi = update(V, phi, dt, dr, scheme)
            time = time + dt
            if it%20 ==0:
                if correction_V:
                    correction = time * V[0]
                else: correction = 0.
                axis[0].plot(R-correction,phi)
                axis[0].set_title(scheme)
                axis[1].scatter(it, np.sum(phi))


    N = 200
    R = np.linspace(-2, 5, N)
    phi = np.zeros_like(R)
    V = np.zeros_like(R)
    phi_sin = np.where(np.abs(R)>1, 0, 1+np.cos(R*np.pi))
    phi_rec = np.where(np.abs(R)>1, 0, 1.)
    Ra = 1.
    dr = R[1]-R[0]
    dt = 0.01

    #phi = update(V, phi, dt, dr, {'Ra':Ra})
    fig, ax = plt.subplots(2,2)
    run(phi_sin, [ax[0,0],ax[0,1]], {'Ra':Ra, 'bc':'_'})
    run(phi_sin, [ax[1,0],ax[1,1]], {'Ra':Ra, 'bc':'dirichlet'})


def test_bc():
    N = 200
    V0 = 1.
    V = V0 * np.ones([N])
    R = np.linspace(-2, 5, N)
    phi = np.zeros_like(R)
    phi_sin = np.where(np.abs(R)>1, 0, 1+np.cos(R*np.pi))
    phi_rec = np.where(np.abs(R)>1, 0, 1.)

    dr = R[1]-R[0]
    dt = 0.5*dr/V0

    fig, ax = plt.subplots()
#	fig, ax = plt.subplots(2, 1, sharex = True)

    def run(init, axis, scheme, correction_V=False):
        
        phi = init
        time = 0.
        axis.plot(R, init, 'k', linewidth=2)
        for it in range(1,100):
            phi_0 = phi
            phi= update(V, phi, dt, dr, scheme)
            time = time + dt
            if it%20 ==0:
                if correction_V:
                    correction = time * V[0]
                else: correction = 0.
                axis.plot(R-correction,phi)
                axis.set_title(scheme)

    run(phi_sin, ax, {'advection':"upwind", "bc":"dirichlet"})

def analytical_solutions():
    """ Test the analytical solutions for the various boundary conditions. """
    options = {'advection':"FLS", 
                'phi0':0.3, 
                'delta':1., 
                's': 1}
    options["delta"] = 1./ np.sqrt(4/3*0.3/0.7**2)

    psi0 = 1- options["phi0"]
    N = 100
    R = np.linspace(0, 1, N+1)
    dr = R[1]-R[0]
    psi = psi0* np.ones(N)
    phi0 = 1-psi0

    fig, ax = plt.subplots(1,2) #cartesian and spherical solutions
    # function from Sumita
    options["BC"] = "V==0"
    options["coordinates"] = "cartesian"
    velocity = velocity_Sumita(1-psi, R, options)
    ax[0].plot(velocity, R[1:-1], label="Sumita_cart_V=0")
    velocity = compaction.analytic_Sumita_cart(phi0, R)
    ax[0].plot(velocity, R[:], '--', label="analytic_Sumita_cart_V=0")
    velocity = velocity_Sumita_spher(1-psi, R, options)
    ax[1].plot(velocity, R[1:-1], label="Sumita_spher_V=0")

    # function from Sramek
    options["coordinates"] = "cartesian"
    velocity = velocity_Sramek(1-psi, R, options)
    ax[0].plot(velocity, R[1:-1], label="Sramek_cart_V=0")
    velocity = compaction.analytic_Sramek_cart(phi0, R, options)
    ax[0].plot(velocity, R[:], '--', label="analytic_Sramek_cart_V=0")

    options["BC"] = "dVdz==0"
    velocity = velocity_Sramek(1-psi, R, options)
    ax[0].plot(velocity, R[1:-1], label="Sramek_cart_dVdz=0")
    velocity = compaction.analytic_Sramek_cart(phi0, R, options)
    ax[0].plot(velocity, R[:], '--', label="analytic_Sramek_cart_V=0")

    options["coordinates"] = "spherical"
    options["BC"] = "V==0"
    velocity = velocity_Sramek(1-psi, R, options)
    ax[1].plot(velocity, R[1:-1], label="Sramek_spher_V=0")
    velocity = compaction.analytic_Sramek_spher(phi0, R[1:], options)
    ax[1].plot(velocity, R[1:], '--', label="analytic_Sramek_spher_V=0")

    ax[0].legend()
    ax[1].legend()

if __name__ == "__main__":


    #test_TDMA()
    Schema()
    advection_point()
    #diffusion()
    advection_gradient_velocity()
    analytical_solutions()
    plt.show()
