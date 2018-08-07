
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mush import *
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
	dt = 0.5*dr/V0

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


def test_velocity():


	options = {'advection':"FLS", \
				'Ra':0., \
				'K0':1., \
				'eta':1., \
				'delta':1., \
				'bc':'',\
				'psi0':0.4,
				's': 1}

	def velocity_analytical_check(ax, options):

		psi0 = options["psi0"]
		N = 1000
		R = np.linspace(0, 1, N-1)
		dr = R[1]-R[0]
		psi = psi0* np.ones(N)

		# from the inversion
		velocity = velocity_Sramek(1-psi, R, options)
		ax.plot(velocity, R, 'o')

		#analytical solution
		#R2 = np.linspace(0, 1, N+1)
		h = np.sqrt(options["delta"]**2 * psi0*(1-psi0)*(1+4/3*(1-psi0)))
		analytical_solution = -options["delta"]**2* psi0*(1-psi0)**2*\
								(1+ np.sinh(1/h)*np.sinh(R/h)/(np.cosh(1/h)+1)-np.cosh(R/h))
		ax.plot(analytical_solution, R, linewidth=2)
		return np.sum(velocity-analytical_solution)**2
		

	fig, ax = plt.subplots()
	for delta in np.linspace(0.1, 1.5, 10):

		options["delta"] = delta
		error = velocity_analytical_check(ax, options)
		print("error : {}".format(error))


if __name__ == "__main__":


	test_TDMA()
	Schema()
	advection_point() 
	#diffusion()
	test_velocity()
	plt.show()


	plt.show()
