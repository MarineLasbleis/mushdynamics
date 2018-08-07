
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mush import *
import tdma




def test_velocity():


	options = {'advection':"FLS", \
				'Ra':0., \
				'K0':1., \
				'eta':1., \
				'delta':1., \
				'bc':'',\
				'psi0':0.1,
				's': 1,
				'grain':1}

	def velocity_analytical_check(ax, options):

		psi0 = options["psi0"]
		N = 100
		R = np.linspace(0, 1, N+1)
		dr = R[1]-R[0]
		psi = psi0* np.ones(N)
		phi0=1-psi0

		# from the inversion
		velocity = velocity_Sumita(1-psi, R, options)		
		#print("velocity = {}".format(velocity))
		ax.plot(velocity, R[1:-1], 'r')
		
		#analytical solution
		x1=np.sqrt((1+phi0)/(1-phi0))/phi0 * np.sqrt(3./4.)
		x2=-x1
		c3=-(phi0**3/((1+phi0)))
		c2=(c3*(np.exp(x1)-1))/(np.exp(x2)-np.exp(x1))
		c1=-c2-c3
		analytical_solution= c1*np.exp(x1*R) + c2*np.exp(x2*R) + c3
		#print("analytical_solution = {}".format(analytical_solution))
		ax.plot(analytical_solution, R, linewidth=2)
		return np.sum(velocity-analytical_solution[1:-1])**2

	_fig, ax = plt.subplots()

	velocity_analytical_check(ax, options)
	#for delta in np.linspace(0.1, 1.5, 10):

	#	options["delta"] = delta
	#	error = velocity_analytical_check(ax, options)
	#	print("error : {}".format(error))


if __name__ == "__main__":

	test_velocity()
	plt.show()

