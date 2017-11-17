""" Reproduction of calculations from Sumita et 1996 on compaction and implications for Earth's inner core """



import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

def Schema():
	print(" <-1 volume-> \n \
|-----.-----|-----.-----|-----.-----|\n\
 0     0     1     1     2     2     3\n \
->          ->          ->          ->\n \
V0          V1          V2          V3 \n \
     phi0        phi1        phi2   \n \
      DP0         DP1         DP2\n ")


Nr = 10 #number of points in space
Nt = 10 #number of points in time




def fluxlimiterscheme(velocity, variable, dr, options={}): 
	""" output the coefficients for the advection scheme using the flux limiter scheme. 

	Coefficients for the scheme are lambda+, lambda-, v+ and v-
	(p and m are used instead of + and m)

	The code return the coefficients for the tridiagonal matrix a, b, c

	The scheme is used to solve the equation 
	D/Dt variable = D/Dr (variable * V)
	where D/Dr means partial derivative 

	Equation 3.54 in Sramek thesis gives:
	DF/Dx (at the point i) ~ 1/dx *(a_1 variable_i-1 + b_i variable_i +c_i variable_i+1)

	"""

	# detect size of the system and initialize variables:
	lambdap, lambdam, vp, vm = np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable)
	_a, _b, _c, _d = np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable)

	try: 
		option = options['advection']
	except KeyError:
		option = 'upwind'

	if option == 'upwind':
		pass #lambdap and lambdam == zeros. 
	elif option == "centered":
		lambdam[:] = np.ones_like(variable)
		lambdap[:] = np.ones_like(variable)
	elif option == "FLS":
		Rp=(variable[1:-2]-variable[0:-3])/(variable[2:-1]-variable[1:-2])
		Rm=(variable[3:]-variable[2:-1])  /(variable[2:-1]-variable[1:-2])
		#minmod
		lambdap[1:-2]=np.fmax(0.,np.minimum(1.,Rp))
		lambdam[1:-2]=np.fmax(0.,np.minimum(1.,Rm))

	else: 
		print("Problem with the choosen scheme for advection. Default is upwind.")

	vp[:-1] = 0.5*(velocity[:-1]+np.abs(velocity[:-1]))
	vm[:-1] = 0.5*(velocity[:-1]-np.abs(velocity[:-1]))
	vp[-1] = 0.
	vm[-1] = 0.

	_a[1:] = -vp[:-1]*(1-lambdap[:-1]/2.) - vm[:-1]*lambdam[:-1]/2.
	_b[1:] =  vp[1:]*(1-lambdap[1:]/2.) + vm[1:]*lambdam[1:]/2. \
				- vm[:-1]*(1-lambdam[:-1]/2.) - vp[:-1]*lambdap[:-1]/2.
	_c[1:] =  vm[1:]*(1-lambdam[1:]/2.) + vp[1:]*lambdap[1:]/2.

	_d[1:-1] = _a[1:-1]*variable[:-2]+_b[1:-1]*variable[1:-1]+_c[1:-1]*variable[2:]

	_a[0], _b[0], _c[0] = 0., 1., 0.

	return _a/(2*dr), _b/(2*dr), _c/(2*dr), _d/(2*dr)


def sourceterm():
	pass #for now, we can pass this! But useful for later. Will go into the "d". 

def CrankNicholson(variable, dr, options):
	""" Calculate the coefficients for the diffusion

	Variable is only here to provide the good shape for the arrays.
	"""

	# TODO : find another name for Ra (not the good name, as it is actually 1/Ra)
	try: 
		Ra = options["Ra"]
	except KeyError:
		Ra = 0.  # Default value is 0.

	_a, _b, _c, _d = np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable)
	diff = Ra/(dr**2)
	_a[:-1] = - diff/2.
	_b[:-1] = diff 
	_c[:-1] = - diff/2.
	_d[1:-1] = (variable[:-2] - 2*variable[1:-1] + variable[2:])/2.

	return _a, _b, _c, _d


def velocity_Sramek(variable, radius, options):  # Sramek thesis p46, equation 3.22

	dr = radius[1]-radius[0]

	try:
		s = options['s']
	except KeyError:
		s=1.
		print("s (sign of density difference) was not defined, please consider defining it for later. Default value is {}".format(s))
	
	try:
		K0 = options['K0']
	except KeyError:
		K0=1.
		print("K0 was not defined, please consider defining it for later. Default value is {}".format(K0))

	try: 
		delta = options["delta"]
	except KeyError:
		delta = .1
		print("Delta (compaction length) was not defined, please consider defining it for later. Default value is {}".format(delta))

	_inter = (1.+4/3*variable)*(1-variable)/variable

	_a, _b, _c, _d = np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable)

	_a[:-1] = _inter[:-1]/dr**2
	_b[:-1] = -1/delta**2/variable[:-1]/variable[1:] \
				-  _inter[:-1]/dr**2\
				-  _inter[1:] /dr**2
	_c[:-1] = _inter[1:]/dr**2
	_d[:-1] = s*(1-np.sqrt(variable[:-1]*variable[1:])) #if buoyancy/density variations, add terms here! s is 1 or -1.

	_a[-1], _b[-1], _c[-1], _d[-1] = 0,1,0,0

	new_velocity = np.zeros_like(variable)
	new_velocity[1:-1] = TDMAsolver(_a[1:-1], _b[1:-1], _c[1:-1], _d[1:-1])
	return new_velocity


def velocity_Sumita(variable, radius, options={}):

	# spherical symmetry
	# cartesian symmetry

	dr = radius[1]-radius[0] #assuming no variations of dr

	try:
		K0 = options['K0']
	except KeyError:
		K0=0.
		print("K0 was not defined, please consider defining it for later. Default value is {}".format(K0))

	try: 
		eta = options["eta"]
	except KeyError:
		eta = 1.
		print("eta was not defined, please consider defining it for later. Default value is {}".format(eta))

	_a, _b, _c, _d = np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable), np.zeros_like(variable)

	_a[:-1] = - K0 *eta /dr**2 * variable[:-1]
	_b[:-1] = 1 + K0*eta/dr**2 * (variable[1:]+variable[:-1])
	_c[:-1] = - K0 *eta /dr**2 * variable[1:]
	_d[:-1] = -K0 * np.sqrt(variable[:-1]*variable[1:])
	
	_a[-1], _b[-1], _c[-1], _d[-1] = 0,1,0,0

	new_velocity = TDMAsolver(_a, _b, _c, _d)

	print(new_velocity)

	return new_velocity


## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
	'''
	TDMA solver, a b c d can be NumPy array type or Python list type.
	refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	'''
	nf = len(a)     # number of equations
	ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array
	for it in range(1, nf):
		mc = ac[it]/bc[it-1]
		bc[it] = bc[it] - mc*cc[it-1] 
		dc[it] = dc[it] - mc*dc[it-1]
	xc = ac
	xc[-1] = dc[-1]/bc[-1]

	for il in range(nf-2, -1, -1):
		xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

	del bc, cc, dc  # delete variables from memory

	return xc



def update(V, phi, dt, dr, options = {'advection':"upwind", 'Ra':0.}):

		a_adv, b_adv, c_adv, d_adv = fluxlimiterscheme(V, phi, dr, options)
		a_diff, b_diff, c_diff, d_diff = CrankNicholson(phi, dr, options)
		_a, _b, _c, _d = a_adv+a_diff, b_adv+b_diff, c_adv+c_diff, d_adv+d_diff
	
		_a = _a*dt
		_b = 1.+_b*dt
		_c = _c*dt
		_d = phi-_d*dt
		#d = boundary_conditions(phi, a, b, c, d)
		# for the boundary conditions, we need to modify d[0] and d[-1]
		phi2 = phi[:]
		phi2 = np.zeros_like(phi)
		_phi = TDMAsolver(_a[1:-1], _b[1:-1], _c[1:-1], _d[1:-1])
		#phi2 = TDMAsolver(_a, _b, _c, _d)
		phi2[1:-1] = _phi
		phi2[0], phi2[-1] = phi[0], phi[-1]
		return phi2


def boundary_conditions(variable, a, b, c, d, option={}):

	d[0]  = d[0] -a[0] *variable[0]
	d[-1] = d[-1]-c[-1]*variable[-1] #Dirichlet

	return d


def Test_TDMA():
	a = np.zeros([10,1])
	b = np.ones([10,1])
	c = np.zeros([10,1])
	d = np.ones([10,1])

	x = TDMAsolver(a, b, c, d)
	print('Test TDMA is OK if value is 0: {}'.format(np.sum(x-1)))	


def advection_point():

	N = 200
	V0 = 1.
	V = V0 * np.ones([N])
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
		for it in range(1,150):
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
		axis.plot(R, init, 'k', linewidth=2)
		for it in range(1,150):
			phi_0 = phi[:]
			phi = update(V, phi, dt, dr, scheme)
			time = time + dt
			if it%20 ==0:
				if correction_V:
					correction = time * V[0]
				else: correction = 0.
				axis.plot(R-correction,phi)
				axis.set_title(scheme)

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
	fig, ax = plt.subplots()
	run(phi_sin, ax, {'Ra':Ra})


def compaction_column():


	options = {'advection':"FLS", 'Ra':0., 'K0':0.1, 'eta':1.}
	options = {'advection':"FLS", 'Ra':0., 'K0':0.1, 'eta':1., 'delta':0.1}

	psi0 = 0.4
	N = 100
	R = np.linspace(0, 1, N)
	dr = R[1]-R[0]
	psi = psi0* np.ones_like(R)
	psi[0] = 1.
	psi[-1] = 0.

	calcul_velocity = velocity_Sramek
	velocity = calcul_velocity(1-psi, R, options)
	v_m = np.amax(np.abs(velocity))
	dt = 0.5*dr/(v_m)

	fig, ax = plt.subplots(1,2, sharey=True)
	ax[0].plot(psi, R)
	ax[1].plot(velocity, R, 'o')

	h = np.sqrt(options["delta"]**2 * psi0*(1-psi0)*(1+4/3*(1-psi0)))
	analytical_solution = -options["delta"]**2* psi0*(1-psi0)**2*\
							(1+ np.sinh(1/h)*np.sinh(R/h)/(np.cosh(1/h)+1)-np.cosh(R/h))
	ax[1].plot(analytical_solution, R, linewidth=2)

	for it in range(1,2):
		psi = update(velocity, psi, dt, dr, options)
		#psi = np.where(psi>0, psi, 0)
		velocity = calcul_velocity(1-psi, R, options)
		v_m = np.amax(np.abs(velocity))
		dt = 0.5*dr/(v_m)
		if it%1==0 :
			ax[0].plot(psi, R)
			ax[1].plot(velocity, R)


	ax[0].set_xlim([0,1])


if __name__ == '__main__':


	#here is the main part of the code
	print('Sumita et al 1996, Geoph. J. Int.')
	Schema()

	#Test_TDMA()
	#advection_point()
	#iffusion()
	compaction_column()

	plt.show()

	