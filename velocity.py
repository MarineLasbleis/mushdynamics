import numpy as np
import scipy as sc
from scipy.sparse import diags
from scipy.sparse.linalg import inv
import matplotlib.pyplot as plt

from tdma import *
from Sumita96 import velocity_Sramek




fig, ax = plt.subplots(2)
delta = 1.
K0 = 1.
s = 1.
psi0 = 0.5

z = np.linspace(0, 1, 200)
dr = z[1]-z[0]
print(dr)
a, b, c, d = np.zeros_like(z), np.zeros_like(z), np.zeros_like(z), np.zeros_like(z)
vel = np.zeros_like(z)

variable = psi0*np.ones_like(z)
#variable[0] = 0.
#variable[-1] = 1.

_inter = (1.+4./3.*variable)*(1.-variable)/variable


a[1:] = _inter[:-1]/dr**2
b[1:] = -1./(delta**2*variable[:-1]*variable[1:]) \
				-  _inter[:-1]/dr**2\
				-  _inter[1:] /dr**2
c[1:] = _inter[1:]/dr**2
d[1:] = s*(1-np.sqrt(variable[:-1]*variable[1:])) #if buoyancy/density variations, add terms here! s is 1 or -1.

a[-1], b[-1], c[-1], d[-1] = 0,1,0,0
a[0], b[0], c[0], d[0] = 0,1,0,0


#vel[1:-1] = TDMAsolver(a[1:-1], b[1:-1], c[1:-1], d[1:-1])
vel = inversion_matrice(a[1:], b, c[:-1], d)
print("length of velocity: {}".format(vel))


vel2 = velocity_Sramek(variable, z, {})


h = np.sqrt(delta**2 * psi0*(1-psi0)*(1+4/3*(psi0)))
analytical_solution = -s*delta**2* psi0**2*(1-psi0)*\
							(1+ np.sinh(1/h)*np.sinh(z/h)/(np.cosh(1/h)+1)-np.cosh(z/h))
ax[1].plot(analytical_solution, z, linewidth=2)
ax[1].plot(vel, z, 'ro')
ax[1].plot(vel2, z, 'g-')
ax[0].plot(variable, z)

print(vel)

plt.show()

