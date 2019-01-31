import numpy as np
import matplotlib.pyplot as plt
from numpy.lib import scimath

def function_n2(zeta):
    """ From analytical solution:

    x = 1/2 ((sqrt(3) sqrt(27 A^2 - 1) + 9 A)^(1/3)/3^(2/3) + 1/(3^(1/3) (sqrt(3) sqrt(27 A^2 - 1) + 9 A)^(1/3)) + 1)
    """
    racine_3 = scimath.power((3**0.5*scimath.power((27.*zeta**2.-1), 0.5)+9.*zeta), 1./3.)
    return np.real(0.5*(racine_3/3.**(2./3.) + 1/racine_3/3.**(1./3.) +1))

def function_n3(zeta):
    return np.sqrt(zeta/3)


fig, ax = plt.subplots()

zeta = np.linspace(0.0,1,20)

ax.plot(zeta, function_n2(zeta)/function_n2(1.))
ax.plot(zeta, function_n3(zeta)/function_n3(1.))


plt.show()
