
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from operator import itemgetter
import pandas as pd
import yaml

from scipy.optimize import leastsq
import mush


def logistic4(x, A, B, C, D):
    """4PL lgoistic equation."""
    # return ((A-D)/(1.0+(C*np.exp(B*x)))) + D
    return ((A-D)/(1.0+((x/C)**B))) + D
def residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D = p
    err = y-logistic4(x, A, B, C, D)
    return err
    # find 1st minimum of function
def peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C,D = p
    return logistic4(x, A, B, C, D)

folder = "/home/marine/ownCloud/Research/Projets/mush/compaction"
file="/exp_1.00_coeff_0.05_radius_1.00/output_20.00.timestep"

data = pd.read_csv(folder+file, sep=" ")
phi = data["porosity"].values[:-1]

R = data["radius"].values[:-1]

min_phi = np.argmin(phi)
print(min_phi, phi[min_phi])

y_values = phi[min_phi:]
radius = R[min_phi:]

p0 = [phi[min_phi], np.abs((phi[-1]-phi[min_phi])/(radius[-1]-radius[0])), (radius[-1]+radius[0])/2 , phi[-1]]
plsq = leastsq(residuals, p0, args=(y_values, radius))
    # get the thickness of the sigmoid
A, B, C, D = plsq[0]
print(A, B, C, D)

fig, ax = plt.subplots()

ax.plot(R, phi,  '.g')
ax.plot(R, logistic4(R, A, B, C, D))
phi_C = logistic4(C, A, B, C, D)
ax.plot(R, phi_C +(R-C)*B*(D-A)/4/C, 'r' )
print(logistic4(C, A, B, C, D))
#ax.set_ylim([0, 1])

print(radius[-1]-C+logistic4(C, A, B, C, D)*4*C/B/(D-A))

plt.show()