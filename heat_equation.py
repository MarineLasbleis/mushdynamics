"""
	This program solves the heat equation
		u_t = u_xx
	with dirichlet boundary condition
		u(0,t) = u(1,t) = 0
	with the Initial Conditions
		u(x,0) = 10*sin( pi*x )
	over the domain x = [0, 1]
 
	The program solves the heat equation using a finite difference method
	where we use a center difference method in space and Crank-Nicolson 
	in time.
"""
 
import scipy as sc
import scipy.sparse as sparse
import scipy.sparse.linalg
import numpy as np
import pylab as pl
#import include.CreateMovie as movie
import matplotlib.pyplot as plt
 
# Number of internal points
N = 10
 
# Calculate Spatial Step-Size
h = 1/(N+1.0)
 
# Create Temporal Step-Size, TFinal, Number of Time-Steps
k = h/2
TFinal = 1
NumOfTimeSteps = 1#int(TFinal/k)
 
# Create grid-points on x axis
x = np.linspace(0,1,N+2)
x = x[1:-1]
 
# Initial Conditions
u = np.transpose(np.mat(10*np.sin(np.pi*x)))
 
# Second-Derivative Matrix
data = np.ones((3, N))
data[1] = -2*data[1]
diags = [-1,0,1]
D2 = sparse.spdiags(data,diags,N,N)/(h**2)
 
# Identity Matrix
I = sparse.identity(N)
 
# Data for each time-step
data = []
 
for i in range(NumOfTimeSteps):
	# Solve the System: (I - k/2*D2) u_new = (I + k/2*D2)*u_old
	A = (I -k/2*D2)
	b = ( I + k/2*D2 )*u
	u = np.transpose(np.mat(sparse.linalg.spsolve( A,  b ) ))

	print(A)
	print(np.mat(A))
	print(b.shape)
	print(u.shape)
 
	# Save Data
	data.append(u)
 


A = sparse.diags([[1,1],[1,1,1],[1,1]],[-1,0,1])
b = np.array([1,1,1])
print(np.transpose(np.mat(sparse.linalg.spsolve( A,  b ) )))

 
# Define the Frame Speed and Movie Length
FPS = 20
MovieLength = 10
 
# Function to plot any given Frame
def plotFunction( frame ):
	plt.plot(x, data[int(NumOfTimeSteps*frame/(FPS*MovieLength))] )
	plt.axis((0,1,0,10.1))
 
# Generate the movie
#movie.CreateMovie(plotFunction, int(MovieLength*FPS), FPS)