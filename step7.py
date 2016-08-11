#!/usr/bin/python


from mpl_toolkits.mplot3d import Axes3D
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
import sympy
from matplotlib import cm

from Solver import Solver, MultiSolver


def profile(solver):
	u = np.ones(solver.nx)
	dx, dy = solver.dx
	u[.5/dy : 1/dy + 1, .5 / dx : 1 / dx + 1] = 2
	return u

# 2D linear convection
def convection(u, dr ,dt):
	un = u.copy()
	dx, dy  = dr
	c = 1.
	dt = 0.2 * dx 
	nu = 0.05
	
	u[1:-1,1:-1] = un[1:-1,1:-1]+nu*dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2]) + nu*dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])
	u[0,:] = 1
	u[-1,:] = 1
	u[:,0] = 1
	u[:,-1] = 1
	return u

def diffuse(ntsteps):
	msolver = MultiSolver( (Solver((0, 2.), 51), Solver((0, 2.), 51) ),  ntsteps,  0.1)
	msolver.initial_conditions = lambda x = None: profile(msolver)

	res  = msolver.solve(convection)

	###Plot Initial Condition
	fig = plt.figure(figsize=(11,7), dpi=100)          ##the figsize parameter can be used to produce different sized images
	ax  = fig.gca(projection='3d')                      
	X, Y = np.meshgrid(*msolver.argument())                            
	surf = ax.plot_surface(X,Y, res, cmap=cm.coolwarm)
	plt.show()

def main():
	for i in range(30): diffuse(100 + 5 * i) 


if __name__ == '__main__':
	main()

