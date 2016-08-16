#!/usr/bin/python


import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np

from solver.solver import Solver

## convection equation in terms of finite differences
def convect_step(un, dx, dt):
	result = un.copy()
	c = 1
	# forward difference 
	result[1:] = un[1:] - c * dt / dx * (un[1:] - un[0:-1])

	#Boundary condition	for forward difference
	result[0] = 1  
	return result
def main():
	for xiterations in [41 + 10 * i for i in range(3)]:
		solver = Solver((0, 2.), xiterations, 25, 0.025)
		if xiterations == 41:
			plt.plot(solver.argument(), solver.initial_conditions())
		plt.plot(solver.argument(), solver.solve(convect_step))
	plt.show()

if __name__ == '__main__':
	main()

