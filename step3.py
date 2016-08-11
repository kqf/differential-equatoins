#!/usr/bin/python


import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np

from Solver import Solver

## diffusion equation in terms of finite differences
def diffusion_eqn(un, dx, dt):
	c, sigma, nu = 1, 0.2, 0.3
	dt = sigma * dx ** 2 / nu
	res = un.copy()
	res[1:-1] =  un[1:-1] + nu * dt / dx **2 * (un[2:] - 2 * un[1:-1] + un[:-2])
	res[0] = 2
	res[-1] = 2
	return res
	# return  un[i] + nu * dt / dx **2 * (un[i + 1] - 2 * un[i] + un[i - 1])

def main():
	for xiterations in [41 + i  * 10 for i in range(10)]:
		solver = Solver((0, 2.), xiterations, 125, 0.025, iter_offset=(1, -1))
		if xiterations == 41:
			plt.plot(solver.argument(), solver.initial_conditions())
		plt.plot(solver.argument(), solver.solve(diffusion_eqn))
	plt.show()

if __name__ == '__main__':
	main()

