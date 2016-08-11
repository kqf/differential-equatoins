#!/usr/bin/python


import ROOT
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
import sympy

from Solver import Solver

def profile():
	x, nu, t = sympy.symbols('x nu t')
	phi = sympy.exp(-(x-4*t)**2/(4*nu*(t+1))) + sympy.exp(-(x-4*t-2*np.pi)**2/(4*nu*(t+1)))
	phiprime = phi.diff(x)
	u = -2 * nu * (phiprime / phi) + 4
	return sympy.utilities.lambdify((t, x, nu), u) 

## burgers equation in terms of finite differences
def burgers_eqn(un, dx, dt):
	nu = 0.07
	dt =  dx * nu
	res = un.copy()

	fff = lambda l, c, r: c - c * dt/dx * (c - l) + nu * (r - 2 * c + l) * dt / dx**2

	a, b, c = un[:-2], un[1:-1], un[2:]
	res[1: -1] = fff(a, b, c)
	res[0]  = fff(un[-1], un[0], un[1])
	res[-1] = fff(un[-2], un[-1], un[0])
	return res

def main():
	for tt in range(100):
		t = 100 + 10 * tt
		solver = Solver((0, 2. * np.pi), 101, t, 0.025, (1, -1), True)
		nu = 0.07
		f = profile()
		initial = np.vectorize(lambda x: f(0, x, nu))
		final = np.vectorize(lambda x: f(t * solver.dx * nu, x, nu))

		solver.initial_conditions = lambda :  initial( solver.argument() )

		plt.figure(figsize=(11,7), dpi=100)

		x = solver.argument()
		plt.plot(x, solver.initial_conditions(), marker='o', lw=2, label='Inital state')
		plt.plot(x, solver.solve(burgers_eqn), marker='o', lw=2, label='Computational')
		plt.plot(x, final(x), marker='o', lw=2, label='Analytical final')

		plt.xlim([0,2*np.pi])
		plt.ylim([0,10])
		plt.legend()
		plt.show()

if __name__ == '__main__':
	main()

