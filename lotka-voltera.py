#!/usr/bin/python


import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np

from Solver import Solver

## convection equation in terms of finite differences
def convect_step(uv, dx, dt):
	un, vn = uv
	a, b = 1. , 1./4.
	c, d = 0.2, 0.6
	u = un + dt * ( a * un  - b * vn * un)
	v = vn + dt * ( c * vn * un  - d * vn)

	return u, v


def main():
	solver = Solver((0, 2.), 0, 10, 0.0001)
	solver.initial_conditions = lambda x = 1 : (10, 6)

	time = np.arange(0, 10, 0.0001)
	ans = solver.solve1d(convect_step)
	x, y =  zip(*ans)
	# print len(time), len(ans)
	plt.plot(x, y)
	# plt.plot(time, ans)

	# plt.plot(solver.argument(), solver.initial_conditions())
	# plt.plot(solver.argument(), solver.solve1(convect_step))
	plt.show()

if __name__ == '__main__':
	main()

