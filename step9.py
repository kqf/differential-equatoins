#!/usr/bin/python


from mpl_toolkits.mplot3d import Axes3D
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
import sympy
from matplotlib import cm

from Solver import Solver, MultiSolver, StationaryMultiSolver


def profile(solver):
	u = np.zeros(solver.nx)
	return u

# 2D  laplace
def laplace(pn, r, dr):
    x, y = r
    dx, dy = dr
    p = pn.copy()
    p[1:-1,1:-1] = (dy**2*(pn[1:-1,2:]+pn[1:-1,0:-2]) + dx**2*(pn[2:,1:-1]+pn[0:-2,1:-1]))/(2*(dx**2+dy**2))

    p[:,0] = 0           ##p = 0 @ x = 0
    p[:,-1] = y          ##p = y @ x = 2
    p[0,:] = p[1,:]      ##dp/dy = 0 @ y = 0
    p[-1,:] = p[-2,:]    ##dp/dy = 0 @ y = 1

    return p

def plot2D(x, y, p):
    fig = plt.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    # ax.set_xlim(0,2)
    # ax.set_ylim(0,1)
    ax.view_init(30,225)

def laplace_solve():
	msolver = StationaryMultiSolver( (Solver((0, 2.), 31), Solver((0, 1.), 31) ) )
	msolver.initial_conditions = lambda x = None: profile(msolver)

	res  = msolver.solve(laplace)
	x, y = msolver.argument()
	plot2D(x, y, res)
	plt.show()

def main():
	laplace_solve()


if __name__ == '__main__':
	main()

