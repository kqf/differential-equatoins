#!/usr/bin/python


from mpl_toolkits.mplot3d import Axes3D
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from Solver import Solver, MultiSolver


def profile(solver):
    u = np.ones(solver.nx)
    v = np.ones(solver.nx)
    dx, dy = solver.dx
    # set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    u[.5 / dy:1 / dy + 1, .5 / dx:1 / dx + 1] = 2
    # set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    v[.5 / dy:1 / dy + 1, .5 / dx:1 / dx + 1] = 2
    return u, v

# 2D convection


def convection(m, dr, dt):
    u, v = m
    vn = v.copy()
    un = u.copy()
    dx, dy = dr
    c = 1.
    dt = 0.2 * dx

    u[1:, 1:] = un[1:, 1:] - \
        (un[1:, 1:] * c * dt / dx * (un[1:, 1:] - un[1:, :-1])
         ) - vn[1:, 1:] * c * dt / dy * (un[1:, 1:] - un[:-1, 1:])
    v[1:, 1:] = vn[1:, 1:] - \
        (un[1:, 1:] * c * dt / dx * (vn[1:, 1:] - vn[1:, :-1])
         ) - vn[1:, 1:] * c * dt / dy * (vn[1:, 1:] - vn[:-1, 1:])

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

    return u, v


def main():
    msolver = MultiSolver(
        (Solver((0, 2.), 101), Solver((0, 2.), 101)), 81, 0.1)
    msolver.initial_conditions = lambda x = None: profile(msolver)

    uinit, vinit = msolver.initial_conditions(msolver.argument())
    u, v = msolver.solve(convection)
    # printu

    # Plot Initial Condition
    # the figsize parameter can be used to produce different sized images
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(*msolver.argument())

    res = u.copy()
    res[uinit > 1.] = uinit[uinit > 1.]
    ax.plot_surface(X, Y, res, cmap=cm.coolwarm)
    # surf = ax.plot_surface(X,Y, uinit, cmap=cm.coolwarm)
    # surf = ax.plot_surface(X,Y, 2. - v, cmap=cm.coolwarm);
    plt.show()


if __name__ == '__main__':
    main()
