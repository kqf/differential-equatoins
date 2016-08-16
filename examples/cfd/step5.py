#!/usr/bin/python


import matplotlib.pyplot as plt
import numpy as np

from solver.solver import Solver, MultiSolver


def profile(solver):
    u = np.ones(solver.nx)
    dx, dy = solver.dx
    u[.5 / dy: 1 / dy + 1, .5 / dx: 1 / dx + 1] = 2
    return u


def convection(u, dr, dt):
    # 2D linear convection
    un = u.copy()
    dx, dy = dr
    c = 1.
    dt = 0.2 * dx

    u[1:, 1:] = un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) - \
        (c * dt / dy * (un[1:, 1:] - un[:-1, 1:]))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    return u


def main():
    msolver = MultiSolver((Solver((0, 2.), 81), Solver((0, 2.), 81)), 100, 0.1)
    msolver.initial_conditions = lambda x = None: profile(msolver)

    init = msolver.initial_conditions(msolver.argument())
    res = msolver.solve(convection)

    # Plot Initial Condition
    # the figsize parameter can be used to produce different sized images
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(*msolver.argument())
    ax.plot_surface(X, Y, init)
    ax.plot_surface(X, Y, res)
    plt.show()


if __name__ == '__main__':
    main()
