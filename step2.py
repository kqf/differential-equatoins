#!/usr/bin/python


import matplotlib.pyplot as plt
from Solver import Solver


# nonlinear convection equation in terms of finite differences


def nlin_convection_step(un, dx, dt):
    result = un.copy()
    # forward difference
    result[1:] = un[1:] - un[1:] * dt / dx * (un[1:] - un[0:-1])

    # Boundary condition for forward difference
    result[0] = 1
    return result


def main():
    for xiterations in [41 + i for i in range(4)]:
        solver = Solver((0, 2.), xiterations, 25, 0.025)
        if xiterations == 41:
            plt.plot(solver.argument(), solver.initial_conditions())
        plt.plot(solver.argument(), solver.solve(nlin_convection_step))

    # plt.yscale('log')
    plt.show()


if __name__ == '__main__':
    main()
