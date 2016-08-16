import numpy as np

# Convention:
# u_{i - 1} = [0: -1]
# u_{i} = [1: -1]
# u_{i + 1} = [2:]


class Solver(object):
    def __init__(self, interval, nx, nt=0, dt=0):
        super(Solver, self).__init__()
        self.interval = interval
        a, b = interval
        self.dx = (b - a) / (nx - 1.)
        self.dt = dt
        self.nx = nx
        self.nt = nt

    def initial_conditions(self, argument=None):
        a, b = self.interval
        u = np.ones(self.nx)
        # setting u = 2 between 0.5 and 1 as per our I.C.s
        u[int(self.dx / 2): int(1 / self.dx + 1)] = 2
        return u

    def argument(self):
        a, b = self.interval
        return np.linspace(a, b, self.nx)

    def solve(self, function):
        u = self.initial_conditions(self.argument())
        for time in range(self.nt):
            u = function(u, self.dx, self.dt)
        return u

    def solve1d(self, function):
        u = []
        u.append(self.initial_conditions(self.argument()))
        for time in np.arange(0, self.nt, self.dt)[1:]:
            u.append(function(u[-1], self.dx, self.dt))
        return u


class MultiSolver(Solver):
    def __init__(self, solvers, nt=0, dt=0):
        super(MultiSolver, self).__init__((0, 0), 0, nt, dt)
        self.solvers = solvers
        self.dx = np.array([s.dx for s in self.solvers])
        self.nx = tuple(s.nx for s in self.solvers)

    def argument(self):
        return np.array([s.argument() for s in self.solvers])

    def initial_conditions(self, argument=None):
        return np.ones(self.nx)


class StationaryMultiSolver(MultiSolver):
    def __init__(self, solvers, default_tol=1e-3):
        super(StationaryMultiSolver, self).__init__(solvers, 0, 0)
        self.default_tol = default_tol

    def solve(self, function, tol=None):
        tol = tol if tol else self.default_tol
        quantity, error = self.initial_conditions(self.argument), 1.

        while error > tol:
            new = function(quantity, self.argument(), self.dx)
            error = self.distance(new, quantity)
            quantity = new
        return quantity

    def distance(self, new, old):
        numerator = np.sum(new - old)
        # avoid zeros in the denumerator ?
        denominator = np.sum(np.abs(old)) + self.default_tol
        return abs(numerator) / denominator
