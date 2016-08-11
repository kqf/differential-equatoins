

import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
import time as ttime

# Convention:
# u_{i - 1} = [0: -1]
# u_{i} = [1: -1]
# u_{i + 1} = [2:]

class Solver(object):
	def __init__(self, interval, nx, nt = 0, dt = 0):
		super(Solver, self).__init__()
		self.interval = interval
		a , b = interval
		self.dx	= (b - a) / (nx - 1.) 
		self.dt = dt
		self.nx = nx
		self.nt = nt

	def initial_conditions(self, argument = None):
		a, b = self.interval
		u = np.ones(self.nx)
		u[0.5/self.dx : 1/self.dx+1] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s
		return u

	def argument(self):
		a, b = self.interval
		return np.linspace(a, b, self.nx)

	def solve(self, function):
		u = self.initial_conditions( self.argument() )
		for time in range(self.nt):
			u = function(u, self.dx, self.dt)
		return u

	def solve1d(self, function):
		u = []
		u.append( self.initial_conditions( self.argument() ) )
		for time in np.arange(0, self.nt, self.dt)[1:]:
			u.append( function(u[-1], self.dx, self.dt) )
		return u

class MultiSolver(Solver):
	def __init__(self, solvers, nt = 0, dt = 0):
		super(MultiSolver, self).__init__((0, 0), 0, nt, dt)
		self.solvers = solvers 
		self.dx = np.array([s.dx for s in self.solvers])
		self.nx = tuple(s.nx for s in self.solvers)

	def argument(self):
		return np.array([s.argument() for s in self.solvers])

	def initial_conditions(self, argument = None):
		return np.ones( self.nx )
		