from solver.solver import Solver

EXPECTED = [0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1]


def test_exponent():
    solver = Solver((0, 1), nx=2, nt=1, dt=0.1)
    solver.initial_conditions = lambda x: 0.1
    assert solver.solve1d(lambda u, dx, dt: u + dx) == EXPECTED
