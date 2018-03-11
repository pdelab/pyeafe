from dolfin import *
import numpy as np
import inspect
import sys
import os

current_frame = inspect.getfile(inspect.currentframe())
current_dir = os.path.dirname(os.path.abspath(current_frame))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)
import assembly as eafe

def run():

    def boundary(x, on_boundary):
        return on_boundary

    diffusivity = 1.0E-2
    def diffusion_expression(x):
        return diffusivity

    def convection_expression(x):
        return np.array([x[1],-x[0], 0.0])

    exact_solution = Expression(
        "sin(2 * DOLFIN_PI * x[0]) * cos(2 * DOLFIN_PI * x[1])",
        degree=4
    )

    right_side_expression = Expression(
        "8 * DOLFIN_PI * DOLFIN_PI * diffusivity\
        * sin(2 * DOLFIN_PI * x[0]) * cos(2 * DOLFIN_PI * x[1])\
        - 2 * DOLFIN_PI * x[1] * cos(2 * DOLFIN_PI * x[0])\
        * cos(2 * DOLFIN_PI * x[1])\
        - 2 * DOLFIN_PI * x[0] * sin(2 * DOLFIN_PI * x[0])\
        * sin(2 * DOLFIN_PI * x[1])",
        degree=2,
        diffusivity=diffusivity
    )

    granularity = 8
    mesh = UnitCubeMesh(granularity, granularity, granularity)

    continuous_pw_linear_space = FunctionSpace(mesh, "Lagrange", 1)
    test_function = TestFunction(continuous_pw_linear_space)
    linear_functional = right_side_expression * test_function * dx

    stiffness_matrix = eafe.eafe_assemble(mesh, diffusion_expression, convection_expression)
    rhs_vector = assemble(linear_functional)

    bc = DirichletBC(continuous_pw_linear_space, exact_solution, boundary)
    bc.apply(stiffness_matrix, rhs_vector)

    solution = Function(continuous_pw_linear_space)
    solver = LUSolver(stiffness_matrix, "default")
    solver.parameters["symmetric"] = False
    solver.solve(solution.vector(), rhs_vector)

    if (errornorm(exact_solution, solution, 'l2', 3) > 2.12E-1):
        print "L2 error increased! Solver failed"
        exit(1)

    if (errornorm(exact_solution, solution, 'H1', 3) > 2.33E-0):
        print "H1 error increased! Solver failed"
        exit(1)

    print "L2 error = ", errornorm(exact_solution, solution, 'l2', 3)
    print "H1 error = ", errornorm(exact_solution, solution, 'H1', 3)
    print "Success!"

if __name__ == "__main__":
    run()
