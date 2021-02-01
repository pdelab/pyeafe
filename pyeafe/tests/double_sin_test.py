#! /usr/bin/env python
from dolfin import (
    DirichletBC,
    Expression,
    Function,
    FunctionSpace,
    LUSolver,
    TestFunction,
    UnitSquareMesh,
    assemble,
    dx,
    errornorm,
)

import numpy as np
import pyeafe


def run():
    def boundary(x, on_boundary):
        return on_boundary

    diffusivity = 1.0e-2

    def diffusion_expression(x):
        return diffusivity

    def convection_expression(x):
        return np.array([x[1], -x[0]])

    right_side_expression = Expression(
        "8 * DOLFIN_PI * DOLFIN_PI * diffusivity\
        * sin(2 * DOLFIN_PI * x[0]) * cos(2 * DOLFIN_PI * x[1])\
        - 2 * DOLFIN_PI * x[1] * cos(2 * DOLFIN_PI * x[0])\
        * cos(2 * DOLFIN_PI * x[1])\
        - 2 * DOLFIN_PI * x[0] * sin(2 * DOLFIN_PI * x[0])\
        * sin(2 * DOLFIN_PI * x[1])",
        degree=2,
        diffusivity=diffusivity,
    )

    granularity = 8
    mesh = UnitSquareMesh(granularity, granularity)

    continuous_pw_linear_space = FunctionSpace(mesh, "Lagrange", 1)
    test_function = TestFunction(continuous_pw_linear_space)
    linear_functional = right_side_expression * test_function * dx

    stiffness_matrix = pyeafe.eafe_assemble(
        mesh, diffusion_expression, convection_expression
    )
    rhs_vector = assemble(linear_functional)

    exact_solution = Expression(
        "sin(2 * DOLFIN_PI * x[0]) * cos(2 * DOLFIN_PI * x[1])", degree=4
    )
    bc = DirichletBC(continuous_pw_linear_space, exact_solution, boundary)
    bc.apply(stiffness_matrix, rhs_vector)

    solution = Function(continuous_pw_linear_space)
    solver = LUSolver(stiffness_matrix, "default")
    solver.parameters["symmetric"] = False
    solver.solve(solution.vector(), rhs_vector)

    l2_error = errornorm(exact_solution, solution, "l2", 3)
    if l2_error > 2.12e-1:
        raise ValueError("L2 error increased! Solver failed")

    h1_error = errornorm(exact_solution, solution, "H1", 3)
    if h1_error > 2.33e-0:
        raise ValueError("H1 error increased! Solver failed")

    print(f"L2 error = {l2_error}")
    print(f"H1 error = {h1_error}")
    print("Success!")


if __name__ == "__main__":
    run()
