#! /usr/bin/env python
from dolfin import (
    DirichletBC,
    Expression,
    Function,
    FunctionSpace,
    LUSolver,
    TestFunction,
    UnitCubeMesh,
    assemble,
    dx,
    errornorm,
)

import pyeafe


def test_double_sin_3d_problem():
    def boundary(x, on_boundary):
        return on_boundary

    diffusivity = 1.0e-2
    diffusion = Expression("diffusivity", degree=0, diffusivity=diffusivity)
    convection = Expression(("x[1]", "-x[0]", "0."), degree=1)
    exact_solution = Expression(
        "sin(2 * DOLFIN_PI * x[0]) * cos(2 * DOLFIN_PI * x[1])", degree=4
    )

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
    mesh = UnitCubeMesh(granularity, granularity, granularity)

    continuous_pw_linear_space = FunctionSpace(mesh, "Lagrange", 1)
    test_function = TestFunction(continuous_pw_linear_space)
    linear_functional = right_side_expression * test_function * dx

    stiffness_matrix = pyeafe.eafe_assemble(mesh, diffusion, convection)
    rhs_vector = assemble(linear_functional)

    bc = DirichletBC(continuous_pw_linear_space, exact_solution, boundary)
    bc.apply(stiffness_matrix, rhs_vector)

    solution = Function(continuous_pw_linear_space)
    solver = LUSolver(stiffness_matrix, "default")
    solver.parameters["symmetric"] = False
    solver.solve(solution.vector(), rhs_vector)

    assert (
        errornorm(exact_solution, solution, "l2", 3) < 2.12e-1
    ), "L2 error increased! Solver failed"

    assert (
        errornorm(exact_solution, solution, "H1", 3) < 2.33e-0
    ), "H1 error increased! Solver failed"
