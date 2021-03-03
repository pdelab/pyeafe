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
    interpolate,
)

import pyeafe


def boundary(x, on_boundary):
    return on_boundary


true_solution = Expression(
    "sin(2 * DOLFIN_PI * x[0]) * cos(2 * DOLFIN_PI * x[1])", degree=4
)
diffusivity = 1.0e-2
rhs = Expression(
    "8 * DOLFIN_PI * DOLFIN_PI * diffusivity\
    * sin(2 * DOLFIN_PI * x[0]) * cos(2 * DOLFIN_PI * x[1])\
    - 2 * DOLFIN_PI * x[1] * cos(2 * DOLFIN_PI * x[0])\
    * cos(2 * DOLFIN_PI * x[1])\
    - 2 * DOLFIN_PI * x[0] * sin(2 * DOLFIN_PI * x[0])\
    * sin(2 * DOLFIN_PI * x[1])",
    degree=2,
    diffusivity=diffusivity,
)


def compute_error(
    stiffness_matrix,
    mesh,
    right_side_expression=rhs,
    exact_solution=true_solution,
    **kwargs
):
    continuous_pw_linear_space = FunctionSpace(mesh, "CG", 1)
    test_function = TestFunction(continuous_pw_linear_space)
    linear_functional = right_side_expression * test_function * dx
    rhs_vector = assemble(linear_functional)

    bc = DirichletBC(continuous_pw_linear_space, exact_solution, boundary)
    bc.apply(stiffness_matrix, rhs_vector)

    computed_solution = Function(continuous_pw_linear_space)
    solver = LUSolver(stiffness_matrix, "default")
    solver.parameters["symmetric"] = False
    solver.solve(computed_solution.vector(), rhs_vector)

    l2_error = errornorm(exact_solution, computed_solution, "l2", 3)
    if l2_error > 2.12e-1:
        raise ValueError("L2-error is larger than expected: " + l2_error)

    h1_error = errornorm(exact_solution, computed_solution, "H1", 3)
    if h1_error > 2.33e-0:
        raise ValueError("H1-error is larger than expected: " + h1_error)


def test_evaluate():
    granularity = 8
    mesh = UnitSquareMesh(granularity, granularity)
    diffusion = Expression("1.e-2", degree=0)
    convection = Expression(("x[1]", "-x[0]"), degree=1)
    eafe = pyeafe.eafe_assemble(mesh, diffusion, convection)
    compute_error(eafe, mesh)

    diffusion_expression = Expression("diffusivity", diffusivity=diffusivity, degree=2)
    convection_expression = Expression(("x[1]", "-x[0]"), degree=2)
    eafe_expression = pyeafe.eafe_assemble(
        mesh, diffusion_expression, convection_expression
    )
    compute_error(eafe_expression, mesh)

    cg = FunctionSpace(mesh, "CG", 1)
    bdm = FunctionSpace(mesh, "BDM", 1)
    diffusion_interpolant = interpolate(diffusion_expression, cg)
    convection_interpolant = interpolate(convection_expression, bdm)
    eafe_function = pyeafe.eafe_assemble(
        mesh, diffusion_interpolant, convection_interpolant
    )
    compute_error(eafe_function, mesh)
