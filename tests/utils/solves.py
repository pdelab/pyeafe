from typing import Optional
from dolfin import (
    assemble,
    dx,
    errornorm,
    DirichletBC,
    Function,
    FunctionSpace,
    LUSolver,
    Mesh,
    TestFunction,
)

from pyeafe import eafe_assemble, Coefficient


def assert_solves(
    mesh: Mesh,
    diffusion: Coefficient,
    convection: Optional[Coefficient],
    reaction: Optional[Coefficient],
    source: Coefficient,
    exact: Coefficient,
    l2_tol: Optional[float] = 1.0e-8,
    h1_tol: Optional[float] = 1.0e-6,
):
    eafe_matrix = eafe_assemble(mesh, diffusion, convection, reaction)

    pw_linears = FunctionSpace(mesh, "Lagrange", 1)
    test_function = TestFunction(pw_linears)
    rhs_vector = assemble(source * test_function * dx)

    bc = DirichletBC(pw_linears, exact, lambda _, on_bdnry: on_bdnry)
    bc.apply(eafe_matrix, rhs_vector)

    solution = Function(pw_linears)
    solver = LUSolver(eafe_matrix, "default")
    solver.parameters["symmetric"] = False
    solver.solve(solution.vector(), rhs_vector)

    l2_err: float = errornorm(exact, solution, "l2", 3)
    assert l2_err <= l2_tol, f"L2 error too large: {l2_err} > {l2_tol}"

    h1_err: float = errornorm(exact, solution, "H1", 3)
    assert h1_err <= h1_tol, f"H1 error too large: {h1_err} > {h1_tol}"
