import numpy as np
from dolfin import (
    interpolate,
    Constant,
    Expression,
    FunctionSpace,
    UnitSquareMesh,
)

from tests.utils.solves import assert_solves
from tests.utils.types import assemble_with_mixed_types

MESH = UnitSquareMesh(16, 16)
EXACT_EXPRESSION = "sin(DOLFIN_PI * x[0]) * sin(DOLFIN_PI * x[1])"
EXACT = Expression(EXACT_EXPRESSION, degree=4)

DIFFUSIVITY = 1.0 / (2.0 * np.pi * np.pi)
DIFFUSION_TERM = f"2. * DOLFIN_PI * DOLFIN_PI * {DIFFUSIVITY} * {EXACT_EXPRESSION}"

CONVECTIVITY = [1.0 / np.pi, 1.0 / np.pi]
CONVECTION_TERM = "+".join(
    [
        "-cos(DOLFIN_PI * x[0]) * sin(DOLFIN_PI * x[1])",
        "-sin(DOLFIN_PI * x[0]) * cos(DOLFIN_PI * x[1])",
    ]
)

REACTIVITY = 1.0
REACTION_TERM = EXACT_EXPRESSION


def test_2d_types():
    mesh = UnitSquareMesh(2, 2)
    CG = FunctionSpace(mesh, "CG", 1)
    BDM = FunctionSpace(mesh, "BDM", 1)
    diffusion = Expression("1. + x[0] * (1. - x[0])", degree=2)
    convection = Expression(("-x[1]", "x[0]"), degree=1)
    reaction = Expression("x[0] * (1. - x[0])", degree=2)

    assemble_with_mixed_types(
        mesh,
        [Constant(1.0), diffusion, interpolate(diffusion, CG)],
        [None, Constant((1.0, -1.0)), convection, interpolate(convection, BDM)],
        [None, Constant(1.0), reaction, interpolate(reaction, CG)],
    )


def test_2d_diffusion():
    assert_solves(
        MESH,
        Constant(DIFFUSIVITY),
        None,
        None,
        Expression(DIFFUSION_TERM, degree=4),
        EXACT,
        5.4e-3,
        2.2e-1,
    )


def test_2d_convection():
    assert_solves(
        MESH,
        Constant(DIFFUSIVITY),
        Constant(CONVECTIVITY),
        None,
        Expression(f"{DIFFUSION_TERM} + {CONVECTION_TERM}", degree=4),
        EXACT,
        1.0e-2,
        2.2e-1,
    )


def test_2d_reaction():
    assert_solves(
        MESH,
        Constant(DIFFUSIVITY),
        None,
        Constant(REACTIVITY),
        Expression(f"{DIFFUSION_TERM} + {REACTION_TERM}", degree=4),
        EXACT,
        6.2e-3,
        2.2e-1,
    )


def test_2d_convection_reaction():
    assert_solves(
        MESH,
        Constant(DIFFUSIVITY),
        Constant(CONVECTIVITY),
        Constant(REACTIVITY),
        Expression(
            f"{DIFFUSION_TERM} + {CONVECTION_TERM} + {REACTION_TERM}",
            degree=4,
        ),
        EXACT,
        8.5e-3,
        2.2e-1,
    )
