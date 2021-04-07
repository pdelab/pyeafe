import numpy as np
from dolfin import (
    interpolate,
    Constant,
    Expression,
    FunctionSpace,
    UnitIntervalMesh,
)

from tests.utils.solves import assert_solves
from tests.utils.types import assemble_with_mixed_types

MESH = UnitIntervalMesh(32)
EXACT_EXPRESSION = "sin(DOLFIN_PI * x[0])"
EXACT = Expression(EXACT_EXPRESSION, degree=4)

DIFFUSIVITY = 1.0 / (np.pi * np.pi)
DIFFUSION_TERM = f"DOLFIN_PI * DOLFIN_PI * {DIFFUSIVITY} * {EXACT_EXPRESSION}"

CONVECTIVITY = [1.0 / np.pi]
CONVECTION_TERM = "-cos(DOLFIN_PI * x[0])"

REACTIVITY = 1.0
REACTION_TERM = EXACT_EXPRESSION


def test_1d_types():
    mesh = UnitIntervalMesh(4)
    CG = FunctionSpace(mesh, "CG", 1)
    DG = FunctionSpace(mesh, "DG", 0)
    diffusion = Expression("1. + x[0] * (1. - x[0])", degree=2)
    convection = Expression("-x[0]", degree=1)
    dg_convection = interpolate(convection, DG)
    reaction = Expression("x[0] * (1. - x[0])", degree=2)

    assemble_with_mixed_types(
        mesh,
        [Constant(1.0), diffusion, interpolate(diffusion, CG)],
        [None, Constant(1.0), convection, dg_convection],
        [None, Constant(1.0), reaction, interpolate(reaction, CG)],
    )


def test_1d_diffusion():
    assert_solves(
        MESH,
        Constant(DIFFUSIVITY),
        None,
        None,
        Expression(DIFFUSION_TERM, degree=4),
        EXACT,
        6.3e-4,
        6.3e-2,
    )


def test_1d_convection():
    assert_solves(
        MESH,
        Constant(DIFFUSIVITY),
        Constant(CONVECTIVITY),
        None,
        Expression(f"{DIFFUSION_TERM} + {CONVECTION_TERM}", degree=4),
        EXACT,
        1.1e-3,
        6.3e-2,
    )


def test_1d_reaction():
    assert_solves(
        MESH,
        Constant(DIFFUSIVITY),
        None,
        Constant(REACTIVITY),
        Expression(f"{DIFFUSION_TERM} + {REACTION_TERM}", degree=4),
        EXACT,
        8.9e-4,
        6.3e-2,
    )


def test_1d_convection_reaction():
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
        1.2e-3,
        6.3e-2,
    )
