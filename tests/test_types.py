from pyeafe import eafe_assemble
import pytest

from dolfin import (
    Constant,
    Expression,
    FunctionSpace,
    # MultiMeshFunction,
    UnitSquareMesh,
    interpolate,
)


MESH = UnitSquareMesh(4, 4)
CG = FunctionSpace(MESH, "CG", 1)
BDM = FunctionSpace(MESH, "BDM", 1)
DIFFUSION = Expression("1. + x[0] * x[0]", degree=2)
CONVECTION = Expression(("-x[1]", "x[0]"), degree=1)
REACTION = Expression("x[0] * x[0]", degree=2)

VALID_DIFFUSIONS = [
    Constant(1.0),
    DIFFUSION,
    interpolate(DIFFUSION, CG),
]

VALID_CONVECTIONS = [
    None,
    Constant((1.0, -1.0)),
    CONVECTION,
    interpolate(CONVECTION, BDM),
]

VALID_REACTIONS = [
    None,
    Constant(1.0),
    REACTION,
    interpolate(REACTION, CG),
]


def test_assemble_with_mixed_types():
    for diffusion in VALID_DIFFUSIONS:
        for convection in VALID_CONVECTIONS:
            for reaction in VALID_REACTIONS:
                eafe_assemble(MESH, diffusion, convection, reaction)


def test_assemble_raises_for_invalid_mesh():
    with pytest.raises(TypeError):
        eafe_assemble(None, DIFFUSION)


def test_assemble_raises_for_invalid_diffusion():
    with pytest.raises(TypeError):
        eafe_assemble(MESH, None)

    with pytest.raises(TypeError):
        eafe_assemble(MESH, 1.0)


def test_assemble_raises_for_invalid_convection():
    with pytest.raises(TypeError):
        eafe_assemble(MESH, DIFFUSION, 5.0)

    with pytest.raises(ValueError):
        eafe_assemble(MESH, DIFFUSION, Expression("1.", degree=1))

    with pytest.raises(ValueError):
        eafe_assemble(MESH, DIFFUSION, Expression(("1."), degree=1))


def test_assemble_raises_for_invalid_reaction():
    with pytest.raises(TypeError):
        eafe_assemble(MESH, DIFFUSION, CONVECTION, 5.0)

    with pytest.raises(TypeError):
        eafe_assemble(MESH, DIFFUSION, reaction=5.0)
