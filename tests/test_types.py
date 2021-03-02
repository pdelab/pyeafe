from dolfin import UnitSquareMesh, Expression
from pyeafe import eafe_assemble
import pytest

MESH = UnitSquareMesh(4, 4)
DIFFUSION = Expression("1. + x[0] * x[0]", degree=2)
CONVECTION = Expression(("-x[1]", "x[0]"), degree=1)
DIFFUSION = Expression("x[0] * x[0]", degree=2)


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
