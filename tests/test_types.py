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

# from dolfin.cpp.function import (
#     Constant as CppConstant,
#     Expression as CppExpression,
#     Function as CppFunction,
#     GenericFunction as CppGenericFunction,
#     MultiMeshFunction as CppMultiMeshFunction,
# )


MESH = UnitSquareMesh(4, 4)
DIFFUSION = Expression("1. + x[0] * x[0]", degree=2)
CONVECTION = Expression(("-x[1]", "x[0]"), degree=1)
REACTION = Expression("x[0] * x[0]", degree=2)


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


def test_handles_contants():
    eafe_assemble(MESH, Constant(1.0))
    # eafe_assemble(MESH, CppConstant(1.0))


def test_handles_expressions():
    eafe_assemble(MESH, Expression("diff", degree=0, diff=1.0))
    # eafe_assemble(MESH, CppExpression("diff", degree=0, diff=1.0))


def test_handles_functions():
    cg = FunctionSpace(MESH, "CG", 1)
    bdm = FunctionSpace(MESH, "BDM", 1)
    diffusion = interpolate(DIFFUSION, cg)
    convection = interpolate(CONVECTION, bdm)

    eafe_assemble(MESH, diffusion, convection)
    # dolfin.cpp.function.Function


# def test_handles_multi_mesh_functions():
#     # dolfin.MultiMeshFunction
#     # dolfin.cpp.function.MultiMeshFunction
#     # pass


# def test_handles_generic_functions():
#     # dolfin.cpp.function.GenericFunction
#     # pass
