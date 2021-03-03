from inspect import getfullargspec
from typing import Callable, Union

import numpy as np
from dolfin import (
    Constant,
    Expression,
    Function,
    # MultiMeshFunction,
    Point,
    Cell,
)

# from dolfin.cpp.function import (
#     Constant as CppConstant,
#     Expression as CppExpression,
#     Function as CppFunction,
#     GenericFunction as CppGenericFunction,
#     MultiMeshFunction as CppMultiMeshFunction,
# )

Coefficient = Union[
    Constant,
    Expression,
    Function,
    # MultiMeshFunction,
    # CppConstant,
    # CppExpression,
    # CppFunction,
    # CppGenericFunction,
    # CppMultiMeshFunction,
]


def validate_coefficient(coefficient: Coefficient) -> bool:
    if issubclass(coefficient.__class__, Constant):
        return True

    if issubclass(coefficient.__class__, Expression):
        return True

    if issubclass(coefficient.__class__, Function):
        return True

    # if issubclass(coefficient.__class__, MultiMeshFunction):
    #     return True
    #
    # if issubclass(coefficient.__class__, CppConstant):
    #     return True
    #
    # if issubclass(coefficient.__class__, CppExpression):
    #     return True
    #
    # if issubclass(coefficient.__class__, CppFunction):
    #     return True
    #
    # if issubclass(coefficient.__class__, CppGenericFunction):
    #     return True
    #
    # if issubclass(coefficient.__class__, CppMultiMeshFunction):
    #     return True

    return False


# TODO: pass in mesh and define `eval_cell`
def create_safe_eval(
    expression: Coefficient,
    value_shape: int,
) -> Callable[[Point, Cell], Union[float, np.array]]:
    """
    Wrap expression to ensure consistent calling for various input coefficients.
    :param expression: pyeafe.Coefficient to be evaluated
    :param value_shape: integer output dimension of function

    :return: Return a function that has a consistent parameter list:
        (point: dolfin.Point, cell: dolfin.Cell) -> Union[np.array, float]
    """

    if not callable(expression):
        return lambda p, c: expression

    if hasattr(expression, "eval_cell") and callable(getattr(expression, "eval_cell")):

        def evaluate(point, cell):
            values = np.empty(value_shape, dtype=np.float_)
            expression.eval_cell(values, point, cell)
            return values

        return evaluate

    arg_count: int = len(getfullargspec(expression).args)
    return lambda point, cell: expression(point) if arg_count == 1 else expression
