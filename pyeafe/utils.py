from inspect import getfullargspec
from typing import Callable, Union

import numpy as np
from dolfin import (
    Constant,
    Expression,
    Function,
    MultiMeshFunction,
    Point,
    Cell,
)

Coefficient = Union[Constant, Expression, Function, MultiMeshFunction]


def validate_coefficient(coefficient: Coefficient) -> bool:
    if issubclass(coefficient.__class__, Constant):
        return True

    if issubclass(coefficient.__class__, Expression):
        return True

    if issubclass(coefficient.__class__, Function):
        return True

    if issubclass(coefficient.__class__, MultiMeshFunction):
        return True

    return False


def ensure_cell_eval(
    coefficient: Coefficient,
    value_shape: int,
) -> Callable[[Point, Cell], Union[float, np.array]]:
    """
    Wrap coefficient to ensure consistent calling for various input coefficients.
    :param coefficient: pyeafe.Coefficient to be evaluated
    :param value_shape: integer output dimension of function

    :return: Return a function that has a uses eval_cell when available.
    """

    if not validate_coefficient(coefficient):
        raise TypeError("Invalid coefficient: Must inherit from pyeafe.Coefficient")

    if hasattr(coefficient, "eval_cell") and callable(
        getattr(coefficient, "eval_cell")
    ):

        def evaluate(point, cell):
            values = np.empty(value_shape, dtype=np.float_)
            coefficient.eval_cell(values, point, cell)
            return values

        return evaluate

    arg_count: int = len(getfullargspec(coefficient).args)
    return lambda point, cell: coefficient(point) if arg_count == 1 else coefficient
