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


def wrap_vector_constants(coefficient: Coefficient, value_shape: int):
    values = np.empty(value_shape, dtype=np.float_)
    coefficient.eval(values, np.zeros(value_shape))

    if values.flatten().shape == (value_shape,):
        return lambda point, cell: values

    raise ValueError("Invalid vector coefficient: dimension must match mesh dimension")


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

    if value_shape > 1:
        if issubclass(coefficient.__class__, Constant):
            return wrap_vector_constants(coefficient, value_shape)

        coefficient_rank: int = coefficient.value_rank()
        if coefficient_rank != 1:
            raise ValueError("Invalid vector coefficient: value_rank must be 1")

        if coefficient.value_dimension(0) != value_shape:
            raise ValueError(
                "Invalid vector coefficient: value_dimension(0) must match mesh spatial dimension"  # noqa E501
            )

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
