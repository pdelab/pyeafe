from dolfin import Expression
from inspect import getfullargspec
from typing import Union

import numpy as np

FlexibleExpression = Union[Expression, callable, np.array, float]


# TODO: pass in mesh and define `eval_cell`
def create_safe_eval(
    expression: FlexibleExpression,
    value_shape: int,
) -> callable:
    """
    Wrap expression to ensure consistent calling for various input coefficients.
    :param expression: dolfin.Expression, python function, or constant value
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
