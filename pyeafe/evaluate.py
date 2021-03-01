from inspect import getfullargspec
import numpy as np


def returnZeroOutput(dimension):
    if dimension < 1:
        print("value_shape must be positive")
        raise ValueError("Cannot safely evaluate required function")

    if dimension == 1:

        def returnZero(point, cell):
            return 0.0

        return returnZero

    def returnZeros(point, cell):
        return np.zeros(dimension)

    return returnZeros


"""
Return a function that has a consistent parameter list
as used in eafe assembly routines to safely evaluate
a previously-defined DOLFIN function or expression
representing a PDE coefficient

Usage:
    create_safe_eval(expression, value_shape)

    - expression: python function, DOLFIN function,
        or DOLFIN expression to be evaluated.if not provided,
        value_shape is used to return a zero scalar value or zero vector
    - value_shape: integer output dimension of function.
        Set to return scalar values by default (value_shape=1).
"""


def create_safe_eval(expression, value_shape):
    if hasattr(expression, "eval_cell") and callable(getattr(expression, "eval_cell")):

        def evaluate(point, cell):
            values = np.empty(value_shape, dtype=np.float_)
            expression.eval_cell(values, point, cell)
            return values

        return evaluate

    if callable(expression):
        signature = getfullargspec(expression)
        if len(signature.args) == 1:

            def eval_with_cell_stub(point, cell):
                return expression(point)

            return eval_with_cell_stub

        return expression

    def return_value(point, cell):
        return expression

    return return_value
