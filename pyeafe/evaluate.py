'''
Safely evaluate a function based by providing default
arguments and a consistent parameter list as used in
eafe assembly routines


Usage:
    create_safe_eval(fn, output_dim, strict)

    - fn: [optional] function to be evaluated,
        if not provided, output_dim is used to return
        a zero scalr value or zero vector
    - output_dim: [optional] integer output dimension of function.
        set to return scalar values by default (output_dim=1).
    - strict: boolean value that throws an error if true and fn==None.
        If this parameter is false, fn can default to zero values
'''

from inspect import getargspec
import numpy as np

__all__ = ['create_safe_eval']


def create_safe_eval(fn=None, output_dim=1, strict=False, **kwargs):
    if (output_dim < 1):
        print("output_dim must be positive")
        raise ValueError("Cannot safely evaluate required function")

    if (fn is None and strict is True):
        print("fn must be provided in strict mode")
        raise ValueError("Cannot safely evaluate required function")

    elif (fn is None and output_dim == 1):
        def returnZero(point, cell):
            return 0.0

        return returnZero

    elif (fn is None and output_dim > 1):
        def returnZeros(point, cell):
            return np.zeros(spatial_dim)

        return returnZeros

    signature = getargspec(fn)
    if (len(signature.args) == 1):
        def safe_fn(point, cell):
            return fn(point)

        return safe_fn

    return fn
