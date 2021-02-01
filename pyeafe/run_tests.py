#! /usr/bin/env python
"""
Tests wrapping module
"""
from __future__ import division
from tests import (
    double_sin_3d_test,
    double_sin_reaction_test,
    double_sin_test,
    evaluation_test,
)


__all__ = ["run_tests"]


def run_tests():
    print("Begining the tests:")
    print("Testing over sin...")
    double_sin_test.run()
    print("Testing over sin reaction...")
    double_sin_reaction_test.run()
    print("Testing over sin 3D...")
    double_sin_3d_test.run()
    print("Testing variously defined coefficients...")
    evaluation_test.run()
    print("Done.")


if __name__ == "__main__":
    run_tests()
