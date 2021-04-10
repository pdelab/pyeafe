# PyEAFE

This package is an implementation of the edge-averaged finite element (EAFE) approximation
to linear convection-diffusion-reaction equations,
which stabilize continuous linear finite element discretizations for convection-domintated problems.
An example of such a problem is:

   &nabla; &sdot; ( &alpha;(x) &nabla;u + &beta;(x) u) + &gamma;(x) = f(x),

where &alpha;(x) is the diffusion term,
&beta;(x) is the convection term,
&gamma;(x) is the reaction term,
and *f* is the source term.

To learn more about EAFE, refer to the EAFE section below or the
[original article](http://www.ams.org/journals/mcom/1999-68-228/S0025-5718-99-01148-5/S0025-5718-99-01148-5.pdf) that introduces the method.

PyEAFE is designed to be easy to use, consisting of a single function, `eafe_assemble`,
that returns a stiffness matrix constructed using the EAFE approximation.
In order for this package to work,
the python implementation of [`dolfin`](https://fenicsproject.org/download/) is required.

The function signature for `eafe_assemble` (following syntax for the `typing` module) is:
```python
def eafe_assemble(
    mesh: Mesh,
    diffusion: pyeafe.Coefficient,
    convection: Optional[pyeafe.Coefficient] = None,
    reaction: Optional[pyeafe.Coefficient] = None,
) -> dolfin.Matrix:
```
where `pyeafe.Coefficient` is compatible with any of the following `dolfin` classes:
  - `Constant`
  - `Expression`
  - `Function`

**Note:**
that these classes are assumed to be imported from the `dolfin` module directly,
and `pyeafe` is not compatible with the `dolfin.cpp.function`.


## Example Usage

Example usage of `pyeafe` and `dolfin` can be found in some Jupyter notebooks found on
[PDE Labs](https://github.com/thepnpsolver/pdelab),
which provides a convergence analysis and examples in its `pyeafe` module.

For completeness, an example is provided below:
```python
from pyeafe import eafe_assemble
from dolfin import (
    assemble,
    dx,
    errornorm,
    Constant,
    DirichletBC,
    Expression,
    Function,
    FunctionSpace,
    LUSolver,
    TestFunction,
    UnitSquareMesh,
)

diffusion = Constant(0.5)
convection = Constant([1.0, 0.0])
source_term = Expression(
    "DOLFIN_PI * DOLFIN_PI * sin(DOLFIN_PI * x[0]) * sin(DOLFIN_PI * x[1]) \
    - DOLFIN_PI * cos(DOLFIN_PI * x[0]) * sin(DOLFIN_PI * x[1])",
    degree=4,
)
exact_solution = Expression(
    "sin(DOLFIN_PI * x[0]) * sin(DOLFIN_PI * x[1])", degree=4,
)

mesh = UnitSquareMesh(16, 16)
pw_linears = FunctionSpace(mesh, "Lagrange", 1)
test_function = TestFunction(pw_linears)

eafe_matrix = eafe_assemble(mesh, diffusion, convection)
rhs_vector = assemble(source_term * test_function * dx)

bc = DirichletBC(pw_linears, exact_solution, lambda _, on_bndry: on_bndry)
bc.apply(eafe_matrix, rhs_vector)

solution = Function(pw_linears)
solver = LUSolver(eafe_matrix, "default")
solver.parameters["symmetric"] = False
solver.solve(solution.vector(), rhs_vector)

l2_err: float = errornorm(exact_solution, solution, "l2", 3)
assert l2_err <= 6e-3, "L2 error too large"

h1_err: float = errornorm(exact_solution, solution, "H1", 3)
assert h1_err <= 2.2e-1, "H1 error too large"
```


## EAFE

Edge-Averaged Finite Elements (EAFE) are stable approximations
to finite element formulations for linear convection-reaction differential equations.
For any Delaunay mesh,
differential equations with a bounded and positive diffusivity coefficient,
finite convection coefficient,
and non-negative and bounded reaction term,
the EAFE discretization yields a monotone stiffness matrix
(see the
[original article](http://www.ams.org/journals/mcom/1999-68-228/S0025-5718-99-01148-5/S0025-5718-99-01148-5.pdf)
for supporting discussions).

Due to the general nature of the EAFE approximation, it can be applied to stabilize a broad variety of convection-dominated PDEs to compute continuous piecewise linear finite element solutions.

## Installation

To install `pyeafe`, run `pip install pyeafe`.
PyEAFE requires the python [dolfin](https://fenicsproject.org/download/)
package from the [FEniCS project](https://fenicsproject.org/) to be installed.
If you do not want to install `dolfin`,
use the docker image found on [PDE Labs](https://github.com/thepnpsolver/pdelab).
