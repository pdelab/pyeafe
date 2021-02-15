# PyEAFE

This package is an implementation of the edge-averaged finite element (EAFE) approximation
to linear convection-diffusion-reaction equations,
which is use to stabilize discretizations for convection-domintated problems.
To learn more about EAFE, refer to the EAFE section below.

PyEAFE is designed to be easy to use, consisting of a single function, `eafe_assemble`,
that returns a stiffness matrix constructed using the EAFE approximation.
In order for this package to work, the [python implementation of `dolfin`](https://fenicsproject.org/download/)
is required.
The function signature for `eafe_assemble` (following syntax for the `typing` module) is:
```python
def eafe_assemble(
  mesh: dolfin.Mesh,
  diff: dolfin.Expression,
  conv: Optional[dolfin.Expression] = None,
  reac: Optional[dolfin.Expression] = None,
) -> dolfin.Matrix:
```

## EAFE

Edge-Averaged Finite Elements (EAFE) are stable approximations to standard finite element formulations for linear convection-diffusion-reaction differential equations.
For any Delaunay mesh, differential equations with a bounded and positive diffusivity coefficient, finite convection coefficient, and non-negative and bounded reaction term, the EAFE discretization yields a monotone stiffness matrix (see [here](http://www.ams.org/journals/mcom/1999-68-228/S0025-5718-99-01148-5/S0025-5718-99-01148-5.pdf)).
Due to the general nature of the EAFE approximation, it can be applied to stabilize a broad variety of convection-dominated PDEs to compute continuous piecewise linear finite element solutions.

## Installation

To install run
```
  pip install pyeafe
```

PyEAFE requires the python [dolfin](https://fenicsproject.org/download/) package
from the [FEniCS project](https://fenicsproject.org/) to be installed.


If you do not want to install `dolfin`, use the docker image found here: [PDELAB](https://github.com/thepnpsolver/pdelab).
