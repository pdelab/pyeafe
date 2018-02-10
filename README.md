# PyEAFE

## Overview

Edge-Averaged Finite Elements (EAFE) are stable approximations to standard finite element formulations for linear convection-diffusion-reaction differential equations.
For any Delaunay mesh, bounded and positive diffusivity coefficient, finite convection coefficient, and non-negative and bounded reaction term, the EAFE discretization yields a monotone stiffness matrix (see [here](http://www.ams.org/journals/mcom/1999-68-228/S0025-5718-99-01148-5/S0025-5718-99-01148-5.pdf)).
Due to the general nature of the EAFE approximation, it can be applied to stabilize a broad variety of convection-dominated PDEs to compute continuous piecewise linear finite element solutions.

**PyEAFE** implements the EAFE approximation for linear convection-diffusion-reaction equations with PDE finite coefficients based on the [Dolfin software package](https://fenicsproject.org/).
For coupled PDEs and nonlinear problems, helper methods exist in the module that allow users to define the PDE coefficients in terms of finite element functions and their derivatives;
this is usefule for applying EAFE approximations to linearized systems resulting from a Newton iteration scheme or solving weakly coupled differential equations.

## Getting started

To ensure that the git hooks are properly being used (and creating a docker image for MacOS), run:
```
  ./scripts/init
```

When contributing, ensure code styling is consistent by running `./scripts/lint`.
This script is automatically run whenever a commit is made.
Use `git commit ... --no-verify` to bypass linting, although this practice is strongly discouraged.

### MacOS

The FEniCS distribution uses [Docker](https://www.docker.com/) for simplicity.
Install Docker by following [these instructions](https://docs.docker.com/docker-for-mac/install/).
To run the project in a docker container, run
```
  ./scripts/start
```
The Docker container will open in a shared directory to the repository's root so that any updates to files in the repository are available from within the Docker container.

To exit the Docker container, simply run the `exit` command.

### Linux systems

If FEniCS is not already installed or there are compatibility issues, follow the steps for MacOS to run FEniCS in a Docker container.
