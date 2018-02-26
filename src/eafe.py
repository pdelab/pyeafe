'''
Assembly of stiffness matrix for a generic
convection diffusion reaction equation:
  diff in L1, conv in [L1]^d, reac in L_infty,

  -div(diff * grad(u) + conv * u) + reac * u.

Integral quantities are approximated by quadrature.

Usage:
    eafe(mesh, diff, conv, reac, boundary)

    - mesh: mesh defining finite element space
    - diff: DOLFIN expression for diffusion
    - conv: DOLFIN expression for convection
    - reac: [optional] DOLFIN expression for reaction
    - boundary: [optional] function taking in spatial
        coordinate and on_boundary boolean to determine
        whether coordinate is on a Dirichlet boundary
        (boundary condition is assumed to be zero)
'''

from __future__ import division
from dolfin import *
import numpy as np
import sys


def bernoulli1(r, diffusion_value):
    eps = 1e-10
    if (np.absolute(r) < diffusion_value * eps):
        return diffusion_value
    elif (r < -diffusion_value * eps):
        return r / np.expm1(r/diffusion_value)
    else:
        return (r * np.exp(-r/diffusion_value)
                / (1 - np.exp(-r/diffusion_value)))


def eafe(mesh, diff, conv, reac=None, boundary=None, **kwargs):

    if (reac is not None):
        print "reaction terms are not yet supported"

    ##################################################
    # Mesh and FEM space
    ##################################################
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v)) * dx

    ##################################################
    # Build the stiffness matrix
    ##################################################
    A = assemble(a)
    A.zero()
    dof_map = V.dofmap()
    dof_coord = V.tabulate_dof_coordinates()

    for cell in cells(mesh):
        local_to_global_map = dof_map.cell_dofs(cell.index())
        # build the local tensor
        local_tensor = assemble_local(a, cell)
        # EAFE: change the local tensor
        # Step 1: Find the point related to dofs
        a0 = np.array([dof_coord[2*local_to_global_map[0]],
                       dof_coord[2*local_to_global_map[0]+1]])
        a1 = np.array([dof_coord[2*local_to_global_map[1]],
                       dof_coord[2*local_to_global_map[1]+1]])
        a2 = np.array([dof_coord[2*local_to_global_map[2]],
                       dof_coord[2*local_to_global_map[2]+1]])
        barycenter = (a0+a1+a2) / 3
        # Step 2: Find the convection by local constant approximation
        beta = conv(barycenter)
        # Step 3: Apply bernoulli function
        diff_val = diff(barycenter)
        b01 = bernoulli1(np.inner(beta, a0-a1), diff_val)
        b10 = bernoulli1(np.inner(beta, a1-a0), diff_val)
        b02 = bernoulli1(np.inner(beta, a0-a2), diff_val)
        b20 = bernoulli1(np.inner(beta, a2-a0), diff_val)
        b12 = bernoulli1(np.inner(beta, a1-a2), diff_val)
        b21 = bernoulli1(np.inner(beta, a2-a1), diff_val)
        # Step 4: Change the local tensor
        local_tensor[0][1] *= b01
        local_tensor[1][0] *= b10
        local_tensor[0][2] *= b02
        local_tensor[2][0] *= b20
        local_tensor[1][2] *= b12
        local_tensor[2][1] *= b21
        local_tensor[0][0] = -local_tensor[1][0] - local_tensor[2][0]
        local_tensor[1][1] = -local_tensor[0][1] - local_tensor[2][1]
        local_tensor[2][2] = -local_tensor[0][2] - local_tensor[1][2]
        # Build the stiffness matrix
        A.add(local_tensor, local_to_global_map, local_to_global_map)
        A.apply("insert")

    if boundary is not None:
        bc = DirichletBC(V, 0.0, boundary)
        bc.apply(A)

    return A
