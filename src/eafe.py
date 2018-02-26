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
    quadrature_degree = parameters["form_compiler"]["quadrature_degree"]
    parameters["form_compiler"]["quadrature_degree"] = 2

    if (reac is None):
        def reac(vertex, cell, **kwargs):
            return 0.0

    ##################################################
    # Mesh and FEM space
    ##################################################
    V = FunctionSpace(mesh, "Lagrange", 1)
    spatial_dim = mesh.topology().dim()
    cell_vertex_count = spatial_dim + 1

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
        local_tensor = assemble_local(a, cell)

        barycenter = np.zeros(spatial_dim)
        cell_vertex = np.empty([cell_vertex_count, spatial_dim])
        for local_dof in range(0, cell_vertex_count):
            vertex_id = spatial_dim * local_to_global_map[local_dof]
            for coord in range(0, spatial_dim):
                cell_vertex[local_dof, coord] = dof_coord[vertex_id + coord]
                barycenter[coord] += (cell_vertex[local_dof, coord]
                                      / cell_vertex_count)

        diff_val = diff(barycenter)
        beta = conv(barycenter)

        for vertex_id in range(0, cell_vertex_count):
            vertex = cell_vertex[vertex_id]

            try:
                local_tensor[vertex_id, vertex_id] = reac(vertex, cell)
            except:
                local_tensor[vertex_id, vertex_id] = reac(barycenter)

            local_tensor[vertex_id, vertex_id] *= (cell.volume()
                                                   / cell_vertex_count)

            for edge_id in range(0, cell_vertex_count):
                if (edge_id == vertex_id):
                    continue

                edge = vertex - cell_vertex[edge_id]
                eafe_weight = bernoulli1(np.inner(beta, edge), diff_val)
                local_tensor[vertex_id, edge_id] *= eafe_weight
                off_diagonal = local_tensor[vertex_id, edge_id]
                local_tensor[vertex_id, vertex_id] -= off_diagonal

        A.add(local_tensor, local_to_global_map, local_to_global_map)
        A.apply("insert")

    if boundary is not None:
        bc = DirichletBC(V, 0.0, boundary)
        bc.apply(A)

    parameters["form_compiler"]["quadrature_degree"] = quadrature_degree
    return A
