"""
Assembly of stiffness matrix for a generic
convection diffusion reaction equation:
  diff in L1, conv in [L1]^d, reac in L_infty,

  -div(diff * grad(u) + conv * u) + reac * u.

Integral quantities are approximated by quadrature.

Usage:
    eafe_assemble(mesh, diff, conv, reac, boundary)

    - mesh: mesh defining finite element space
    - diff: DOLFIN expression for diffusion
    - conv: DOLFIN expression for convection
    - reac: [optional] DOLFIN expression for reaction
    - boundary: [optional] function taking in spatial
        coordinate and on_boundary boolean to determine
        whether coordinate is on a Dirichlet boundary
        (boundary condition is assumed to be zero)
"""

from __future__ import division
import numpy as np
import logging

from dolfin import (
    DirichletBC,
    FunctionSpace,
    TestFunction,
    TrialFunction,
    assemble,
    assemble_local,
    cells,
    dx,
    grad,
    inner,
    parameters,
)

from evaluate import create_safe_eval

__all__ = ["eafe_assemble"]


def bernoulli(r):
    if np.absolute(r) < 1e-10:
        return 1.0
    elif r < 0.0:
        return r / np.expm1(r)
    else:
        return r * np.exp(-r) / (1 - np.exp(-r))


def eafe_assemble(mesh, diff, conv=None, reac=None, boundary=None, **kwargs):
    logging.getLogger("FFC").setLevel(logging.WARNING)
    quadrature_degree = parameters["form_compiler"]["quadrature_degree"]
    parameters["form_compiler"]["quadrature_degree"] = 2
    spatial_dim = mesh.topology().dim()
    cell_vertex_count = spatial_dim + 1

    safe_diff = create_safe_eval(diff, 1, True)
    safe_conv = create_safe_eval(conv, spatial_dim)
    safe_reac = create_safe_eval(reac, 1)

    def edge_harmonic(start, edge, cell):
        midpt = start + 0.5 * edge
        return safe_diff(midpt, cell)

    def edge_psi(start, edge, cell):
        midpt = start + 0.5 * edge
        diffusion = safe_diff(midpt, cell)
        convection = safe_conv(midpt, cell)
        midpt_approx = np.inner(convection, edge)
        return bernoulli(-midpt_approx / diffusion)

    def lumped_reac(vertex, cell):
        reaction = safe_reac(vertex, cell)
        return reaction * cell.volume() / cell_vertex_count

    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v)) * dx
    A = assemble(a)
    A.zero()
    dof_map = V.dofmap()
    dof_coord = V.tabulate_dof_coordinates()

    for cell in cells(mesh):
        local_to_global_map = dof_map.cell_dofs(cell.index())
        local_tensor = assemble_local(a, cell)

        cell_vertex = np.empty([cell_vertex_count, spatial_dim])
        for local_dof in range(0, cell_vertex_count):
            dof_id = spatial_dim * local_to_global_map[local_dof]
            for coord in range(0, spatial_dim):
                cell_vertex[local_dof, coord] = dof_coord[dof_id + coord]

        for vertex_id in range(0, cell_vertex_count):
            vertex = cell_vertex[vertex_id]
            local_tensor[vertex_id, vertex_id] = lumped_reac(vertex, cell)
            for edge_id in range(0, cell_vertex_count):
                if edge_id == vertex_id:
                    continue

                edge = cell_vertex[edge_id] - vertex
                harmonic = edge_harmonic(vertex, edge, cell)
                psi = edge_psi(vertex, edge, cell)
                local_tensor[vertex_id, edge_id] *= harmonic * psi
                local_tensor[vertex_id, vertex_id] -= local_tensor[vertex_id, edge_id]

        A.add(local_tensor, local_to_global_map, local_to_global_map)
        A.apply("insert")

    if boundary is not None:
        bc = DirichletBC(V, 0.0, boundary)
        bc.apply(A)

    parameters["form_compiler"]["quadrature_degree"] = quadrature_degree
    return A
