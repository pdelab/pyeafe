""" The primary function of the library: eafe_assemble """

from __future__ import division
import numpy as np
import logging
import dolfin
from petsc4py import PETSc
from dolfin import (
    FunctionSpace,
    Mesh,
    TestFunction,
    TrialFunction,
    as_backend_type,
    assemble,
    assemble_local,
    cells,
    dx,
    grad,
    inner,
    parameters,
)
from typing import Callable, Optional

from pyeafe.utils import Coefficient, ensure_edge_eval


def bernoulli(r: float) -> float:
    """ evaluate the Bernoulli function at the given value """

    if np.absolute(r) < 1e-10:
        return 1.0

    if r < 0.0:
        return r / np.expm1(r)

    return r * np.exp(-r) / (1 - np.exp(-r))


def define_edge_advection(
    dim: int, diffusion: Coefficient, convection: Optional[Coefficient] = None
) -> Callable[[dolfin.Vertex, dolfin.Edge, dolfin.Cell], float]:
    """
    Get most efficient method for computing the Edge-Averaged flux value.

    :param dim: integer dimension of the convection term
    :param diffusion: Coefficient for diffusivity (should be positive)
    :param convection: [Optional] Coefficient for convection term
                        (must evaluate to an array of size `dim`)

    :return: Method for computing flux contribution from a cell across an edge:
        Callable[[dolfin.Vertex, dolfin.Edge, dolfin.Cell], float]
    """

    diffusion_with_edge_eval = ensure_edge_eval(diffusion, 1)
    if convection is None:

        def edge_harmonic(start, edge, cell):
            midpt = start + 0.5 * edge
            return diffusion_with_edge_eval(midpt, cell)

        return edge_harmonic

    convection_with_edge_eval = ensure_edge_eval(convection, dim)
    if dim > 1:
        conv_rank: int = convection.value_rank()
        if conv_rank != 1:
            raise ValueError("Invalid convection parameter: value_rank must be 1")

        if convection.value_dimension(0) != dim:
            raise ValueError(
                "Invalid convection parameter: value_dimension(0) must match mesh spatial dimension"  # noqa E501
            )

    def edge_psi(start, edge, cell):
        midpt = start + 0.5 * edge
        diff = diffusion_with_edge_eval(midpt, cell)
        conv = convection_with_edge_eval(midpt, cell)
        midpt_approx = np.inner(conv, edge)
        return diff * bernoulli(-midpt_approx / diff)

    return edge_psi


def define_mass_lumping(
    cell_vertex_count: int,
    reaction: Optional[Coefficient] = None,
) -> callable:
    """
    Return the most efficient method for providing the values on the diagonal
    pertaining to the optionally given reaction coefficient.

    :param cell_vertex_count: number of vertices in a given cell
    :param reaction: [Optional] Coefficient for reaction term.

    :return: Method for geting a cell's contribution to the mass-lumped reaction
        stiffness matrix:
        (vertex: dolfin.Vertex, cell: dolfin.Cell) -> float
    """

    if reaction is None:
        return lambda v, c: 0.0

    reaction_with_edge_eval = ensure_edge_eval(reaction, 1)

    def lumped_reac(vertex, cell):
        return reaction_with_edge_eval(vertex, cell) * cell.volume() / cell_vertex_count

    return lumped_reac


def eafe_assemble(
    mesh: Mesh,
    diffusion: Coefficient,
    convection: Optional[Coefficient] = None,
    reaction: Optional[Coefficient] = None,
) -> dolfin.Matrix:
    """
    Assembly of stiffness matrix for a generic linear elliptic equation:
    for `u` a continuous piecewise-linear finite element function,
      -div(diffusion * grad(u) + convection * u) + reaction * u = source

    Edge-integral quantities are approximated by a midpoint quadrature rule.

    :param mesh: Mesh defining finite element space
    :param diff: Coefficient for diffusion coefficient
    :param conv: [Optional] Coefficient for convection coefficient
    :param reac: [Optional] Coefficient for reaction coefficient

    :return: dolfin.Matrix pertaining to EAFE stiffness matrix
    """

    if not issubclass(mesh.__class__, Mesh):
        raise TypeError("Invalid mesh parameter: must inherit from dolfin.Mesh")

    logging.getLogger("FFC").setLevel(logging.WARNING)
    quadrature_degree = parameters["form_compiler"]["quadrature_degree"]
    parameters["form_compiler"]["quadrature_degree"] = 2

    spatial_dim: int = mesh.topology().dim()
    cell_vertex_count: int = spatial_dim + 1

    try:
        edge_advection = define_edge_advection(spatial_dim, diffusion, convection)
        lumped_reaction = define_mass_lumping(cell_vertex_count, reaction)
    except Exception as e:
        parameters["form_compiler"]["quadrature_degree"] = quadrature_degree
        raise e

    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v)) * dx
    A: dolfin.Matrix = assemble(a)
    mat: dolfin.PETScMatrix = as_backend_type(A).mat()

    A.zero()
    dof_map = V.dofmap()
    dof_coord = V.tabulate_dof_coordinates()

    for cell in cells(mesh):
        local_tensor = assemble_local(a, cell)
        local_to_global_map = dof_map.cell_dofs(cell.index())
        cell_vertex = np.empty([cell_vertex_count, spatial_dim])

        for local_dof in range(0, cell_vertex_count):
            dof_id = local_to_global_map[local_dof]
            for coord in range(0, spatial_dim):
                cell_vertex[local_dof, coord] = dof_coord[dof_id, coord]

        for vertex_id in range(0, cell_vertex_count):
            vertex = cell_vertex[vertex_id]
            local_tensor[vertex_id, vertex_id] = lumped_reaction(vertex, cell)

            for edge_id in range(0, cell_vertex_count):
                if edge_id == vertex_id:
                    continue

                edge = cell_vertex[edge_id] - vertex
                local_tensor[vertex_id, edge_id] *= edge_advection(vertex, edge, cell)
                local_tensor[vertex_id, vertex_id] -= local_tensor[vertex_id, edge_id]

        local_to_global_list = local_to_global_map.tolist()
        mat.setValuesLocal(
            local_to_global_list,
            local_to_global_list,
            local_tensor,
            PETSc.InsertMode.ADD,
        )

    A.apply("insert")
    mat.assemble()

    parameters["form_compiler"]["quadrature_degree"] = quadrature_degree
    return A
