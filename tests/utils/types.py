from dolfin import Mesh
from typing import List, Optional

from pyeafe import eafe_assemble, Coefficient


def assemble_with_mixed_types(
    mesh: Mesh,
    diffusions: List[Coefficient],
    convections: List[Optional[Coefficient]],
    reactions: List[Optional[Coefficient]],
):
    for diffusion in diffusions:
        for convection in convections:
            for reaction in reactions:
                eafe_assemble(mesh, diffusion, convection, reaction)
