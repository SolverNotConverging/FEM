"""Resource checks on candidate meshes before assembling their operators."""

from fem.common.validation import positive_integer


class RefinementLimit(RuntimeError):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def resource_limits(max_elements, max_dofs):
    return tuple(None if value is None else positive_integer(value, name)
                 for value, name in ((max_elements, "max_elements"), (max_dofs, "max_dofs")))


def check_mesh_policy(owner, mesh, previous_count, physical=None):
    quality, growth = getattr(owner, "_candidate_policy", (None, None))
    if growth is not None and mesh.nelements > previous_count * growth:
        raise RefinementLimit("growth_limit", f"Candidate grew from {previous_count} to {mesh.nelements} cells.")
    if quality is not None:
        from .policy import mesh_quality
        minimum = mesh_quality(mesh if physical is None else physical)["minimum_quality"]
        if minimum < quality:
            raise RefinementLimit("quality_limit", f"Candidate minimum quality {minimum:g} is below {quality:g}.")


def check_candidate(owner, mesh, order, mixed=False, *, physical=None):
    limits = getattr(owner, "_adaptive_limits", None)
    if limits is None:
        return
    previous_count = getattr(owner, "_previous_cell_count", None)
    if previous_count is not None:
        check_mesh_policy(owner, mesh, previous_count, physical)
    elements, dofs = limits
    if elements is not None and mesh.nelements > elements:
        raise RefinementLimit("element_limit", f"Candidate has {mesh.nelements} elements (budget {elements}).")
    if dofs is None:
        return
    count = field_dofs(mesh, order, mixed)
    if count > dofs:
        raise RefinementLimit("dof_limit", f"Candidate has {count} field DOFs (budget {dofs}).")


def field_dofs(mesh, order, mixed=False):
    if mesh.dim() == 1:
        count = mesh.nvertices + (order - 1) * mesh.nelements
    elif mesh.dim() == 2:
        count = mesh.nvertices + (order - 1) * mesh.nfacets
        if mixed:
            count += order * mesh.nfacets + 2 * (order - 1) * mesh.nelements
    else:
        count = order * mesh.nedges + 2 * (order - 1) * mesh.nfacets
    return count
