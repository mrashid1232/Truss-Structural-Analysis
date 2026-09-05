"""
truss_solver.py

A 2D truss structural analysis engine using the Direct Stiffness Method
(the matrix formulation of finite element analysis for pin-jointed frames).

Given a truss's geometry, member properties, supports and applied loads,
this computes:
    - nodal displacements
    - support reactions
    - axial force in every member (tension positive, compression negative)
    - axial stress in every member and a factor of safety against a
      user-supplied allowable stress

Unlike hand methods (method of joints / method of sections), the stiffness
method handles statically indeterminate trusses as well as determinate
ones, since it solves compatibility and equilibrium simultaneously via
the global stiffness matrix rather than requiring a staticall determinate
free-body chain.

Author: Mahbubur Rashid
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field


@dataclass
class Node:
    x: float  # m
    y: float  # m


@dataclass
class Element:
    node_i: int
    node_j: int
    E: float   # Young's modulus, Pa
    A: float   # cross-sectional area, m^2

    def length(self, nodes: list[Node]) -> float:
        ni, nj = nodes[self.node_i], nodes[self.node_j]
        return float(np.hypot(nj.x - ni.x, nj.y - ni.y))

    def direction_cosines(self, nodes: list[Node]) -> tuple[float, float]:
        ni, nj = nodes[self.node_i], nodes[self.node_j]
        L = self.length(nodes)
        return (nj.x - ni.x) / L, (nj.y - ni.y) / L


class TrussModel:
    """A 2D pin-jointed truss ready for direct-stiffness-method analysis."""

    def __init__(self):
        self.nodes: list[Node] = []
        self.elements: list[Element] = []
        self.fixed_dofs: set[int] = set()      # global dof indices that are restrained
        self.loads: dict[int, float] = {}       # global dof index -> applied force (N)

    # ---------- model building ----------

    def add_node(self, x: float, y: float) -> int:
        self.nodes.append(Node(x, y))
        return len(self.nodes) - 1

    def add_element(self, node_i: int, node_j: int, E: float, A: float) -> int:
        self.elements.append(Element(node_i, node_j, E, A))
        return len(self.elements) - 1

    def fix(self, node: int, x: bool = True, y: bool = True) -> None:
        """Restrain a node's x and/or y displacement (pin = both True, roller = one True)."""
        if x:
            self.fixed_dofs.add(2 * node)
        if y:
            self.fixed_dofs.add(2 * node + 1)

    def add_load(self, node: int, fx: float = 0.0, fy: float = 0.0) -> None:
        self.loads[2 * node] = self.loads.get(2 * node, 0.0) + fx
        self.loads[2 * node + 1] = self.loads.get(2 * node + 1, 0.0) + fy

    # ---------- assembly & solve ----------

    def _global_stiffness(self) -> np.ndarray:
        n_dof = 2 * len(self.nodes)
        K = np.zeros((n_dof, n_dof))
        for el in self.elements:
            L = el.length(self.nodes)
            c, s = el.direction_cosines(self.nodes)
            k_axial = el.E * el.A / L
            k_local = k_axial * np.array([
                [ c*c,  c*s, -c*c, -c*s],
                [ c*s,  s*s, -c*s, -s*s],
                [-c*c, -c*s,  c*c,  c*s],
                [-c*s, -s*s,  c*s,  s*s],
            ])
            dofs = [2*el.node_i, 2*el.node_i+1, 2*el.node_j, 2*el.node_j+1]
            for a in range(4):
                for b in range(4):
                    K[dofs[a], dofs[b]] += k_local[a, b]
        return K

    def solve(self) -> "TrussResult":
        n_nodes = len(self.nodes)
        n_dof = 2 * n_nodes
        K = self._global_stiffness()

        F = np.zeros(n_dof)
        for dof, val in self.loads.items():
            F[dof] = val

        free = np.array(sorted(set(range(n_dof)) - self.fixed_dofs))
        fixed = np.array(sorted(self.fixed_dofs))

        if len(free) == 0:
            raise ValueError("Truss has no free degrees of freedom.")

        K_ff = K[np.ix_(free, free)]
        F_f = F[free]

        # Basic stability check: singular K_ff means a mechanism (unstable truss)
        if np.linalg.matrix_rank(K_ff) < K_ff.shape[0]:
            raise ValueError(
                "Global stiffness matrix is singular on the free DOFs — "
                "the truss is a mechanism (under-restrained or a node is "
                "not adequately triangulated). Check supports/members."
            )

        U = np.zeros(n_dof)
        U[free] = np.linalg.solve(K_ff, F_f)

        reactions_full = K @ U - F
        reactions = {dof: reactions_full[dof] for dof in fixed}

        member_forces = []
        member_stresses = []
        for el in self.elements:
            L = el.length(self.nodes)
            c, s = el.direction_cosines(self.nodes)
            dofs = [2*el.node_i, 2*el.node_i+1, 2*el.node_j, 2*el.node_j+1]
            u_local = U[dofs]
            # axial force from elongation: F = (EA/L) * (elongation)
            elongation = -c*u_local[0] - s*u_local[1] + c*u_local[2] + s*u_local[3]
            force = el.E * el.A / L * elongation
            member_forces.append(force)
            member_stresses.append(force / el.A)

        return TrussResult(
            displacements=U,
            reactions=reactions,
            member_forces=np.array(member_forces),
            member_stresses=np.array(member_stresses),
            model=self,
        )


@dataclass
class TrussResult:
    displacements: np.ndarray
    reactions: dict[int, float]
    member_forces: np.ndarray
    member_stresses: np.ndarray
    model: TrussModel

    def factor_of_safety(self, allowable_stress: float) -> np.ndarray:
        """FOS per member against a given allowable stress (Pa), tension and compression alike.
        Zero-force members (e.g. redundant members under a particular load case)
        report FOS as +inf, which is the physically correct answer, not an error."""
        with np.errstate(divide="ignore"):
            return allowable_stress / np.abs(self.member_stresses)

    def summary(self, allowable_stress: float | None = None) -> str:
        lines = ["Member forces (tension +ve, compression -ve):"]
        for i, (f, s) in enumerate(zip(self.member_forces, self.member_stresses)):
            state = "tension" if f > 0 else ("compression" if f < 0 else "zero-force")
            line = f"  Member {i:2d} ({self.model.elements[i].node_i}->{self.model.elements[i].node_j}): " \
                   f"{f:8.2f} N  ({state}), stress = {s/1e6:6.2f} MPa"
            if allowable_stress:
                with np.errstate(divide="ignore"):
                    fos = allowable_stress / abs(s) if s != 0 else float("inf")
                line += f",  FOS = {fos:5.2f}" if fos != float("inf") else ",  FOS =   inf"
            lines.append(line)
        lines.append("\nReactions:")
        for dof, val in sorted(self.reactions.items()):
            node, axis = divmod(dof, 2)
            lines.append(f"  Node {node} {'Fx' if axis == 0 else 'Fy'}: {val:8.2f} N")
        return "\n".join(lines)
