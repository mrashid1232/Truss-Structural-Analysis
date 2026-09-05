"""
crane_boom_example.py

Applies truss_solver.py to a representative cantilevered Pratt-truss boom,
in the same spirit as the truss-based crane boom I designed and built
physically for my Table-Top Crane project (see CV). This is an idealised
example geometry for demonstrating the solver, not a literal digital
twin of the physical build.

Boom: 4-bay cantilevered Pratt truss, 400 mm reach, 100 mm depth,
fixed at the wall end, carrying a tip load representing the crane's
target 3 kg maximum lift (worst case of the 1-3 kg design range).

Run:
    python crane_boom_example.py

Outputs:
    plots/boom_deformed_shape.png
    plots/boom_member_forces.png
    (results table printed to console)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from truss_solver import TrussModel

# ---------------------------------------------------------------
# 1. Assumed material & section (representative, not measured values
#    from the physical build -- stated explicitly for transparency)
# ---------------------------------------------------------------
E_TIMBER = 11e9          # Pa, representative softwood dowel along grain
DOWEL_DIAMETER = 0.008   # m (8 mm)
A_DOWEL = np.pi * (DOWEL_DIAMETER / 2) ** 2
ALLOWABLE_STRESS = 15e6  # Pa, representative allowable stress for timber dowel

TIP_MASS_KG = 3.0        # worst case of the 1-3 kg design range
G = 9.81
TIP_LOAD_N = TIP_MASS_KG * G

BAY_LENGTH = 0.100  # m
BOOM_DEPTH = 0.100  # m
N_BAYS = 4


def build_model() -> tuple[TrussModel, dict]:
    model = TrussModel()

    # Bottom chord: B0..B4 (B4 is the free tip)
    bottom = [model.add_node(i * BAY_LENGTH, 0.0) for i in range(N_BAYS + 1)]
    # Top chord: T0..T3 (tapers to a point at the tip, typical crane-boom shape)
    top = [model.add_node(i * BAY_LENGTH, BOOM_DEPTH) for i in range(N_BAYS)]

    # Chords
    for i in range(N_BAYS):
        model.add_element(bottom[i], bottom[i + 1], E_TIMBER, A_DOWEL)
    for i in range(N_BAYS - 1):
        model.add_element(top[i], top[i + 1], E_TIMBER, A_DOWEL)

    # Verticals
    for i in range(N_BAYS):
        model.add_element(bottom[i], top[i], E_TIMBER, A_DOWEL)

    # Diagonals (Pratt pattern: sloping down towards the free end)
    for i in range(N_BAYS - 1):
        model.add_element(top[i], bottom[i + 1], E_TIMBER, A_DOWEL)
    model.add_element(top[N_BAYS - 1], bottom[N_BAYS], E_TIMBER, A_DOWEL)

    # Wall-mounted base: fully fixed (both chords pinned into the wall)
    model.fix(bottom[0], x=True, y=True)
    model.fix(top[0], x=True, y=True)

    # Tip load: downward force from the lifted mass, at the free tip node
    model.add_load(bottom[N_BAYS], fx=0.0, fy=-TIP_LOAD_N)

    node_positions = {"bottom": bottom, "top": top}
    return model, node_positions


def plot_deformed_shape(model, result, scale=200, path="plots/boom_deformed_shape.png"):
    fig, ax = plt.subplots(figsize=(8, 4))

    def draw(nodes, disp_scale, color, label, lw):
        for el in model.elements:
            ni, nj = model.nodes[el.node_i], model.nodes[el.node_j]
            if disp_scale:
                ui = result.displacements[2*el.node_i:2*el.node_i+2] * disp_scale
                uj = result.displacements[2*el.node_j:2*el.node_j+2] * disp_scale
                x = [ni.x + ui[0], nj.x + uj[0]]
                y = [ni.y + ui[1], nj.y + uj[1]]
            else:
                x = [ni.x, nj.x]
                y = [ni.y, nj.y]
            ax.plot(x, y, color=color, linewidth=lw, zorder=2 if disp_scale else 1)

    draw(model.nodes, 0, "0.75", "Undeformed", 1.5)
    draw(model.nodes, scale, "#1F3864", f"Deformed (x{scale})", 2.2)

    ax.scatter([n.x for n in model.nodes], [n.y for n in model.nodes],
               color="#1F3864", zorder=3, s=15)
    ax.set_aspect("equal")
    ax.set_title(f"Cantilevered crane boom — tip load {TIP_LOAD_N:.1f} N "
                 f"({TIP_MASS_KG:.0f} kg), displacement magnified {scale}x")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.legend(["Undeformed", f"Deformed (x{scale})"], loc="lower left")
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_member_forces(model, result, path="plots/boom_member_forces.png"):
    fig, ax = plt.subplots(figsize=(9, 4.5))

    forces = result.member_forces
    max_abs = np.max(np.abs(forces))
    segments = []
    colors = []
    for el, f in zip(model.elements, forces):
        ni, nj = model.nodes[el.node_i], model.nodes[el.node_j]
        segments.append([(ni.x, ni.y), (nj.x, nj.y)])
        colors.append(f)

    lc = LineCollection(segments, cmap="coolwarm", linewidths=4)
    lc.set_array(np.array(colors))
    lc.set_clim(-max_abs, max_abs)
    ax.add_collection(lc)
    ax.scatter([n.x for n in model.nodes], [n.y for n in model.nodes],
               color="black", zorder=3, s=15)

    cbar = fig.colorbar(lc, ax=ax, orientation="vertical", fraction=0.04, pad=0.03)
    cbar.set_label("Member axial force (N)\nblue = compression, red = tension", fontsize=9)

    ax.set_aspect("equal")
    ax.set_title("Member forces under 3 kg tip load")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    model, _ = build_model()
    result = model.solve()

    print(f"Tip load: {TIP_MASS_KG:.1f} kg -> {TIP_LOAD_N:.2f} N\n")
    print(result.summary(allowable_stress=ALLOWABLE_STRESS))

    fos_all = result.factor_of_safety(ALLOWABLE_STRESS)
    worst = int(np.argmin(fos_all))
    print(f"\nCritical member: {worst} "
          f"(node {model.elements[worst].node_i} -> {model.elements[worst].node_j}), "
          f"FOS = {fos_all[worst]:.2f}")

    plot_deformed_shape(model, result)
    plot_member_forces(model, result)
    print("\nPlots written to plots/boom_deformed_shape.png and plots/boom_member_forces.png")


if __name__ == "__main__":
    main()
