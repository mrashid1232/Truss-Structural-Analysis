"""
test_solver.py

Validates TrussModel/solve() against a classic statically-determinate
3-bar "A-frame" truss with a hand-calculated (method of joints) solution.

Geometry:
    A = (0, 0)  -- pin support
    B = (4, 0)  -- roller support (vertical reaction only)
    C = (2, 2)  -- apex, 1000 N downward point load

Hand calculation (method of joints):
    Reactions:      R_A = (0, 500) N,   R_B = (0, 500) N
    Member AB:      +500.00 N   (tension)
    Members AC, BC: -707.11 N  each (compression)

If this test passes, the solver's assembly, boundary-condition handling,
and member-force recovery are all verified against an independent,
hand-checkable result.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from truss_solver import TrussModel

# arbitrary but consistent E, A -- result is independent of these for a
# statically determinate truss (only affects displacements, not forces)
E = 200e9   # Pa (steel)
A = 0.0005  # m^2


def test_three_bar_truss_matches_hand_calculation():
    model = TrussModel()
    a = model.add_node(0.0, 0.0)
    b = model.add_node(4.0, 0.0)
    c = model.add_node(2.0, 2.0)

    model.add_element(a, c, E, A)   # member 0: AC
    model.add_element(b, c, E, A)   # member 1: BC
    model.add_element(a, b, E, A)   # member 2: AB

    model.fix(a, x=True, y=True)    # pin
    model.fix(b, x=False, y=True)   # roller

    model.add_load(c, fx=0.0, fy=-1000.0)

    result = model.solve()

    f_ac, f_bc, f_ab = result.member_forces

    assert np.isclose(f_ab, 500.0, atol=0.5), f"AB expected +500 N, got {f_ab:.2f} N"
    assert np.isclose(f_ac, -707.11, atol=0.5), f"AC expected -707.11 N, got {f_ac:.2f} N"
    assert np.isclose(f_bc, -707.11, atol=0.5), f"BC expected -707.11 N, got {f_bc:.2f} N"

    r_ay = result.reactions[2 * a + 1]
    r_by = result.reactions[2 * b + 1]
    assert np.isclose(r_ay, 500.0, atol=0.5)
    assert np.isclose(r_by, 500.0, atol=0.5)

    print("PASS: solver matches hand-calculated 3-bar truss to within 0.5 N")
    print(f"  AB (tension):        {f_ab:8.2f} N   (expected +500.00 N)")
    print(f"  AC (compression):    {f_ac:8.2f} N   (expected -707.11 N)")
    print(f"  BC (compression):    {f_bc:8.2f} N   (expected -707.11 N)")
    print(f"  Reaction at A (Fy):  {r_ay:8.2f} N   (expected +500.00 N)")
    print(f"  Reaction at B (Fy):  {r_by:8.2f} N   (expected +500.00 N)")


if __name__ == "__main__":
    test_three_bar_truss_matches_hand_calculation()
