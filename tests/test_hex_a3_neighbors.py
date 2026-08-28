"""Exact aperture-3 neighbours (A0): py4dggs hex_a3 vs pydggal getZoneNeighbors.

The I3H packed value equals pydggal's DGGRSZone int, so neighbours are compared
by EXACT int-set equality — the strongest possible bar (no rounding, any seed).
"""
import random
import pytest
from py4dggs import ISEA3H, IGEO7, Zone
from py4dggs.topologies.hex_a3 import move5x6_vertex_v2, _i3h_get_neighbors
from py4dggs.indexings.i3h import unpack_i3h
from _pydggal_oracle import (requires_pydggal, oracle_grid, oracle_zone_int,
                             neighbour_ints_of_int)

pytestmark = requires_pydggal


def _py4dggs_nb(val):
    return set(_i3h_get_neighbors(*unpack_i3h(val)))


# --- Task 1: the oracle helper -------------------------------------------- #
def test_oracle_helper_returns_int_set():
    g = oracle_grid("ISEA3H")
    val = oracle_zone_int(g, 0.0, 0.0, 5)
    nb = neighbour_ints_of_int(g, val)
    assert nb is not None and all(isinstance(n, int) for n in nb) and val not in nb


# --- Task 2: the move5x6Vertex2 (crossEarly v2) primitive ------------------ #
def test_move5x6_v2_noncrossing_is_plain_offset():
    # interior point, small offset: no interruption crossing -> v = c + d exactly
    assert move5x6_vertex_v2(2.4, 2.5, 0.1, 0.05, True) == [2.5, 2.55]


def test_move5x6_v2_snaps_to_poles():
    # A crossing offset landing on a polar diagonal snaps to the pole vertex.
    # These are the real offsets an even-level topLeft/bottomRight neighbour of a
    # pole-adjacent cell produces (centroid (2/3,0) -> "North" (1,0); the mirror
    # case -> "South" (4,6)).
    assert move5x6_vertex_v2(2 / 3, 0.0, 0.0, -1 / 3, True) == [1, 0]
    assert move5x6_vertex_v2(2 + 1 / 3, 4.0, 1 / 3, 2 / 3, True) == [4, 6]


# --- Task 3: the getNeighbor/getNeighbors port (module level) -------------- #
@pytest.mark.parametrize("lat,lon,res", [
    (0.0, 0.0, 5), (48.0, 11.2, 6), (89.99, 0.0, 3), (-89.99, 0.0, 4),
    (-58.4, -168.8, 12),   # the cell the geometric k-ring MISSES a neighbour for
    (-58.28, 10.0, 8), (30.0, 108.8, 9),
])
def test_hex_a3_neighbors_exact_vs_pydggal(lat, lon, res):
    g = oracle_grid("ISEA3H")
    val = oracle_zone_int(g, lat, lon, res)
    assert _py4dggs_nb(val) == neighbour_ints_of_int(g, val)


# --- Task 4: the Grid.neighbors dispatch ----------------------------------- #
def test_grid_dispatch_uses_topology_override_for_3h():
    g = oracle_grid("ISEA3H")
    val = oracle_zone_int(g, -58.4, -168.8, 12)   # k-ring got this WRONG before A0
    got = {z.value for z in Zone(ISEA3H, val).neighbors}
    assert got == neighbour_ints_of_int(g, val)


def test_grid_dispatch_leaves_aperture7_on_edge_kring():
    # aperture-7 has NO topology.neighbors override -> keeps the edge k-ring
    assert getattr(IGEO7.topology, "neighbors", None) is None


# --- Regression guard: the risky-cell classes, exact, any seed ------------- #
def test_risky_cells_exact_all_3h_grids():
    """Pentagons, icosahedron apices (+-58.28), interruption-seam longitudes,
    and high resolutions — the exact classes the edge k-ring overshoots on."""
    from py4dggs import IVEA3H, RTEA3H
    grids = {"ISEA3H": ISEA3H, "IVEA3H": IVEA3H, "RTEA3H": RTEA3H}
    risky = [(-58.4, -168.8, 12), (-58.28, 10.0, 8), (58.28, 11.2, 6),
             (89.99, 0.0, 3), (-89.99, 0.0, 4), (0.0, 11.2, 5),
             (0.0, -168.8, 7), (30.0, 108.8, 9), (-30.0, -71.2, 11)]
    for name, grid in grids.items():
        g = oracle_grid(name)
        for lat, lon, res in risky:
            val = oracle_zone_int(g, lat, lon, res)
            got = {z.value for z in Zone(grid, val).neighbors}
            assert got == neighbour_ints_of_int(g, val), f"{name} ({lat},{lon},r{res})"


def test_high_resolution_exact_all_3h_grids():
    """Exact neighbours at very high resolution (19-33, up to the max), where
    d/do3 shrink toward the 1e-11/1e-10 epsilon scale — the domain the fuzz
    tests (res <= 18) do not reach. Confirms no epsilon-scale branch breaks."""
    from py4dggs import IVEA3H, RTEA3H
    grids = {"ISEA3H": ISEA3H, "IVEA3H": IVEA3H, "RTEA3H": RTEA3H}
    for name, grid in grids.items():
        g = oracle_grid(name)
        rng = random.Random(31337)
        n_fail = 0
        for _ in range(300):
            lat = rng.uniform(-89.0, 89.0); lon = rng.uniform(-180.0, 180.0)
            res = rng.randint(19, 33)
            val = oracle_zone_int(g, lat, lon, res)
            onb = neighbour_ints_of_int(g, val)
            if onb is None:
                continue
            if {z.value for z in Zone(grid, val).neighbors} != onb:
                n_fail += 1
        assert n_fail == 0, f"{name}: {n_fail}/300 high-res neighbour mismatches"
