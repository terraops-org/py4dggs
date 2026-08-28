"""Live differential fuzz: py4dggs.ISEA3H vs pydggal ISEA3H (point-keyed geometry).

Neighbours are compared by EXACT int-set equality (I3H packed value == pydggal
DGGRSZone int) since A0 gave hex_a3 an exact topological `neighbors` override.
"""
import random
import pytest
from py4dggs import ISEA3H, Zone
from _pydggal_oracle import (requires_pydggal, oracle_grid, oracle_zone_int,
                             centroid_of_int, vertices_of_int, neighbour_ints_of_int)

pytestmark = requires_pydggal
TOL = 1e-7

def _srt(pairs, nd=4):
    return sorted((round(a, nd), round(b, nd)) for a, b in pairs)

def test_fuzz_geometry_matches_pydggal():
    g = oracle_grid("ISEA3H"); rng = random.Random(20260702)
    c_fail = v_fail = n_fail = 0; nb_checked = 0
    for _ in range(500):
        lat = rng.uniform(-88, 88); lon = rng.uniform(-180, 180); res = rng.randint(1, 11)
        val = oracle_zone_int(g, lat, lon, res)
        z = Zone(ISEA3H, val)                    # DGGAL-layout int == our value
        # centroid
        olat, olon = centroid_of_int(g, val); c = z.centroid
        if abs(c.lat - olat) > TOL or abs(c.lon - olon) > TOL: c_fail += 1
        # vertices (set + count)
        dv = [(p.lat, p.lon) for p in z.vertices]; ov = vertices_of_int(g, val)
        if len(dv) != len(ov) or _srt(dv) != _srt(ov): v_fail += 1
        # neighbours (EXACT int set)
        onb = neighbour_ints_of_int(g, val)
        if onb is not None:
            nb_checked += 1
            if {n.value for n in z.neighbors} != onb: n_fail += 1
    assert c_fail == 0, f"{c_fail}/500 centroid mismatches"
    assert v_fail == 0, f"{v_fail}/500 vertex mismatches"
    assert nb_checked > 400 and n_fail == 0, f"{n_fail}/{nb_checked} neighbour mismatches"

def test_neighbours_exact_any_seed():
    """A0's exact neighbours match pydggal on a DIFFERENT seed than the run above
    (the old geometric k-ring passed only by fixed-seed luck)."""
    g = oracle_grid("ISEA3H"); rng = random.Random(99991)
    n_fail = 0; nb_checked = 0
    for _ in range(1500):
        lat = rng.uniform(-89.5, 89.5); lon = rng.uniform(-180, 180); res = rng.randint(1, 18)
        val = oracle_zone_int(g, lat, lon, res)
        onb = neighbour_ints_of_int(g, val)
        if onb is None:
            continue
        nb_checked += 1
        if {n.value for n in Zone(ISEA3H, val).neighbors} != onb:
            n_fail += 1
    assert nb_checked > 1200 and n_fail == 0, f"{n_fail}/{nb_checked} neighbour mismatches"
