"""Live differential fuzz: py4dggs.RTEA3H vs pydggal RTEA3H (point-keyed geometry).

Mirrors test_ivea3h_fuzz.py / test_isea3h_fuzz.py. Before A0 the neighbour
assertion here was TOLERANT, because py4dggs.Grid.neighbors used a grid-agnostic
GEOMETRIC edge-crossing k-ring that overshot ~0.01% of aperture-3 cells (at
apices / interruption seams / pentagons / polar rows). A0 (sub-project A, first
slice) gave hex_a3 an exact topological `neighbors` override (a faithful port of
DGGAL's I3HZone::getNeighbor/getNeighbors), so RTEA3H neighbours now match
pydggal EXACTLY for every cell -- hence this file was reverted to STRICT
int-set equality (the acceptance bar Jorge set), matching ISEA3H/IVEA3H.
"""
import random
import pytest
from py4dggs import RTEA3H, Zone
from _pydggal_oracle import (requires_pydggal, oracle_grid, oracle_zone_int,
                             centroid_of_int, vertices_of_int, neighbour_ints_of_int)

pytestmark = requires_pydggal
TOL = 1e-7

def _srt(pairs, nd=4):
    return sorted((round(a, nd), round(b, nd)) for a, b in pairs)

def test_fuzz_geometry_matches_pydggal():
    g = oracle_grid("RTEA3H"); rng = random.Random(20260703)
    c_fail = v_fail = n_fail = 0; nb_checked = 0
    for _ in range(500):
        lat = rng.uniform(-88, 88); lon = rng.uniform(-180, 180); res = rng.randint(1, 11)
        val = oracle_zone_int(g, lat, lon, res)
        z = Zone(RTEA3H, val)                    # DGGAL-layout int == our value
        # centroid (exact -- no tolerance)
        olat, olon = centroid_of_int(g, val); c = z.centroid
        if abs(c.lat - olat) > TOL or abs(c.lon - olon) > TOL: c_fail += 1
        # vertices (exact -- no tolerance; set + count)
        dv = [(p.lat, p.lon) for p in z.vertices]; ov = vertices_of_int(g, val)
        if len(dv) != len(ov) or _srt(dv) != _srt(ov): v_fail += 1
        # neighbours (EXACT int set -- strict since A0; see module docstring)
        onb = neighbour_ints_of_int(g, val)
        if onb is not None:
            nb_checked += 1
            if {n.value for n in z.neighbors} != onb: n_fail += 1
    assert c_fail == 0, f"{c_fail}/500 centroid mismatches"
    assert v_fail == 0, f"{v_fail}/500 vertex mismatches"
    assert nb_checked > 400 and n_fail == 0, f"{n_fail}/{nb_checked} neighbour mismatches"

def test_neighbours_exact_any_seed():
    """A0's exact neighbours match pydggal on a DIFFERENT seed than the run above."""
    g = oracle_grid("RTEA3H"); rng = random.Random(77771)
    n_fail = 0; nb_checked = 0
    for _ in range(1500):
        lat = rng.uniform(-89.5, 89.5); lon = rng.uniform(-180, 180); res = rng.randint(1, 18)
        val = oracle_zone_int(g, lat, lon, res)
        onb = neighbour_ints_of_int(g, val)
        if onb is None:
            continue
        nb_checked += 1
        if {n.value for n in Zone(RTEA3H, val).neighbors} != onb:
            n_fail += 1
    assert nb_checked > 1200 and n_fail == 0, f"{n_fail}/{nb_checked} neighbour mismatches"
