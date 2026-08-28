"""IGEO7 vertex/neighbour correctness vs pydggal (DGGAL ground truth).

The Phase-0 oracle (igeo7-py) carries a real bug: interruption-spanning cells get
degenerate (0,0) / wrong "far-side" boundary vertices (the eC's i1/i2 frame
reconciliation in addNonPolarBaseVertices, RI7H.ec:1991-2004, was dropped in the
JS/Python port). Because py4dggs was bit-identical to igeo7-py, it inherited the bug.

These tests verify IGEO7 directly against pydggal — the DGGAL engine itself, which
is correct — for the interruption-spanning cells (where igeo7-py is wrong) and a
broad sweep. Vertices/neighbours are the ONLY outputs where py4dggs deliberately
diverges from igeo7-py; text_id / centroid / parent / children stay bit-identical
to igeo7-py (see test_isea7h_fuzz.py, which superseded the old igeo7-py-dependent
test_grid_igeo7.py -- pydggal is now the sole oracle for every grid in this repo).
"""
import pytest

from py4dggs import IGEO7
from _pydggal_oracle import requires_pydggal, oracle_grid, vertices_of, neighbours_of

# NOTE: requires_pydggal is applied per-test to the oracle-comparing tests only.
# test_no_degenerate_vertices_anywhere uses ONLY py4dggs (no pydggal) and must run in
# every CI — it is the primary always-on guard against a regression that
# re-introduces the interruption bug (degenerate vertices) when dggal is absent.


def _sorted_round(pairs, nd=4):
    return sorted((round(a, nd), round(b, nd)) for a, b in pairs)


def _sweep_zones(resolutions=(2, 3, 4), step=4.0):
    """Distinct IGEO7 zones over a dense lat/lon sweep (hits interruption cells)."""
    seen = {}
    lat = -89.0
    while lat <= 89.0:
        lon = -180.0
        while lon < 180.0:
            for res in resolutions:
                z = IGEO7.zone_from_geo(lat, lon, res)
                seen.setdefault(z.text_id, z)
            lon += step
        lat += step
    return seen


# Cells confirmed (2026-07-02) to carry the igeo7-py degeneracy bug: their far-side
# vertices came out (0,0) or on the wrong face. These MUST now match pydggal.
KNOWN_INTERRUPTION_CELLS = [
    "0604", "06004", "06040", "060400", "06044", "0641", "06410", "06414",
    "064145", "0645", "06454", "064545", "110002", "1102", "11025", "1124",
]


def _oracle():
    return oracle_grid("ISEA7H_Z7")


@requires_pydggal
@pytest.mark.parametrize("tid", KNOWN_INTERRUPTION_CELLS)
def test_known_interruption_cell_vertices_match_pydggal(tid):
    g = _oracle()
    py4dggs_v = [(v.lat, v.lon) for v in IGEO7.zone_from_text(tid).vertices]
    # no degenerate (0,0) vertices
    assert not any(abs(la) < 1e-12 and abs(lo) < 1e-12 for la, lo in py4dggs_v), \
        f"{tid}: degenerate (0,0) vertex in {py4dggs_v}"
    assert _sorted_round(py4dggs_v) == _sorted_round(vertices_of(g, tid)), \
        f"{tid}: py4dggs vertices != pydggal"


def test_no_degenerate_vertices_anywhere():
    """Across a broad sweep, NO IGEO7 cell may have a (0,0) degenerate vertex."""
    zones = _sweep_zones()
    bad = []
    for tid, z in zones.items():
        if any(abs(v.lat) < 1e-12 and abs(v.lon) < 1e-12 for v in z.vertices):
            bad.append(tid)
    assert not bad, f"{len(bad)}/{len(zones)} cells still degenerate: {sorted(bad)[:20]}"


@requires_pydggal
def test_sweep_vertices_match_pydggal():
    """Every swept IGEO7 cell's vertices match pydggal (order-insensitive, ~11 m)."""
    g = _oracle()
    zones = _sweep_zones()
    fails = []
    for tid, z in sorted(zones.items()):
        py4dggs_v = _sorted_round([(v.lat, v.lon) for v in z.vertices])
        ora_v = _sorted_round(vertices_of(g, tid))
        if py4dggs_v != ora_v:
            fails.append(tid)
    assert not fails, f"{len(fails)}/{len(zones)} cells' vertices != pydggal: {fails[:20]}"


@requires_pydggal
def test_known_interruption_cell_neighbours_match_pydggal():
    """The k-ring cascade: fixing vertices must fix neighbours for these cells too."""
    g = _oracle()
    checked = 0
    for tid in KNOWN_INTERRUPTION_CELLS:
        ora = neighbours_of(g, tid)
        if ora is None:  # pydggal int64 overflow for high-bit zones — skip
            continue
        got = {n.text_id for n in IGEO7.zone_from_text(tid).neighbors}
        assert got == ora, f"{tid}: py4dggs neighbours {got} != pydggal {ora}"
        checked += 1
    assert checked > 0, "no interruption cell neighbours could be checked"
