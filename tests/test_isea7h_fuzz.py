"""Task: live differential fuzz -- py4dggs.IGEO7 vs the pydggal ISEA7H_Z7 oracle.

IGEO7 is the flagship grid (ISEA projection + aperture-7 hex topology + Z7
indexing, the canonical 11.2 deg / authalic config) -- historically it was
verified against a frozen Python port (`igeo7-py`, a separate, single-grid
sibling repo) instead of pydggal directly. That sibling-repo dev dependency
was removed (2026-07-07): pydggal is now the sole oracle for every grid in
this library, IGEO7 included, matching the convention already used by
IVEA7H/RTEA7H/ISEA3H/IVEA3H/RTEA3H (see test_rtea7h_fuzz.py, established
2026-07-01). igeo7-py's own real vertex/neighbour bug (dropped eC i1/i2
interruption-frame reconciliation) is why vertices/neighbours were always
checked against pydggal, not igeo7-py, even before this change -- see
test_igeo7_vertices_pydggal.py for the dedicated interruption-cell sweep;
this file's vertices/neighbours checks are the same broad-random-fuzz
treatment given to the other 7H grids, for consistency.

Hierarchy (parent/children) is congruent for Z7 (drop/append one base-7
digit), so it is checked for an EXACT match, not a tolerance -- unlike
centroid/vertices, which use pydggal's C-engine float agreement (~1e-9 deg,
not bit-identical).
"""
import random

import pytest
from _pydggal_oracle import requires_pydggal, oracle_grid, geopoint, centroid_of, vertices_of, neighbours_of
from py4dggs import IGEO7

TOL = 1e-9  # degrees; comfortably above the observed ~1e-12 py4dggs-vs-pydggal gap


def _degenerate(tid):
    z = IGEO7.zone_from_text(tid)
    return any((p.lat, p.lon) == (0.0, 0.0) for p in z.vertices)


@requires_pydggal
def test_igeo7_centroid_matches_pydggal():
    """The strong correctness gate: cell centroid within TOL of pydggal."""
    g = oracle_grid("ISEA7H_Z7")
    rng = random.Random(20260707)
    n = 0
    for _ in range(500):
        lat = rng.uniform(-89.9, 89.9); lon = rng.uniform(-180.0, 180.0); res = rng.randint(1, 12)
        tid = g.getZoneTextID(g.getZoneFromWGS84Centroid(res, geopoint(lat, lon)))
        olat, olon = centroid_of(g, tid)
        c = IGEO7.zone_from_text(tid).centroid
        assert abs(c.lat - olat) < TOL and abs(c.lon - olon) < TOL, f"centroid drift at {tid}"
        n += 1
    assert n == 500


@requires_pydggal
def test_igeo7_forward_matches_pydggal():
    """forward cell id matches pydggal except rare boundary tie-breaks (>=99%)."""
    g = oracle_grid("ISEA7H_Z7")
    rng = random.Random(20260707)
    n = ok = 0
    for _ in range(1000):
        lat = rng.uniform(-89.9, 89.9); lon = rng.uniform(-180.0, 180.0); res = rng.randint(1, 12)
        want = g.getZoneTextID(g.getZoneFromWGS84Centroid(res, geopoint(lat, lon)))
        got = IGEO7.zone_from_geo(lat, lon, res).text_id
        n += 1; ok += (got == want)
    assert ok / n >= 0.99, f"IGEO7 forward match rate {ok}/{n} below 99%"


@requires_pydggal
def test_igeo7_vertices_match_pydggal():
    """Vertices match pydggal as an unordered set within TOL, skipping degenerate
    (0,0)-vertex cells (see test_igeo7_vertices_pydggal.py for the dedicated
    interruption-cell correctness sweep)."""
    g = oracle_grid("ISEA7H_Z7")
    rng = random.Random(20260707)
    checked = skipped = 0
    for _ in range(500):
        lat = rng.uniform(-89.9, 89.9); lon = rng.uniform(-180.0, 180.0); res = rng.randint(1, 12)
        tid = g.getZoneTextID(g.getZoneFromWGS84Centroid(res, geopoint(lat, lon)))
        if _degenerate(tid):
            skipped += 1; continue
        ov = vertices_of(g, tid)
        dv = [(p.lat, p.lon) for p in IGEO7.zone_from_text(tid).vertices]
        assert len(dv) == len(ov), f"vertex count at {tid}"
        for d in dv:
            assert min(abs(d[0] - o[0]) + abs(d[1] - o[1]) for o in ov) < TOL, f"vertex drift at {tid}"
        checked += 1
    assert checked > 400, f"too few non-degenerate cells checked ({checked}); skipped {skipped}"


@requires_pydggal
def test_igeo7_neighbours_match_pydggal():
    """Neighbour text-id set matches pydggal, skipping degenerate cells and any
    zone the pydggal binding overflows on."""
    g = oracle_grid("ISEA7H_Z7")
    rng = random.Random(20260707)
    checked = 0
    for _ in range(500):
        lat = rng.uniform(-89.9, 89.9); lon = rng.uniform(-180.0, 180.0); res = rng.randint(1, 12)
        tid = g.getZoneTextID(g.getZoneFromWGS84Centroid(res, geopoint(lat, lon)))
        if _degenerate(tid):
            continue
        want = neighbours_of(g, tid)
        got = {z.text_id for z in IGEO7.zone_from_text(tid).neighbors}
        assert got == want, f"k-ring drift at {tid}: {got} vs {want}"
        checked += 1
    # Historically ~45% of Z7 zones (base cells 8-11, bit 63 set) were skipped
    # here: pydggal 0.0.6's getZoneNeighbors wrapper overflows on them. The
    # oracle now goes around that at the C-ABI layer (_pydggal_oracle._zone_array),
    # so NO cell is skipped for binding reasons and base cells 8-11 are verified
    # like any other. The threshold is deliberately just under the non-degenerate
    # count: if guard-and-skip ever creeps back, coverage drops by ~45% and this
    # fails instead of silently shrinking.
    assert checked > 450


def test_igeo7_hierarchy_self_consistent():
    """Z7 hierarchy is congruent BY DEFINITION (parent = drop last direction
    digit, children = append digit 0-6) -- there is no pydggal oracle for it:
    pydggal's own getZoneParent/getZoneChildren return DGGAL's *geometric*
    Z7 hierarchy (non-congruent, 2 parents / 13 children for interior cells --
    see documentation/dggs-py-port-lessons.md finding #21 in the sibling
    ut.IGEO7 repo), a genuinely different relation. igeo7-spec's own golden
    hierarchy table generator says so explicitly ("getZoneChildren/
    getZoneParents -> NON-CONGRUENT geometric hierarchy -- NOT used"; it
    builds parent/children by string truncation/append instead, exactly what
    py4dggs does). So the real correctness gate is the frozen golden-table
    conformance (test_conformance.py::test_hierarchy) plus this: every child
    py4dggs produces must round-trip back to its exact parent (self-consistency,
    no oracle required)."""
    rng = random.Random(20260707)
    for _ in range(300):
        lat = rng.uniform(-89.9, 89.9); lon = rng.uniform(-180.0, 180.0); res = rng.randint(1, 11)
        z = IGEO7.zone_from_geo(lat, lon, res)
        for child in z.children:
            assert child.parent.text_id == z.text_id, (
                f"{child.text_id}: parent {child.parent.text_id!r} != {z.text_id!r}"
            )
        assert len(z.children) == (6 if z.is_pentagon else 7)


@requires_pydggal
def test_igeo7_pentagon_kring_matches_pydggal():
    """Directly exercise known non-polar pentagons at several resolutions (the
    risky path for the 5-edge k-ring); children are checked structurally, not
    against pydggal (see test_igeo7_hierarchy_self_consistent for why)."""
    g = oracle_grid("ISEA7H_Z7")
    for tid in ("060", "070", "0100", "01000", "010000"):
        z = IGEO7.zone_from_text(tid)
        assert z.is_pentagon, f"{tid}: should be a pentagon"

        zneigh = {n.text_id for n in z.neighbors}
        oneigh = neighbours_of(g, tid)
        if oneigh is not None:
            assert zneigh == oneigh, f"{tid}: neighbors {zneigh} != oracle {oneigh}"
        assert len(zneigh) == 5, f"{tid}: pentagon should have 5 neighbors, got {len(zneigh)}"

        zch = {c.text_id for c in z.children}
        assert all(c.parent.text_id == tid for c in z.children), f"{tid}: child parent round-trip"
        assert len(zch) == 6, f"{tid}: pentagon should have 6 children, got {len(zch)}"


@requires_pydggal
def test_igeo7_polar_pentagons_match_pydggal():
    """Base-0 (north pole) and base-11 (south pole) pentagons exercise the
    engine's special polar branches."""
    g = oracle_grid("ISEA7H_Z7")
    for tid in ("000", "0000", "11000", "110000"):
        z = IGEO7.zone_from_text(tid)
        assert z.is_pentagon, f"{tid}: must be a pentagon"

        got_nb = {n.text_id for n in z.neighbors}
        want_nb = neighbours_of(g, tid)
        if want_nb is not None:
            assert got_nb == want_nb, f"{tid}: neighbor mismatch"
        assert len(got_nb) == 5, f"{tid}: pentagon must have 5 neighbors"

        assert all(c.parent.text_id == tid for c in z.children), f"{tid}: child parent round-trip"
        assert len(z.children) == 6, f"{tid}: pentagon must have 6 children"


def test_igeo7_zone_is_immutable_and_hashable():
    """No oracle needed -- a plain structural property of Zone."""
    z = IGEO7.zone_from_geo(38.7223, -9.1393, 5)
    assert z in {z}
    with pytest.raises(AttributeError):
        z.foo = 1
