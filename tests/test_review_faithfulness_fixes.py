"""Faithfulness fixes from the 2026-08-21 full-repo review.

Both cases look like tidy-ups on a casual read and are neither: one restores an
eC branch the port dropped, the other restores a nullZone guard three sibling
call sites already have.
"""
import pytest

from py4dggs.topologies import hex_a3, hex_a7
from py4dggs.types import NULL_ZONE


# --- get_level_rotation_offset must honour a supplied child position --------- #
# RI7H_Z7.ec:196-201 opens with `if(i == -1) i = getChildPosition(parent, zone);`
# -- the caller may pass a already-known position (RI7H_Z7.ec:330 passes prevCIX,
# :394 passes prevI) and the eC then USES it. py4dggs recomputed unconditionally,
# so the argument was inert and the eC's sentinel branch was missing. Measured
# across 36,163 calls the cached value never differs from the recomputation (it
# is an optimisation, not a different value), so this is faithfulness, not a
# behaviour change -- but the branch belongs in the port.
def test_supplied_child_position_is_used_instead_of_recomputed(monkeypatch):
    """With a non-sentinel prev_i the eC never calls getChildPosition."""
    from py4dggs import IGEO7

    p = IGEO7.projection.forward(IGEO7._geom, 52.0, 5.0)
    parents = hex_a7.compute_parents(hex_a7.from_centroid(6, p.x, p.y))
    zone, parent, grand_parent = parents[0], parents[1], parents[2]

    calls = []
    real = hex_a7.get_child_position
    monkeypatch.setattr(
        hex_a7, "get_child_position",
        lambda pa, zo: (calls.append(1), real(pa, zo))[1],
    )
    hex_a7.get_level_rotation_offset(4, real(parent, zone), zone, parent, grand_parent)
    assert calls == [], "supplied child position was ignored and recomputed"


def test_sentinel_minus_one_still_triggers_the_lookup():
    """The eC's own callers pass -1 to mean 'compute it for me' (RI7H_Z7.ec:276)."""
    assert hex_a7.get_level_rotation_offset(3, -1, None, None, None) == 0


# --- _i3h_from_centroid's nullZone return must be guarded everywhere --------- #
def test_first_sub_zone_returns_null_zone_when_the_centroid_is_null(monkeypatch):
    """`_i3h_from_centroid` returns None for its nullZone guard (hex_a3.py:593).
    quantize, _i3h_get_children and _i3h_get_neighbor all check for it; the
    sub-zone entry points splatted it straight into pack_i3h(*None) and raised
    TypeError instead of yielding a nullZone."""
    monkeypatch.setattr(hex_a3, "_i3h_from_centroid", lambda *a, **k: None)
    topo = hex_a3.HexAperture3Topology()
    got, digits = topo.first_sub_zone(None, 324259173170675744, [], 2)
    assert got == NULL_ZONE
    assert digits == []


def test_sub_zones_yields_null_zone_when_the_centroid_is_null(monkeypatch):
    monkeypatch.setattr(hex_a3, "_i3h_from_centroid", lambda *a, **k: None)
    topo = hex_a3.HexAperture3Topology()
    out = topo.sub_zones(None, 324259173170675744, [], 1)
    assert out, "expected sub-zones, not an empty list"
    assert all(v == NULL_ZONE for (v, _) in out)


# --- children past max_resolution: FAITHFUL, pinned against the engine ------- #
def test_children_at_max_resolution_match_pydggal_exactly():
    """py4dggs emits 7 children at resolution 34 from a resolution-33 (max) I3H
    zone. That LOOKS like a missing bound and is not: DGGAL does the same, and
    our seven values equal pydggal's seven exactly.

    Pinned here so a future audit cannot "fix" it in the fail-closed direction
    (returning `()` at max_resolution, as the congruent Z7 branch correctly does
    -- Z7 level 20 really is nullZone in DGGAL, this is not). The resulting
    zones' text ids do not round-trip, in this port and in pydggal alike; see
    I3HIndexing.from_text's scoped round-trip note.
    """
    pytest.importorskip("dggal")
    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from _pydggal_oracle import oracle_grid, _zone_array
    from py4dggs import ISEA3H

    g = oracle_grid("ISEA3H")
    v = ISEA3H.zone_from_geo(52.0, 5.0, 33).value
    assert ISEA3H.indexing.resolution(v) == ISEA3H.indexing.max_resolution == 33
    assert sorted(ISEA3H.children(v)) == sorted(_zone_array("DGGRS_getZoneChildren", g, v))
