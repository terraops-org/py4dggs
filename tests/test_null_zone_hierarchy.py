"""I3H NULL_ZONE hierarchy mirrors pydggal exactly.

Policy (2026-07-29, extended to hierarchy 2026-08-21): for nullZone cases py4dggs
mirrors pydggal rather than inventing its own error signalling. `centroid` /
`text_id` were covered then; `parents()` / `children()` were not, and they
diverged -- returning ordinary-looking zones manufactured from the "no such zone"
sentinel, including a *resolution-0* zone from `pack_i3h(31 + 1, 15, 0, 0)`.

The generic path cannot reproduce DGGAL's answers here: NULL_ZONE is all-ones,
and bits 62-63 lie OUTSIDE all four bitfields (levelI9R:5:57, rootRhombus:4:53,
rhombusIX:51:2, subHex:2:0), so any value rebuilt from an unpacked tuple loses
them. The mirror is therefore a lookup, pinned against the live oracle here so it
cannot drift silently.

NULL_ZONE is genuinely reachable: `zone_from_geo` within ~0.1 deg of a pole at the
coarsest resolutions hits it about 1 in 10^4 near-pole points.
"""
import pytest

from py4dggs import ISEA3H
from py4dggs.types import NULL_ZONE

pytest.importorskip("dggal")


def _oracle():
    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from _pydggal_oracle import oracle_grid, _zone_array
    return oracle_grid("ISEA3H"), _zone_array


def test_null_zone_parents_match_pydggal_exactly():
    g, _zone_array = _oracle()
    assert list(ISEA3H.parents(NULL_ZONE)) == _zone_array("DGGRS_getZoneParents", g, NULL_ZONE)


def test_null_zone_children_match_pydggal_exactly():
    """Mirrored verbatim, including DGGAL's duplicate entries and the fact that
    NULL_ZONE is listed as its own child -- 'mirror exactly' means exactly."""
    g, _zone_array = _oracle()
    assert list(ISEA3H.children(NULL_ZONE)) == _zone_array("DGGRS_getZoneChildren", g, NULL_ZONE)


def test_null_zone_children_no_longer_fabricate_a_resolution_zero_zone():
    """The specific regression: pack_i3h(31 + 1, 15, 0, 0) unpacked as (0, 15, 0, 0)."""
    assert 4746794007248502784 not in ISEA3H.children(NULL_ZONE)


def test_ordinary_zone_hierarchy_is_untouched():
    g, _zone_array = _oracle()
    v = ISEA3H.zone_from_geo(52.0, 5.0, 6).value
    assert sorted(ISEA3H.parents(v)) == sorted(_zone_array("DGGRS_getZoneParents", g, v))
    assert sorted(ISEA3H.children(v)) == sorted(_zone_array("DGGRS_getZoneChildren", g, v))
