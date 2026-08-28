"""Bounds and identity fixes from the 2026-08-21 full-repo review.

Each test encodes one reviewed failure scenario. Grouped here rather than spread
across the topical files because they share a theme: values that are *reachable*
but unbounded, and comparisons that drop the grid they belong to.
"""
import pytest

from py4dggs import IGEO7, ISEA3H, IVEA3H
from py4dggs.indexings.i3h import pack_i3h, unpack_i3h
from py4dggs.types import InvalidZoneError


# --- pack_i3h: eC bitfield semantics ---------------------------------------- #
# RI3H.ec:855-861 declares I3HZone as bitfields (`uint levelI9R:5:57`, etc.), and
# C bitfields TRUNCATE per field on assignment -- writing 16 to a 4-bit field
# stores 0, it does not bleed into the neighbouring field. py4dggs masked only
# `rix`, so an over-range level/root/sub_hex corrupted an ADJACENT field and
# fabricated a plausible-looking zone.
def test_pack_i3h_truncates_root_instead_of_bleeding_into_level():
    assert unpack_i3h(pack_i3h(0, 16, 0, 0)) == (0, 0, 0, 0)


def test_pack_i3h_truncates_sub_hex_instead_of_bleeding_into_rix():
    assert unpack_i3h(pack_i3h(0, 0, 0, 4)) == (0, 0, 0, 0)


def test_pack_i3h_truncates_level_without_disturbing_root():
    assert unpack_i3h(pack_i3h(32, 4, 0, 0)) == (0, 4, 0, 0)


# --- sub_zones cardinality bound (the DoS) ---------------------------------- #
def test_sub_zones_refuses_a_depth_whose_result_cannot_be_materialised():
    """count_sub_zones is a cheap closed form; sub_zones builds a tuple. Depth 33
    on a res-0 I3H cell counts 4,632,550,579,746,406 sub-zones and passed every
    guard, so the call was an unbounded allocation, not merely a slow one."""
    v = ISEA3H.zone_from_geo(52.0, 5.0, 0).value
    assert ISEA3H.count_sub_zones(v, 33) > 10**15  # still answerable: closed form
    with pytest.raises(InvalidZoneError, match="sub-zones"):
        ISEA3H.sub_zones(v, 33)


def test_sub_zones_still_allows_ordinary_depths():
    v = ISEA3H.zone_from_geo(52.0, 5.0, 4).value
    assert len(ISEA3H.sub_zones(v, 2)) == ISEA3H.count_sub_zones(v, 2)


# --- sub_zone_index depth 0 -------------------------------------------------- #
def test_sub_zone_index_of_a_zone_with_itself_is_zero():
    """sub_zones(v, 0) == (v,), count_sub_zones(v, 0) == 1 and
    sub_zone_at_index(v, 0, 0) == v, so index(v, v) must be 0, not -1."""
    v = ISEA3H.zone_from_geo(52.0, 5.0, 4).value
    assert ISEA3H.sub_zones(v, 0) == (v,)
    assert ISEA3H.sub_zone_index(v, v) == 0


# --- Zone.sub_zone_index must not cross grids -------------------------------- #
def test_sub_zone_index_rejects_a_zone_from_another_grid():
    """ISEA3H/IVEA3H/RTEA3H share an identical I3H packing, so a raw int compare
    cannot tell them apart -- Zone.__eq__ checks grid identity, this must too."""
    a = ISEA3H.zone_from_geo(52.0, 5.0, 4)
    native = a.sub_zones(2)[3]
    foreign = IVEA3H.zone_from_text(native.text_id)
    assert a.sub_zone_index(native) == 3
    assert a.sub_zone_index(foreign) == -1


def test_sub_zone_index_rejects_a_zone_from_a_different_grid_family():
    a = ISEA3H.zone_from_geo(52.0, 5.0, 4)
    assert a.sub_zone_index(IGEO7.zone_from_geo(52.0, 5.0, 6)) == -1
