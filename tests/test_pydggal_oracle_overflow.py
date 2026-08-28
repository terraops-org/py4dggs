"""The pydggal binding overflows on ~45% of Z7 zones -- the oracle must not skip them.

`dggal==0.0.6`'s high-level `Array` copy marshals a C `uint64` DGGRSZone through a
*signed* int64 (``ecrt.py TA(): u.i64 = a``), so ``getZoneNeighbors`` /
``getZoneParents`` / ``getZoneChildren`` raise ``OverflowError`` ("int too big to
convert") for every Z7 zone with bit 63 set -- base cells 8-11, ~45% of zones
(measured 134/300 at res 5). This is a *wrapper* bug: ``libdggal.so`` itself fills a
correct C uint64 array, so the values are recoverable through the binding's own cffi
connector (``dggal.lib`` + ``dggal.ffi``). See findings #11 and #22 in
``documentation/dggs-py-port-lessons.md`` in the sibling ``ut.IGEO7`` repo.

Until now ``_pydggal_oracle`` guard-and-skipped those zones (``except OverflowError:
return None``) and the fuzz tests dropped them (``if want is None: continue``), so
base cells 8-11 -- "exactly where prior bugs hid", per that same catalog -- were never
verified against the oracle at all. These tests pin the recovery.
"""
import pytest

from _pydggal_oracle import (
    requires_pydggal,
    oracle_grid,
    neighbour_ints_of_int,
    parents_of_int,
    children_of_int,
)

# A res-5 ISEA7H_Z7 zone in the base-cell 8-11 range (bit 63 set), i.e. one the
# 0.0.6 wrapper cannot marshal. Frozen as a literal so the test pins the exact
# regression rather than depending on a fuzz draw.
OVERFLOWING_ZONE = 9830935774084726783


@requires_pydggal
def test_wrapper_still_overflows_on_high_bit_zones():
    """Pin the upstream bug itself: if a future dggal fixes it, this fails loudly
    and the cffi bypass can be reconsidered (0.0.6 is the latest release as of
    2026-08-21, so there is nothing to upgrade to yet)."""
    g = oracle_grid("ISEA7H_Z7")
    assert OVERFLOWING_ZONE >> 63 == 1, "test fixture must have bit 63 set"
    with pytest.raises(OverflowError):
        g.getZoneNeighbors(OVERFLOWING_ZONE)


@requires_pydggal
def test_oracle_recovers_neighbours_the_wrapper_cannot_marshal():
    g = oracle_grid("ISEA7H_Z7")
    got = neighbour_ints_of_int(g, OVERFLOWING_ZONE)
    assert got is not None, "oracle still guard-and-skips zones the wrapper overflows on"
    assert len(got) == 6
    assert all(v >> 63 in (0, 1) for v in got)


@requires_pydggal
def test_oracle_recovers_hierarchy_the_wrapper_cannot_marshal():
    g = oracle_grid("ISEA7H_Z7")
    kids = children_of_int(g, OVERFLOWING_ZONE)
    rents = parents_of_int(g, OVERFLOWING_ZONE)
    assert kids, "children oracle still blocked by the binding overflow"
    assert rents, "parents oracle still blocked by the binding overflow"
