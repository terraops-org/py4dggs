import pytest
from py4dggs.indexings.i3h import I3HIndexing, pack_i3h, unpack_i3h

IDX = I3HIndexing()

def test_pack_unpack_roundtrip():
    for li9r, root, rix, sh in [(0,0,0,0),(2,2,35,2),(15,11,(3**15)**2 - 1,3),(5,10,0,1)]:
        v = pack_i3h(li9r, root, rix, sh)
        assert unpack_i3h(v) == (li9r, root, rix, sh)

def test_matches_dggal_layout():
    # Lisbon res5 zone int from pydggal (verified 2026-07-02) decodes to the fields.
    assert unpack_i3h(306244774661193870) == (2, 2, 35, 2)
    assert pack_i3h(2, 2, 35, 2) == 306244774661193870

def test_derived_fields():
    v = pack_i3h(2, 2, 35, 2)             # level 5, hexagon
    assert IDX.resolution(v) == 5
    assert IDX.base_cell(v) == 2
    assert not IDX.is_pentagon(v)
    assert IDX.is_pentagon(pack_i3h(3, 4, 0, 1))    # rix==0, subHex<=1 -> pentagon
    assert IDX.is_pentagon(pack_i3h(3, 10, 0, 0))   # north pole
    assert not IDX.is_pentagon(pack_i3h(3, 4, 0, 2))  # rix==0 but subHex==2 -> hex

def test_encode_decode_identity():
    v = pack_i3h(2, 2, 35, 2)
    assert IDX.encode(v, []) == v
    assert IDX.decode(v) == (v, [])

def test_canonical_text_roundtrip():
    v = pack_i3h(2, 2, 35, 2)
    t = IDX.to_text(v)
    assert t == "C2-23-C"
    assert IDX.from_text(t) == v

def test_hierarchy_not_on_indexing():
    # A2: the I3H hierarchy is geometric (topology/Grid-level), NOT a digit-path on
    # Indexing -- so I3HIndexing intentionally has no parent/num_children/child_digits.
    for m in ("parent", "num_children", "child_digits"):
        assert not hasattr(IDX, m)
