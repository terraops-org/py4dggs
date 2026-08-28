# tests/test_z7.py
"""Z7Indexing bit-packing: self-consistency tests. The congruent digit-path
scheme (parent = drop last digit, children = append digit 0-6) is correct BY
DEFINITION, not something a second implementation needs to confirm -- there
is no pydggal oracle for it anyway (pydggal's getZoneParent/getZoneChildren
return DGGAL's *geometric* Z7 hierarchy, a different, non-congruent relation
-- see test_isea7h_fuzz.py::test_igeo7_hierarchy_self_consistent). Discrete
correctness vs DGGAL (forward/centroid/vertices/neighbours) is covered live
in test_isea7h_fuzz.py and via frozen golden tables in test_conformance.py."""
import pytest
from py4dggs import IGEO7
from py4dggs.indexings.z7 import Z7Indexing
from py4dggs.types import InvalidZoneError

Z7 = Z7Indexing()


def test_z7_max_resolution_matches_dggal():
    """Audit A1 (root cause): DGGAL supports Z7 levels 0..19 (pydggal
    getMaxDGGRSZoneLevel() == 19); level 20 is representable in the 20-slot
    packing but is nullZone in DGGAL. max_resolution was off by one (20), which
    made zone_from_geo(..., 20) crash in from_7h. It must be 19."""
    assert Z7.max_resolution == 19


def test_zone_from_geo_past_max_resolution_rejected():
    """Audit A1: requesting a resolution above max is a clean InvalidZoneError,
    not a crash (previously res 20 raised IndexError deep in from_7h)."""
    with pytest.raises(InvalidZoneError):
        IGEO7.zone_from_geo(38.7223, -9.1393, Z7.max_resolution + 1)


def test_children_at_max_resolution_is_empty():
    """Audit A1: a max-resolution Z7 zone has no children -- you cannot refine
    past ``max_resolution`` (matches pydggal getZoneChildren count == 0 at level
    19). Grid.children previously crashed with IndexError trying to append a
    digit; the congruent digit-path now returns () at the boundary, mirroring
    parent() -> None at resolution 0. The level just below max still has 7
    (congruent digit hierarchy)."""
    z = IGEO7.zone_from_geo(38.7223, -9.1393, Z7.max_resolution)
    assert z.resolution == Z7.max_resolution
    assert z.children == ()
    z_below = IGEO7.zone_from_geo(38.7223, -9.1393, Z7.max_resolution - 1)
    assert len(z_below.children) == 7


def test_encode_decode_roundtrip():
    base, digits = 3, [6, 5, 3, 1]
    v = Z7.encode(base, digits)
    assert Z7.decode(v) == (base, digits)


def test_resolution_and_pentagon():
    v = Z7.encode(0, [0, 0])                              # centre path => pentagon
    assert Z7.resolution(v) == 2
    assert Z7.is_pentagon(v) is True
    h = Z7.encode(0, [6, 4, 1, 5, 6])                     # hexagon
    assert Z7.is_pentagon(h) is False
    assert Z7.num_children(h) == 7 and Z7.num_children(v) == 6


def test_parent():
    v = Z7.encode(6, [4, 1, 5])
    assert Z7.decode(Z7.parent(v)) == (6, [4, 1])
    assert Z7.parent(Z7.encode(3, [])) is None   # base cell (res 0) has no parent


def test_text_roundtrip():
    v = Z7.encode(0, [6, 4, 1, 5, 6])
    assert Z7.to_text(v) == "0064156"
    assert Z7.from_text("0064156") == v


def test_from_text_rejects_noncanonical_pentagon_child():
    # base 0 is a north pentagon; digit 2 is the deleted child -> invalid
    with pytest.raises(InvalidZoneError):
        Z7.from_text("002")


def test_from_text_rejects_oversized_id():
    with pytest.raises(InvalidZoneError):
        Z7.from_text("00" + "1" * 21)   # 21 digits > the 20-slot Z7 packing


def test_z7_fuzz_self_consistent():
    """Seeded fuzz: encode/decode/text-id round-trip identity across hundreds
    of random ids -- pure integer bit-packing needs no external oracle."""
    import random
    rng = random.Random(20260630)
    for _ in range(500):
        base = rng.randint(0, 11)
        res = rng.randint(0, 12)
        digits = [rng.randint(0, 6) for _ in range(res)]
        v = Z7.encode(base, digits)
        assert Z7.decode(v) == (base, digits)
        assert Z7.resolution(v) == res

        text = Z7.to_text(v)
        try:
            rt = Z7.from_text(text)
            # encode() bypasses the canonical-child check, so some fuzzed ids
            # are deliberately non-canonical pentagon paths; if from_text
            # accepts the text it must round-trip back to the same value.
            assert rt == v
        except InvalidZoneError:
            pass  # non-canonical pentagon path; from_text correctly rejects it
