"""A1 rhombic text-id (getZoneID/fromZoneID) for the I3H family."""
import pytest
from py4dggs.indexings.i3h import I3HIndexing, pack_i3h
from py4dggs.types import InvalidZoneError

IDX = I3HIndexing()

# (levelI9R, root, rhombusIX, subHex) -> canonical text-id (verified vs pydggal)
PAIRS = [
    (pack_i3h(2, 2, 35, 2), "C2-23-C"),
    (pack_i3h(2, 4, 8, 3),  "C4-8-D"),
    (pack_i3h(2, 4, 30, 0), "C4-1E-A"),
    (pack_i3h(0, 4, 0, 1),  "A4-0-B"),
    (pack_i3h(1, 0, 1, 2),  "B0-1-C"),
    (pack_i3h(2, 10, 0, 0), "CA-0-A"),   # "North" pole (root 10 -> 'A')
    (pack_i3h(2, 11, 0, 0), "CB-0-A"),   # "South" pole (root 11 -> 'B')
]


@pytest.mark.parametrize("value,text", PAIRS)
def test_to_text(value, text):
    assert IDX.to_text(value) == text


@pytest.mark.parametrize("value,text", PAIRS)
def test_from_text_roundtrip(value, text):
    assert IDX.from_text(text) == value


@pytest.mark.parametrize("bad", [
    "not-an-id", "C2-23-Z",   # malformed / bad sub-hex
    "ZZ-0-A",                 # malformed: 'Z' is not a hex digit for the ix field
    "R0-0-A",                 # out-of-range level (R = levelI9R 17 > 16, iLRCFromLRtI rejects)
    "CC-0-A",                 # root 12: parses via the polar branch but validate rejects (root > 11)
    "c4-8-D", "C04-8-D",      # non-canonical spelling (lowercase / leading zero)
    "C2--A", "",              # malformed
])
def test_from_text_rejects(bad):
    # Audit A2: I3H must raise InvalidZoneError (as Z7 does), not a bare
    # ValueError, so callers can `except InvalidZoneError` uniformly across
    # every grid. InvalidZoneError subclasses ValueError, so this is stricter.
    with pytest.raises(InvalidZoneError):
        IDX.from_text(bad)


# --- differential fuzz vs pydggal (skips cleanly if dggal absent) ---------- #
import random
from py4dggs import ISEA3H, IVEA3H, RTEA3H
from _pydggal_oracle import requires_pydggal, oracle_grid, geopoint


@requires_pydggal
@pytest.mark.parametrize("name,grid", [("ISEA3H", ISEA3H), ("IVEA3H", IVEA3H), ("RTEA3H", RTEA3H)])
def test_textid_exact_vs_pydggal(name, grid):
    g = oracle_grid(name)
    rng = random.Random(20260703)
    enc_fail = dec_fail = 0
    for _ in range(3000):
        lat = rng.uniform(-89.9, 89.9); lon = rng.uniform(-180, 180); res = rng.randint(0, 33)
        z = g.getZoneFromWGS84Centroid(res, geopoint(lat, lon)); zi = int(z)
        ptid = g.getZoneTextID(z)
        if grid.indexing.to_text(zi) != ptid:
            enc_fail += 1
        if grid.indexing.from_text(ptid) != zi:
            dec_fail += 1
    assert enc_fail == 0 and dec_fail == 0, f"{name}: {enc_fail} enc / {dec_fail} dec"


@requires_pydggal
@pytest.mark.parametrize("name,grid", [("ISEA3H", ISEA3H), ("IVEA3H", IVEA3H), ("RTEA3H", RTEA3H)])
def test_pole_textids_roundtrip(name, grid):
    for l9r in range(0, 17):
        for pole in (10, 11):
            for sh in (0, 1):
                val = pack_i3h(l9r, pole, 0, sh)
                assert grid.indexing.from_text(grid.indexing.to_text(val)) == val
