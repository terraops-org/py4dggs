"""DGGAL's `nullZone` — addresses that are representable but have no geometry.

Two independent cases, both verified against pydggal:

  * **Z7 level 20** (`hex_a7`): the 64-bit packing has 20 direction-digit slots,
    so level 20 is *representable*, but `to7H` (RI7H_Z7.ec:348-353) "does not
    support level 20 zones" and yields nullZone -> centroid (0,0), six (0,0)
    vertices. Before `Topology.is_null_geometry` existed, this port computed a
    plausible-looking real coordinate for such a zone.

  * **I3H NULL_ZONE** (`hex_a3`/`i3h`): `zone_from_geo` genuinely returns the
    sentinel for ordinary lat/lon very near a pole at the coarsest resolutions.
    DGGAL prints it as "(null)"; this port used to stringify the sentinel's
    all-ones bit fields into a normal-looking id, losing DGGAL's only signal
    that the zone does not exist.
"""
import pytest

from py4dggs import IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H, RTEA3H, Zone
from py4dggs.interfaces import Indexing, Topology
from py4dggs.indexings.i3h import I3HIndexing
from py4dggs.indexings.z7 import Z7Indexing
from py4dggs.topologies.hex_a3 import HexAperture3Topology
from py4dggs.topologies.hex_a7 import HexAperture7Topology
from py4dggs.types import NULL_TEXT, NULL_ZONE

from _pydggal_oracle import oracle_grid, requires_pydggal

Z7_GRIDS = [IGEO7, IVEA7H, RTEA7H]
I3H_GRIDS = [ISEA3H, IVEA3H, RTEA3H]

# A level-20 (== max_resolution + 1) Z7 text id, hexagon- and pentagon-path.
LVL20_HEX = "05" + "3" * 20
LVL20_PENT = "05" + "0" * 20

# An ordinary lat/lon that DGGAL itself quantizes to nullZone on the I3H grids.
NULL_POINT = (89.99, -169.0, 0)


# --------------------------------------------------------------------------- #
# Z7 level 20: representable, but no geometry
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("grid", Z7_GRIDS, ids=lambda g: g.name)
@pytest.mark.parametrize("text", [LVL20_HEX, LVL20_PENT], ids=["hex-path", "pent-path"])
def test_z7_level20_has_null_geometry(grid, text):
    """Level 20 degenerates to DGGAL's zero geometry rather than inventing a point."""
    z = grid.zone_from_text(text)
    assert z.resolution == 20 == grid.indexing.max_resolution + 1
    assert z.centroid == (0.0, 0.0)
    # Six zero vertices even on a pentagon path -- DGGAL does not drop one here.
    assert z.vertices == tuple([(0.0, 0.0)] * 6)


@pytest.mark.parametrize("grid", Z7_GRIDS, ids=lambda g: g.name)
def test_z7_level20_has_no_neighbours(grid):
    """Regression: the k-ring used to reflect the degenerate vertices back
    through quantize and raise a bare IndexError. pydggal returns [] here."""
    z = grid.zone_from_text(LVL20_HEX)
    assert z.neighbors == ()
    assert z.disk(1) == frozenset({z})


@requires_pydggal
@pytest.mark.parametrize(
    "grid, oracle_name",
    [(IGEO7, "ISEA7H_Z7"), (IVEA7H, "IVEA7H_Z7"), (RTEA7H, "RTEA7H_Z7")],
    ids=lambda x: getattr(x, "name", x),
)
def test_z7_level20_neighbours_match_pydggal(grid, oracle_name):
    g = oracle_grid(oracle_name)
    z_oracle = g.getZoneFromTextID(LVL20_HEX)
    assert [int(n) for n in g.getZoneNeighbors(z_oracle)] == []
    assert grid.zone_from_text(LVL20_HEX).neighbors == ()


@pytest.mark.parametrize("grid", Z7_GRIDS, ids=lambda g: g.name)
def test_z7_level19_is_unaffected(grid):
    """The guard must fire at 20 and not one level early -- level 19 is a real zone."""
    z = grid.zone_from_text("05" + "3" * 19)
    assert z.resolution == 19 == grid.indexing.max_resolution
    assert z.centroid != (0.0, 0.0)
    assert len(z.vertices) == 6
    # ...and a pentagon path at 19 still yields 5 vertices, so we have not
    # flattened the pentagon logic on the way past.
    pent = grid.zone_from_text("05" + "0" * 19)
    assert pent.is_pentagon and len(pent.vertices) == 5


@requires_pydggal
@pytest.mark.parametrize(
    "grid, oracle_name",
    [(IGEO7, "ISEA7H_Z7"), (IVEA7H, "IVEA7H_Z7"), (RTEA7H, "RTEA7H_Z7")],
    ids=lambda x: getattr(x, "name", x),
)
@pytest.mark.parametrize("text", [LVL20_HEX, LVL20_PENT], ids=["hex-path", "pent-path"])
def test_z7_level20_matches_pydggal(grid, oracle_name, text):
    """The whole point: match the engine, not our own idea of what null means."""
    g = oracle_grid(oracle_name)
    z_oracle = g.getZoneFromTextID(text)
    ours = grid.zone_from_text(text)

    # The packed int stays bit-identical -- we degenerate the geometry, not the id.
    assert ours.value == int(z_oracle)
    assert g.getZoneLevel(z_oracle) == ours.resolution == 20

    c = g.getZoneWGS84Centroid(z_oracle)
    assert (float(c.lat), float(c.lon)) == ours.centroid
    v = g.getZoneWGS84Vertices(z_oracle)
    assert tuple((float(p.lat), float(p.lon)) for p in v) == ours.vertices


# --------------------------------------------------------------------------- #
# I3H nullZone: reachable from ordinary lat/lon, must not print as a real id
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("grid", I3H_GRIDS, ids=lambda g: g.name)
def test_i3h_null_zone_is_reachable_and_labelled(grid):
    """A near-pole point at res 0 quantizes to the sentinel; it must say so."""
    z = grid.zone_from_geo(*NULL_POINT)
    assert z.value == NULL_ZONE
    assert z.text_id == NULL_TEXT == "(null)"


@pytest.mark.parametrize("grid", I3H_GRIDS, ids=lambda g: g.name)
def test_i3h_null_text_round_trips(grid):
    """from_text(to_text(v)) == v must hold for every value to_text can emit."""
    assert grid.indexing.to_text(NULL_ZONE) == NULL_TEXT
    assert grid.indexing.from_text(NULL_TEXT) == NULL_ZONE
    assert grid.zone_from_text(NULL_TEXT).value == NULL_ZONE


@pytest.mark.parametrize("grid", I3H_GRIDS, ids=lambda g: g.name)
def test_i3h_null_text_is_not_confusable_with_a_real_id(grid):
    """Regression: the sentinel used to stringify as `F-7FFFFFFFFFFFF-D."""
    assert grid.indexing.to_text(NULL_ZONE) != "`F-7FFFFFFFFFFFF-D"
    # A genuine zone must still get a genuine id.
    real = grid.zone_from_geo(52.0, 5.0, 5)
    assert real.text_id != NULL_TEXT


@pytest.mark.parametrize("bad", ["null", "", "(NULL)", "(null) ", "junk"])
def test_i3h_only_the_canonical_spelling_maps_to_null(bad):
    """DGGAL returns nullZone for any unparseable id; this port keeps raising.

    Documented deliberate deviation -- see I3HIndexing.from_text.
    """
    from py4dggs.types import InvalidZoneError

    with pytest.raises(InvalidZoneError):
        ISEA3H.indexing.from_text(bad)


@requires_pydggal
@pytest.mark.parametrize(
    "grid, oracle_name",
    [(ISEA3H, "ISEA3H"), (IVEA3H, "IVEA3H"), (RTEA3H, "RTEA3H")],
    ids=lambda x: getattr(x, "name", x),
)
def test_i3h_null_zone_matches_pydggal(grid, oracle_name):
    """Same sentinel, same text id, and (still) the same meaningless centroid."""
    from _pydggal_oracle import geopoint

    g = oracle_grid(oracle_name)
    lat, lon, res = NULL_POINT
    z_oracle = g.getZoneFromWGS84Centroid(res, geopoint(lat, lon))
    ours = grid.zone_from_geo(lat, lon, res)

    assert int(z_oracle) == ours.value == NULL_ZONE
    assert g.getZoneTextID(z_oracle) == ours.text_id == NULL_TEXT
    # The centroid stays garbage in BOTH engines -- deliberately not "fixed",
    # since matching DGGAL is the contract.
    c = g.getZoneWGS84Centroid(z_oracle)
    assert float(c.lat) == pytest.approx(ours.centroid.lat, abs=1e-9)
    assert float(c.lon) == pytest.approx(ours.centroid.lon, abs=1e-9)


# --------------------------------------------------------------------------- #
# Protocol conformance -- the declared contract must hold for what we ship
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "impl", [Z7Indexing(), I3HIndexing()], ids=["Z7Indexing", "I3HIndexing"]
)
def test_every_concrete_indexing_satisfies_the_protocol(impl):
    """Regression: I3HIndexing used to fail this, because the Protocol declared
    the congruent-hierarchy methods (parent/num_children/child_digits) as
    required even though the non-congruent I3H hierarchy lives on the topology."""
    assert isinstance(impl, Indexing)
    assert isinstance(impl.max_resolution, int)


@pytest.mark.parametrize(
    "impl",
    [HexAperture7Topology(), HexAperture3Topology()],
    ids=["HexAperture7Topology", "HexAperture3Topology"],
)
def test_every_concrete_topology_satisfies_the_protocol(impl):
    assert isinstance(impl, Topology)
    assert isinstance(impl.aperture, int)


def test_null_geometry_is_an_optional_capability():
    """hex_a7 declares it (level 20); hex_a3 must NOT -- the I3H sentinel keeps
    DGGAL's meaningless-but-nonzero geometry, so an override would diverge."""
    assert hasattr(HexAperture7Topology(), "is_null_geometry")
    assert not hasattr(HexAperture3Topology(), "is_null_geometry")
    # And the sentinel's centroid on an I3H grid is therefore NOT (0,0).
    assert Zone(ISEA3H, NULL_ZONE).centroid != (0.0, 0.0)
