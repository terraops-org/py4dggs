"""`py4dggs.geojson` - RFC 7946 output for cell boundaries.

`Grid.vertices` is a faithful mirror of DGGAL's `getZoneWGS84Vertices`: raw
boundary points, longitudes in [-180, 180], no dateline handling. That is
correct and deliberately NOT changed here (the conformance suite asserts
bit-equality against pydggal). This module is the presentation layer above it:

  * a cell straddling the antimeridian has vertices on both sides, so a naive
    ring spans ~360 deg of longitude and renders as a sliver wrapped the wrong
    way round the globe. It must become a two-part `MultiPolygon`.
  * a cell enclosing a pole winds +/-360 deg in longitude but carries no vertex
    at lat +/-90, so the pole has to be closed over explicitly.

These tests assert *properties* (closure, coordinate order, hemisphere
containment, winding, the centroid staying inside) rather than frozen
coordinates, so they keep their meaning if the split algorithm is refined.
"""
import pytest

from py4dggs import IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H, RTEA3H
from py4dggs.geojson import feature_collection, zone_feature, zone_geometry
from py4dggs.types import InvalidZoneError

ALL_GRIDS = [IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H, RTEA3H]

# Known exercisers, found by probing and cross-checked against pydggal (both
# engines return identical vertices for these cells).
ORDINARY = (52.0, 5.0, 5)        # Netherlands, nowhere near a seam
ANTIMERIDIAN = (-20.0, 179.99, 5)
# Polar cells are located with `pole_cell()` below rather than by a lat/lon
# probe -- see that helper for why quantizing the exact pole is unreliable.


# --- helpers (test-only; the module itself stays dependency-free and small) --- #

def rings_of(geom):
    """Every exterior ring in a Polygon/MultiPolygon geometry."""
    if geom["type"] == "Polygon":
        return [geom["coordinates"][0]]
    return [poly[0] for poly in geom["coordinates"]]


def signed_area(ring):
    """Planar shoelace. Positive == counterclockwise == RFC 7946 exterior."""
    return 0.5 * sum(
        ring[i][0] * ring[i + 1][1] - ring[i + 1][0] * ring[i][1]
        for i in range(len(ring) - 1)
    )


def point_in_ring(lon, lat, ring):
    """Ray casting on the closed ring."""
    inside = False
    for i in range(len(ring) - 1):
        x1, y1 = ring[i]
        x2, y2 = ring[i + 1]
        if (y1 > lat) != (y2 > lat):
            xin = x1 + (lat - y1) * (x2 - x1) / (y2 - y1)
            if lon < xin:
                inside = not inside
    return inside


def point_in_geometry(lon, lat, geom):
    return any(point_in_ring(lon, lat, r) for r in rings_of(geom))


def winding_of(zone):
    """Total longitude turn around a cell's raw boundary, computed here rather
    than imported, so these tests don't check the module against its own logic.
    ~+/-360 means the cell encircles a pole."""
    lons = [v.lon for v in zone.vertices]
    total = 0.0
    for i in range(len(lons)):
        step = lons[(i + 1) % len(lons)] - lons[i]
        step -= 360.0 * round(step / 360.0)
        total += step
    return total


def pole_cell(grid, res, pole=90.0):
    """The cell that actually encloses the given pole at `res`.

    `zone_from_geo(+/-90, 0, res)` cannot be trusted to return it: quantizing
    the exact pole is a boundary tie-break (the same exact-singularity class the
    conformance suites xfail), so it often yields an adjacent cell instead -
    on RTEA7H, at every resolution tried. The pole's true cell is that cell or
    one of its neighbours; exactly one of them winds +/-360.
    """
    seed = grid.zone_from_geo(pole, 0.0, res)
    for candidate in [seed] + list(seed.neighbors):
        if abs(abs(winding_of(candidate)) - 360.0) < 1.0:
            return candidate
    pytest.skip(f"{grid.name}: no cell enclosing pole {pole} near the seed at res {res}")


def assert_coordinates_in_range(geom, label):
    for ring in rings_of(geom):
        for lon, lat in ring:
            assert -180.0 <= lon <= 180.0, f"{label}: lon {lon} out of range"
            assert -90.0 <= lat <= 90.0, f"{label}: lat {lat} out of range"


# --- structure ------------------------------------------------------------- #

def test_ordinary_cell_is_a_polygon():
    z = IGEO7.zone_from_geo(*ORDINARY)
    geom = zone_geometry(z)
    assert geom["type"] == "Polygon"
    assert len(geom["coordinates"]) == 1, "a cell has no holes"


def test_coordinates_are_lon_lat_not_lat_lon():
    """The single most likely bug in this module: GeoJSON is [x, y] = [lon, lat],
    the reverse of GeoPoint(lat, lon)."""
    z = IGEO7.zone_from_geo(*ORDINARY)          # lat ~52, lon ~5 -- unambiguous
    ring = rings_of(zone_geometry(z))[0]
    for lon, lat in ring:
        assert 4.0 < lon < 6.0, f"first element must be longitude, got {lon}"
        assert 51.0 < lat < 53.0, f"second element must be latitude, got {lat}"


def test_ring_is_closed():
    for grid in ALL_GRIDS:
        z = grid.zone_from_geo(*ORDINARY)
        for ring in rings_of(zone_geometry(z)):
            assert ring[0] == ring[-1], f"{grid.name}: ring must repeat its first point"
            assert len(ring) >= 4, "a closed triangle is the minimum"


# --- antimeridian ---------------------------------------------------------- #

def test_antimeridian_cell_becomes_a_multipolygon():
    for grid in ALL_GRIDS:
        z = grid.zone_from_geo(*ANTIMERIDIAN)
        raw = z.vertices
        span = max(v.lon for v in raw) - min(v.lon for v in raw)
        assert span > 180, f"{grid.name}: probe cell is not actually a straddler"
        assert zone_geometry(z)["type"] == "MultiPolygon"


def test_split_parts_each_stay_within_one_hemisphere():
    for grid in ALL_GRIDS:
        z = grid.zone_from_geo(*ANTIMERIDIAN)
        rings = rings_of(zone_geometry(z))
        assert len(rings) == 2, f"{grid.name}: a convex cell cuts into exactly 2 parts"
        for ring in rings:
            lons = [p[0] for p in ring]
            assert max(lons) - min(lons) < 180, (
                f"{grid.name}: part still wraps the wrong way ({max(lons) - min(lons)} deg)"
            )


def test_split_parts_stay_in_valid_coordinate_range():
    """A wrong translation after the cut is exactly what pushes a part past
    -180; the generic range test never reaches a split cell."""
    for grid in ALL_GRIDS:
        z = grid.zone_from_geo(*ANTIMERIDIAN)
        assert_coordinates_in_range(zone_geometry(z), f"{grid.name} split")


def test_split_parts_are_not_degenerate_slivers():
    for grid in ALL_GRIDS:
        z = grid.zone_from_geo(*ANTIMERIDIAN)
        for ring in rings_of(zone_geometry(z)):
            distinct = {(round(x, 9), round(y, 9)) for x, y in ring}
            assert len(distinct) >= 3, f"{grid.name}: part collapsed to {len(distinct)} points"
            assert abs(signed_area(ring)) > 1e-12, f"{grid.name}: part has no area"


def test_ordinary_cell_is_not_split():
    for grid in ALL_GRIDS:
        z = grid.zone_from_geo(*ORDINARY)
        assert zone_geometry(z)["type"] == "Polygon", f"{grid.name}: split a cell that never crosses"


# --- poles ----------------------------------------------------------------- #

@pytest.mark.parametrize("pole", [90.0, -90.0], ids=["north", "south"])
def test_polar_cap_is_closed_over_the_pole(pole):
    for grid in ALL_GRIDS:
        z = pole_cell(grid, 6, pole)
        geom = zone_geometry(z)
        ring = rings_of(geom)[0]
        corners = [i for i, (_, lat) in enumerate(ring) if abs(lat - pole) < 1e-9]
        assert len(corners) == 2, (
            f"{grid.name}: a cell enclosing pole {pole} needs exactly two lat={pole} corners"
        )
        assert abs(corners[0] - corners[1]) == 1, (
            f"{grid.name}: the two pole corners must be adjacent, so the seam is one edge; "
            f"got indices {corners}"
        )
        lons = [p[0] for p in ring]
        assert min(lons) <= -180.0 + 1e-9 and max(lons) >= 180.0 - 1e-9, (
            f"{grid.name}: polar cap must span the full longitude range"
        )
        assert_coordinates_in_range(geom, f"{grid.name} pole {pole}")


@pytest.mark.parametrize("pole", [90.0, -90.0], ids=["north", "south"])
def test_polar_cap_keeps_every_original_vertex(pole):
    """Regression: an early cut-selection bug kept only the vertices past the
    cut, silently collapsing a 6-vertex cap to a triangle."""
    for grid in ALL_GRIDS:
        z = pole_cell(grid, 6, pole)
        ring = rings_of(zone_geometry(z))[0]
        for v in z.vertices:
            assert any(
                abs(lat - v.lat) < 1e-9 and abs(abs(lon) - abs(v.lon)) < 1e-9
                for lon, lat in ring
            ), f"{grid.name}: original vertex {v} vanished from the closed cap"


@pytest.mark.parametrize("pole", [90.0, -90.0], ids=["north", "south"])
def test_polar_pentagon_cap_is_closed_over_the_pole(pole):
    """The res-0 cap is a 5-vertex pentagon, not a hexagon - a separate path."""
    for grid in ALL_GRIDS:
        z = pole_cell(grid, 0, pole)
        assert len(z.vertices) == 5, f"{grid.name}: expected a pentagon at res 0"
        geom = zone_geometry(z)
        assert any(abs(lat - pole) < 1e-9 for _, lat in rings_of(geom)[0])
        assert_coordinates_in_range(geom, f"{grid.name} pentagon pole {pole}")


def test_ordinary_cell_gets_no_pole_vertex():
    for grid in ALL_GRIDS:
        z = grid.zone_from_geo(*ORDINARY)
        for ring in rings_of(zone_geometry(z)):
            assert not any(abs(abs(lat) - 90.0) < 1e-9 for _, lat in ring), (
                f"{grid.name}: invented a pole vertex on an ordinary cell"
            )


# --- RFC 7946 conformance -------------------------------------------------- #

def test_exterior_rings_are_counterclockwise():
    """RFC 7946 s3.1.6: exterior rings follow the right-hand rule."""
    for grid in ALL_GRIDS:
        cells = {
            "ordinary": grid.zone_from_geo(*ORDINARY),
            "antimeridian": grid.zone_from_geo(*ANTIMERIDIAN),
            "polar cap": pole_cell(grid, 6),
            "polar pentagon": pole_cell(grid, 0),
        }
        for label, z in cells.items():
            for ring in rings_of(zone_geometry(z)):
                assert signed_area(ring) > 0, (
                    f"{grid.name} {label}: exterior ring is clockwise"
                )


def test_all_coordinates_are_in_valid_range():
    for grid in ALL_GRIDS:
        for lat in (-89.0, -45.0, 0.0, 45.0, 89.0):
            for lon in (-179.9, -90.0, 0.0, 90.0, 179.9):
                z = grid.zone_from_geo(lat, lon, 4)
                for ring in rings_of(zone_geometry(z)):
                    for x, y in ring:
                        assert -180.0 <= x <= 180.0, f"{grid.name}: lon {x} out of range"
                        assert -90.0 <= y <= 90.0, f"{grid.name}: lat {y} out of range"


def test_centroid_stays_inside_its_own_cell_geometry():
    """The strong property check: whatever the split does, the cell's own
    centroid must remain inside the rendered shape."""
    for grid in ALL_GRIDS:
        for lat in (-70.0, -30.0, 0.0, 30.0, 70.0):
            for lon in (-179.5, -60.0, 0.0, 60.0, 179.5):
                z = grid.zone_from_geo(lat, lon, 5)
                c = z.centroid
                geom = zone_geometry(z)
                assert point_in_geometry(c.lon, c.lat, geom), (
                    f"{grid.name}: centroid {c} fell outside its own cell "
                    f"({geom['type']}, probe lat={lat} lon={lon})"
                )


# --- Feature / FeatureCollection wrappers ---------------------------------- #

def test_zone_feature_has_geojson_feature_shape():
    z = IGEO7.zone_from_geo(*ORDINARY)
    f = zone_feature(z)
    assert f["type"] == "Feature"
    assert f["geometry"] == zone_geometry(z)
    assert f["properties"]["zone"] == z.text_id


def test_zone_feature_properties_are_overridable():
    z = IGEO7.zone_from_geo(*ORDINARY)
    f = zone_feature(z, properties={"value": 42})
    assert f["properties"] == {"value": 42}, "caller's properties replace the default"


def test_feature_collection_wraps_many_zones():
    z = IGEO7.zone_from_geo(*ORDINARY)
    fc = feature_collection([z] + list(z.neighbors))
    assert fc["type"] == "FeatureCollection"
    assert len(fc["features"]) == 7
    assert all(f["type"] == "Feature" for f in fc["features"])


def test_feature_collection_accepts_any_iterable():
    z = IGEO7.zone_from_geo(*ORDINARY)
    fc = feature_collection(iter([z]))
    assert len(fc["features"]) == 1


# --- null zones ------------------------------------------------------------ #

# Z7 level 20 is representable in the 20-slot packing but has no geometry:
# DGGAL yields six (0,0) vertices. Same ids test_null_zone.py pins.
LVL20_HEX = "05" + "3" * 20
LVL20_PENT = "05" + "0" * 20


@pytest.mark.parametrize("text", [LVL20_HEX, LVL20_PENT], ids=["hex-path", "pent-path"])
@pytest.mark.parametrize("grid", [IGEO7, IVEA7H, RTEA7H], ids=lambda g: g.name)
def test_null_geometry_zone_raises(grid, text):
    """Emitting a degenerate (0,0) polygon would silently corrupt a caller's
    output, so a zone with no geometry raises instead - matching the library's
    InvalidZoneError convention."""
    z = grid.zone_from_text(text)
    assert z.vertices == tuple([(0.0, 0.0)] * 6), "precondition: this is a null zone"
    with pytest.raises(InvalidZoneError):
        zone_geometry(z)


def test_null_geometry_zone_raises_from_feature_helpers_too():
    z = IGEO7.zone_from_text(LVL20_HEX)
    with pytest.raises(InvalidZoneError):
        zone_feature(z)
    with pytest.raises(InvalidZoneError):
        feature_collection([z])


def test_feature_collection_properties_may_be_computed_per_zone():
    """A FeatureCollection whose features all carry identical properties is
    near-useless; the common case is properties derived from each zone."""
    zones = [IGEO7.zone_from_geo(*ORDINARY)]
    zones += list(zones[0].neighbors)
    fc = feature_collection(zones, properties=lambda z: {"id": z.text_id, "res": z.resolution})
    assert [f["properties"]["id"] for f in fc["features"]] == [z.text_id for z in zones]
    assert all(f["properties"]["res"] == 5 for f in fc["features"])


def test_zone_feature_explicit_empty_properties_stay_empty():
    """`properties={}` is a deliberate 'no properties', distinct from omitting
    the argument. Pins the `is None` check against a later falsy 'cleanup'."""
    z = IGEO7.zone_from_geo(*ORDINARY)
    assert zone_feature(z, properties={})["properties"] == {}
    assert zone_feature(z)["properties"] == {"zone": z.text_id}
