# src/py4dggs/geojson.py
"""RFC 7946 GeoJSON output for cell boundaries.

`Grid.vertices` mirrors DGGAL's ``getZoneWGS84Vertices`` exactly: the raw
boundary points, longitudes in [-180, 180], in the engine's own order. That is
the library's correctness contract (the conformance suite asserts bit-equality
against pydggal), so **nothing here changes it**. This module sits above it and
solves the two problems that raw contract leaves to the caller:

* **Antimeridian.** A cell straddling +/-180 has vertices on both sides, so the
  naive ring spans ~360 deg of longitude. Any planar consumer -- GeoJSON
  readers, Shapely, Leaflet, PostGIS -- draws a sliver wrapped the wrong way
  round the globe. Such a cell is cut at the antimeridian into a two-part
  ``MultiPolygon``, per RFC 7946 s3.1.9.

* **Poles.** A cell enclosing a pole winds +/-360 deg in longitude but carries
  no vertex at lat +/-90, so it cannot be closed as a planar ring at all. The
  boundary is re-anchored to run from lon -180 to +180 and then closed over the
  pole itself.

Every emitted exterior ring is closed (first point repeated) and wound
counterclockwise, per RFC 7946 s3.1.6.

The cut interpolates latitude linearly in longitude at the +/-180 meridian.
That is what mainstream antimeridian splitters do and is accurate to well under
a metre at any DGGS cell size; it is not a great-circle intersection.
"""
from __future__ import annotations

import math
from typing import Iterable

from py4dggs.types import InvalidZoneError

__all__ = ["zone_geometry", "zone_feature", "feature_collection"]

_FULL_TURN = 360.0
_HALF_TURN = 180.0
# A ring whose total longitude turn is within this of +/-360 deg encloses a pole.
_POLE_WINDING_TOLERANCE = 1.0


# --- longitude bookkeeping ------------------------------------------------- #

def _unwrap(ring: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Make longitudes continuous, so consecutive steps are never > 180 deg.

    The result may leave [-180, 180]; that is the point. It turns a ring that
    "jumps" across the antimeridian into a single connected path that can be
    reasoned about with ordinary planar geometry.
    """
    out = [ring[0]]
    for lon, lat in ring[1:]:
        prev = out[-1][0]
        delta = lon - prev
        delta -= _FULL_TURN * round(delta / _FULL_TURN)
        out.append((prev + delta, lat))
    return out


def _winding(unwrapped: list[tuple[float, float]]) -> float:
    """Total longitude turn around the closed ring, including the closing edge.

    ~0 for an ordinary cell; ~+/-360 when the ring encircles a pole.
    """
    first, last = unwrapped[0][0], unwrapped[-1][0]
    closing = first - last
    closing -= _FULL_TURN * round(closing / _FULL_TURN)
    return (last - first) + closing


def _oriented(ring: list[tuple[float, float]]) -> tuple[list[tuple[float, float]], float]:
    """Unwrap, and flip the traversal so the winding is non-negative."""
    unwrapped = _unwrap(ring)
    turn = _winding(unwrapped)
    if turn < 0:
        unwrapped = _unwrap(list(reversed(ring)))
        turn = -turn
    return unwrapped, turn


def _wrap_lon(lon: float) -> float:
    """Fold a longitude into [-180, 180], leaving the poles' +/-180 intact."""
    if -_HALF_TURN <= lon <= _HALF_TURN:
        return lon
    return lon - _FULL_TURN * round(lon / _FULL_TURN)


# --- ring hygiene ---------------------------------------------------------- #

def _signed_area(ring: list[tuple[float, float]]) -> float:
    """Planar shoelace over a closed ring. Positive == counterclockwise."""
    return 0.5 * sum(
        ring[i][0] * ring[i + 1][1] - ring[i + 1][0] * ring[i][1]
        for i in range(len(ring) - 1)
    )


def _close_and_orient(ring: list[tuple[float, float]]) -> list[list[float]]:
    """Repeat the first point and force counterclockwise winding (RFC 7946
    s3.1.6). Applied to every emitted ring, so no earlier step has to reason
    about traversal direction."""
    closed = list(ring)
    if closed[0] != closed[-1]:
        closed.append(closed[0])
    if _signed_area(closed) < 0:
        closed.reverse()
    return [[lon, lat] for lon, lat in closed]


def _is_degenerate(ring: list[tuple[float, float]]) -> bool:
    """A clipped part that collapsed to a line or a point carries no area."""
    if len(ring) < 3:
        return True
    distinct = {(round(lon, 9), round(lat, 9)) for lon, lat in ring}
    if len(distinct) < 3:
        return True
    closed = list(ring) + [ring[0]]
    return abs(_signed_area(closed)) < 1e-12


# --- the antimeridian cut -------------------------------------------------- #

def _clip_half_plane(
    ring: list[tuple[float, float]], meridian: float, keep_below: bool
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip of a ring against a vertical line in lon.

    `keep_below` keeps the lon <= meridian side, else the lon >= side. Latitude
    at the cut is interpolated linearly in longitude.
    """
    out: list[tuple[float, float]] = []
    n = len(ring)
    for i in range(n):
        lon1, lat1 = ring[i]
        lon2, lat2 = ring[(i + 1) % n]
        in1 = lon1 <= meridian if keep_below else lon1 >= meridian
        in2 = lon2 <= meridian if keep_below else lon2 >= meridian
        if in1:
            out.append((lon1, lat1))
        if in1 != in2 and lon2 != lon1:
            t = (meridian - lon1) / (lon2 - lon1)
            out.append((meridian, lat1 + t * (lat2 - lat1)))
    return out


def _split_at_antimeridian(
    unwrapped: list[tuple[float, float]], meridian: float
) -> list[list[tuple[float, float]]]:
    """Cut at +180 or -180 and bring the far side back into [-180, 180]."""
    below = _clip_half_plane(unwrapped, meridian, keep_below=True)
    above = _clip_half_plane(unwrapped, meridian, keep_below=False)
    # Exactly one side lies outside [-180, 180] and must be translated back.
    shift_below, shift_above = (0.0, -_FULL_TURN) if meridian > 0 else (_FULL_TURN, 0.0)
    parts = [
        [(lon + shift_below, lat) for lon, lat in below],
        [(lon + shift_above, lat) for lon, lat in above],
    ]
    return [p for p in parts if p and not _is_degenerate(p)]


# --- the polar cap --------------------------------------------------------- #

def _close_over_pole(
    unwrapped: list[tuple[float, float]], turn: float
) -> list[tuple[float, float]]:
    """Re-anchor a pole-encircling boundary to run -180 -> +180, then close it
    across the pole.

    `unwrapped` is already oriented so longitude increases along the ring, and
    `turn` is ~+360. The ring is re-cut at whichever antimeridian crossing it
    contains, so the emitted boundary starts at lon -180 and ends at +180; the
    two pole corners then seal it.
    """
    start_lon = unwrapped[0][0]
    # The path spans [start_lon, start_lon + 360). Find the antimeridian inside it.
    k = math.ceil((start_lon + _HALF_TURN) / _FULL_TURN)
    cut = -_HALF_TURN + _FULL_TURN * k

    # Walk the closed turn: the ring plus its own first point, one turn on.
    path = unwrapped + [(start_lon + _FULL_TURN, unwrapped[0][1])]
    lat_at_cut = None
    for i in range(len(path) - 1):
        lon1, lat1 = path[i]
        lon2, lat2 = path[i + 1]
        if lon1 <= cut <= lon2 and lon2 != lon1:
            lat_at_cut = lat1 + (cut - lon1) / (lon2 - lon1) * (lat2 - lat1)
            break
    if lat_at_cut is None:                       # pragma: no cover - defensive
        lat_at_cut = path[0][1]

    # Rotate the cyclic ring to begin at the first vertex past the cut, lifting
    # the wrapped-around tail by one turn so longitudes ascend across [cut,
    # cut + 360). Rotating (rather than filtering by longitude) is what keeps
    # every vertex: the cut rarely sits at the ring's own starting longitude.
    count = len(unwrapped)
    start = next((i for i, (lon, _) in enumerate(unwrapped) if lon >= cut), count)
    ordered = []
    for step in range(count):
        index = start + step
        lon, lat = unwrapped[index % count]
        ordered.append((lon + (_FULL_TURN if index >= count else 0.0), lat))

    shift = -_FULL_TURN * k                      # a whole number of turns: a true re-wrap
    boundary = [(cut + shift, lat_at_cut)]
    boundary += [(lon + shift, lat) for lon, lat in ordered]
    boundary.append((cut + _FULL_TURN + shift, lat_at_cut))

    pole_lat = 90.0 if sum(lat for _, lat in unwrapped) > 0 else -90.0
    boundary.append((_HALF_TURN, pole_lat))
    boundary.append((-_HALF_TURN, pole_lat))
    return boundary


# --- public API ------------------------------------------------------------ #

def _rings_for(zone) -> list[list[list[float]]]:
    vertices = zone.vertices
    if not vertices or all(v.lat == 0.0 and v.lon == 0.0 for v in vertices):
        raise InvalidZoneError(
            f"zone {zone.value} has no geometry (DGGAL nullZone); nothing to render"
        )

    ring = [(v.lon, v.lat) for v in vertices]
    unwrapped, turn = _oriented(ring)

    if abs(turn - _FULL_TURN) < _POLE_WINDING_TOLERANCE:
        return [_close_and_orient(_close_over_pole(unwrapped, turn))]

    lons = [lon for lon, _ in unwrapped]
    for meridian in (_HALF_TURN, -_HALF_TURN):
        if min(lons) < meridian < max(lons):
            parts = _split_at_antimeridian(unwrapped, meridian)
            if len(parts) > 1:
                return [_close_and_orient(p) for p in parts]
            unwrapped = parts[0] if parts else unwrapped
            break

    return [_close_and_orient([(_wrap_lon(lon), lat) for lon, lat in unwrapped])]


def zone_geometry(zone) -> dict:
    """The zone's boundary as an RFC 7946 geometry.

    Returns a ``Polygon`` normally, or a two-part ``MultiPolygon`` when the cell
    straddles the antimeridian. A cell enclosing a pole comes back as a
    ``Polygon`` closed over the pole. Raises `InvalidZoneError` for a zone with
    no geometry (DGGAL's nullZone), rather than emitting a degenerate shape.
    """
    rings = _rings_for(zone)
    if len(rings) == 1:
        return {"type": "Polygon", "coordinates": [rings[0]]}
    return {"type": "MultiPolygon", "coordinates": [[r] for r in rings]}


def zone_feature(zone, properties: dict | None = None) -> dict:
    """The zone as a GeoJSON Feature. Defaults to carrying its text id; pass
    `properties` to replace that wholesale."""
    return {
        "type": "Feature",
        "geometry": zone_geometry(zone),
        "properties": {"zone": zone.text_id} if properties is None else properties,
    }


def feature_collection(zones: Iterable, properties=None) -> dict:
    """A FeatureCollection over any iterable of zones.

    `properties` may be a dict applied to every feature, or a callable taking a
    zone and returning its properties.
    """
    resolve = properties if callable(properties) else lambda _zone: properties
    return {
        "type": "FeatureCollection",
        "features": [zone_feature(zone, resolve(zone)) for zone in zones],
    }
