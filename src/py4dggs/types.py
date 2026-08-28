# src/py4dggs/types.py
"""Shared value types for the DGGS reference."""
from __future__ import annotations
from dataclasses import dataclass
from typing import NamedTuple


class GeoPoint(NamedTuple):
    """A geographic coordinate in decimal degrees (WGS84 lat/lon)."""
    lat: float
    lon: float


class PlanarPoint(NamedTuple):
    """A point in a projection's planar space. `face` is the icosahedron face
    id for icosahedral projections (or a projection-specific sentinel, e.g. -1
    for "derive from x,y"); x,y are that projection's planar coords (the 5x6
    oblique grid for ISEA/IVEA/RTEA)."""
    face: int
    x: float
    y: float


@dataclass(frozen=True)
class GridConfig:
    """Orientation + datum handling for an icosahedral grid (mirrors DGGAL).
    Defaults are the canonical IGEO7 values."""
    orientation_lat_deg: float = 31.7174744114611
    orientation_lon_deg: float = -11.20
    azimuth_deg: float = 0.0
    authalic: bool = True


class InvalidZoneError(ValueError):
    """Raised for malformed or non-canonical zone ids / out-of-range resolution."""


NULL_ZONE = (1 << 64) - 1
"""DGGAL's ``nullZone`` sentinel (dggrs.ec:9) — "no such zone".

Shared by the indexings and topologies so the sentinel has ONE definition. DGGAL
returns it from a quantization whose guards reject the point; it is a real,
reachable outcome of ``zone_from_geo`` for the I3H grids very near the poles at
the coarsest resolutions (both DGGAL and this port produce it, and both then
produce meaningless geometry for it — see ``NULL_TEXT``)."""

NULL_TEXT = "(null)"
"""What DGGAL's ``getZoneTextID`` prints for :data:`NULL_ZONE`.

This is DGGAL's *only* signal that a zone is the null sentinel, so the text-id
serialisation mirrors it exactly rather than stringifying the sentinel's bit
fields into a normal-looking id."""
