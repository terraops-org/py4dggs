"""pydggal (DGGAL Python binding) test oracle for the Phase-1 grids.

pydggal IS the DGGAL ground truth (there is no readable Python oracle for the
new grids). Skips cleanly if `dggal` is not installed.

Binding notes discovered 2026-07-01 (dggal==0.0.6 / ecrt==0.0.6):
  - getZoneWGS84Centroid / getZoneWGS84Vertices return GeoPoint whose .lat/.lon
    are `ecrt.Degrees`, NOT float -> coerce with float() (full ~17-digit double).
  - getZoneNeighbors/getZoneParents/getZoneChildren raise OverflowError in the
    binding for ~45% of Z7 zones (base cells 8-11, bit 63 set): the wrapper
    marshals a C uint64 through a signed int64. The ENGINE is fine, so we call
    libdggal.so directly through dggal's own cffi connector -- see _zone_array.
    Callers no longer need to guard, and no zone is skipped.
  - py4dggs geometry agrees with pydggal to ~1e-12 deg on centroids (NOT bit-
    identical: the C engine and our port differ in the last few bits), so float
    comparisons use a tight tolerance, not `==`. Discrete outputs (textId,
    neighbour text-id sets) DO match exactly.
"""
import importlib.util
import pytest

HAVE_PYDGGAL = importlib.util.find_spec("dggal") is not None
requires_pydggal = pytest.mark.skipif(not HAVE_PYDGGAL, reason="dggal (pydggal) not installed")

_APP = None
def _app():
    global _APP
    if _APP is None:
        import dggal
        from dggal import Application, pydggal_setup
        _APP = Application(appGlobals=vars(dggal))
        pydggal_setup(_APP)
    return _APP

def oracle_grid(name: str):
    """Return a pydggal DGGRS handle by class name, e.g. 'ISEA7H_Z7' / 'IVEA7H_Z7'."""
    _app()
    import dggal
    return getattr(dggal, name)()

# --- the int64 bypass (lessons #11/#22) ------------------------------------- #
# dggal 0.0.6's high-level Array copy marshals a C uint64 DGGRSZone through a
# SIGNED int64 (ecrt.py TA(): u.i64 = a), so getZoneNeighbors/Parents/Children
# raise OverflowError for every Z7 zone with bit 63 set -- base cells 8-11,
# ~45% of zones (measured 134/300 at res 5). The engine is fine: libdggal.so
# fills a correct C uint64 array, and cffi hands each element back as an
# arbitrary-precision Python int. So we skip the broken wrapper and call the C
# function through the binding's own connector. Guard-and-skip (the previous
# approach) silently dropped base cells 8-11 from every neighbour/hierarchy
# check -- exactly where prior porting bugs hid.
_ZONE_ARRAY_CAP = {  # eC's own getMaxNeighbors/getMaxParents/getMaxChildren
    "DGGRS_getZoneNeighbors": 6,
    "DGGRS_getZoneParents": 3,
    "DGGRS_getZoneChildren": 13,
}


def _zone_array(fn_name: str, grid, value: int) -> list[int]:
    """Call a DGGRS zone-array C function directly, returning the ints.

    `getZoneNeighbors` takes a trailing nbTypes out-parameter, which we do not
    need; NULL is accepted (the eC signature marks it optional).
    """
    import dggal
    lib, ffi = dggal.lib, dggal.ffi
    arr = ffi.new(f"eC_DGGRSZone[{_ZONE_ARRAY_CAP[fn_name]}]")
    fn = getattr(lib, fn_name)
    if fn_name == "DGGRS_getZoneNeighbors":
        n = fn(grid.impl, value, arr, ffi.NULL)
    else:
        n = fn(grid.impl, value, arr)
    return [int(arr[i]) for i in range(n)]


def geopoint(lat: float, lon: float):
    from dggal import GeoPoint
    g = GeoPoint(); g.lat = lat; g.lon = lon
    return g

def forward_textid(grid, lat: float, lon: float, res: int) -> str:
    return grid.getZoneTextID(grid.getZoneFromWGS84Centroid(res, geopoint(lat, lon)))

def centroid_of(grid, textid: str) -> tuple[float, float]:
    c = grid.getZoneWGS84Centroid(grid.getZoneFromTextID(textid))
    return float(c.lat), float(c.lon)

def vertices_of(grid, textid: str) -> list[tuple[float, float]]:
    z = grid.getZoneFromTextID(textid)
    return [(float(p.lat), float(p.lon)) for p in grid.getZoneWGS84Vertices(z)]

def neighbours_of(grid, textid: str) -> set[str]:
    """Neighbour text-id set. Never None: goes around the binding's int64 limit
    via :func:`_zone_array`, so base cells 8-11 are covered like any other."""
    z = grid.getZoneFromTextID(textid)
    return {grid.getZoneTextID(n) for n in _zone_array("DGGRS_getZoneNeighbors", grid, z)}

def oracle_zone_int(grid, lat: float, lon: float, res: int) -> int:
    return int(grid.getZoneFromWGS84Centroid(res, geopoint(lat, lon)))

def centroid_of_int(grid, value: int) -> tuple[float, float]:
    c = grid.getZoneWGS84Centroid(value)
    return float(c.lat), float(c.lon)

def vertices_of_int(grid, value: int) -> list[tuple[float, float]]:
    return [(float(p.lat), float(p.lon)) for p in grid.getZoneWGS84Vertices(value)]

def neighbour_centroids_of_int(grid, value: int, nd: int = 5) -> set:
    """Set of neighbour centroids (rounded). Never None -- see :func:`_zone_array`."""
    ns = _zone_array("DGGRS_getZoneNeighbors", grid, value)
    out = set()
    for n in ns:
        c = grid.getZoneWGS84Centroid(n)
        out.add((round(float(c.lat), nd), round(float(c.lon), nd)))
    return out


def neighbour_ints_of_int(grid, value: int) -> set:
    """Neighbour zones as raw ints (== our packed I3H value). Never None -- I3H
    never overflowed, and Z7 no longer does either (see :func:`_zone_array`)."""
    return set(_zone_array("DGGRS_getZoneNeighbors", grid, value))


# --- hierarchy oracle (I3H does not overflow) -------------------------------- #
def parents_of_int(grid, value: int) -> set:
    return set(_zone_array("DGGRS_getZoneParents", grid, value))


def children_of_int(grid, value: int) -> set:
    return set(_zone_array("DGGRS_getZoneChildren", grid, value))


def centroid_parent_of_int(grid, value: int):
    v = int(grid.getZoneCentroidParent(value))
    return None if v == (1 << 64) - 1 else v


def centroid_child_of_int(grid, value: int) -> int:
    return int(grid.getZoneCentroidChild(value))


def is_centroid_child_of_int(grid, value: int) -> bool:
    return bool(grid.isZoneCentroidChild(value))


# --- sub-zones oracle (I3H does not overflow) -------------------------------- #
def subzone_ints_of_int(grid, value: int, depth: int) -> list[int]:
    """Ordered sub-zone ints at `depth` below `value` (`getSubZones`)."""
    zones = grid.getSubZones(value, depth)
    return [int(zones[i]) for i in range(zones.count)]


def count_subzones_of_int(grid, value: int, depth: int) -> int:
    return grid.countSubZones(value, depth)
