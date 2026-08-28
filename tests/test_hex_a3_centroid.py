import pytest
from py4dggs.topologies.hex_a3 import HexAperture3Topology
from py4dggs.projections.isea import ISEAProjection
from py4dggs.types import GridConfig
from py4dggs.indexings.i3h import pack_i3h, unpack_i3h
from _pydggal_oracle import requires_pydggal, oracle_grid, oracle_zone_int, centroid_of_int

PROJ = ISEAProjection(); GEOM = PROJ.build_geometry(GridConfig())
TOPO = HexAperture3Topology()

def _py4dggs_centroid(value):
    li9r, root, rix, sh = unpack_i3h(value)
    p = TOPO.planar_centroid(GEOM, value, [])
    c = PROJ.inverse(GEOM, p)
    return c.lat, c.lon

@requires_pydggal
def test_centroid_matches_pydggal_cities_and_poles():
    g = oracle_grid("ISEA3H")
    pts = [("Lisbon",38.7223,-9.1393),("Tokyo",35.6762,139.6503),("Sydney",-33.8688,151.2093),
           ("NorthPole",90.0,0.0),("SouthPole",-90.0,0.0),("Antimeridian",0.0,179.9)]
    fails = []
    for name, lat, lon in pts:
        for res in range(0, 8):
            v = oracle_zone_int(g, lat, lon, res)
            olat, olon = centroid_of_int(g, v)
            dlat, dlon = _py4dggs_centroid(v)
            if abs(dlat - olat) > 1e-7 or abs(dlon - olon) > 1e-7:
                fails.append(f"{name} r{res} ({unpack_i3h(v)}): py4dggs=({dlat:.6f},{dlon:.6f}) vs ({olat:.6f},{olon:.6f})")
    assert not fails, "centroid mismatches:\n" + "\n".join(fails[:15])

@requires_pydggal
def test_centroid_matches_pydggal_explicit_pole_cells():
    # The brief's pole test above uses (90,0)/(-90,0), which quantize to root-0
    # cells and never exercise the actual pole PENTAGONS (root 10/11). Build
    # those explicitly: for a spread of aperture-3 levels and both sub-hex
    # selectors that occur at the pole (0 = A/even, 1 = B/odd-centroid), check
    # the py4dggs centroid against the pydggal oracle. These apex-pentagon
    # centroids sit at ~(58.4 deg, ...) -- the icosahedron apexes, NOT the
    # geographic poles -- that is expected (root 10/11 are "North"/"South" only
    # by DGGAL's own comment; the true poles are ordinary root-0..9 cells).
    g = oracle_grid("ISEA3H")
    fails = []
    for li9r in (1, 2, 3):
        for sh in (0, 1):
            for root, label in ((10, "north"), (11, "south")):
                v = pack_i3h(li9r, root, 0, sh)
                olat, olon = centroid_of_int(g, v)
                dlat, dlon = _py4dggs_centroid(v)
                if abs(dlat - olat) > 1e-7 or abs(dlon - olon) > 1e-7:
                    fails.append(
                        f"{label} apex li9r={li9r} sh={sh} ({unpack_i3h(v)}): "
                        f"py4dggs=({dlat:.6f},{dlon:.6f}) vs ({olat:.6f},{olon:.6f})"
                    )
    assert not fails, "pole-cell centroid mismatches:\n" + "\n".join(fails[:15])
