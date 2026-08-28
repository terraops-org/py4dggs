import pytest
from py4dggs.topologies.hex_a3 import HexAperture3Topology
from py4dggs.projections.isea import ISEAProjection
from py4dggs.types import GridConfig
from py4dggs.indexings.i3h import pack_i3h, unpack_i3h
from _pydggal_oracle import requires_pydggal, oracle_grid, oracle_zone_int, vertices_of_int

PROJ = ISEAProjection(); GEOM = PROJ.build_geometry(GridConfig())
TOPO = HexAperture3Topology()

def _srt(pairs, nd=4):
    return sorted((round(a, nd), round(b, nd)) for a, b in pairs)

def _py4dggs_vertices(value):
    ps = TOPO.planar_vertices(GEOM, value, [])
    return [ (PROJ.inverse(GEOM, p).lat, PROJ.inverse(GEOM, p).lon) for p in ps ]

@requires_pydggal
def test_vertices_match_pydggal():
    g = oracle_grid("ISEA3H")
    pts = [("Lisbon",38.7223,-9.1393),("Nairobi",-1.2921,36.8219),("Reykjavik",64.1466,-21.9426),
           ("NorthPole",90.0,0.0),("SouthPole",-90.0,0.0)]
    fails = []
    for name, lat, lon in pts:
        for res in range(0, 8):
            v = oracle_zone_int(g, lat, lon, res)
            dv = _py4dggs_vertices(v)
            ov = vertices_of_int(g, v)
            # exact vertex count (pentagon vs hexagon) + set match within ~11 m
            if len(dv) != len(ov) or _srt(dv) != _srt(ov):
                fails.append(f"{name} r{res} ({unpack_i3h(v)}): n={len(dv)}/{len(ov)}")
    assert not fails, "vertex mismatches: " + "; ".join(fails[:15])

@requires_pydggal
def test_vertices_match_pydggal_explicit_pole_cells():
    # The test above uses (90,0)/(-90,0), which quantize to root-0 cells, NOT
    # the pole pentagons (root 10/11). Build those explicitly to exercise the
    # 5-vertex pole fans in getVertices (RI3H.ec:1693-1759).
    g = oracle_grid("ISEA3H")
    fails = []
    for li9r in (1, 2, 3):
        for sh in (0, 1):
            for root, label in ((10, "north"), (11, "south")):
                v = pack_i3h(li9r, root, 0, sh)
                dv = _py4dggs_vertices(v)
                ov = vertices_of_int(g, v)
                if len(dv) != 5 or len(ov) != 5 or _srt(dv) != _srt(ov):
                    fails.append(
                        f"{label} apex li9r={li9r} sh={sh} ({unpack_i3h(v)}): "
                        f"n={len(dv)}/{len(ov)} py4dggs={_srt(dv)} oracle={_srt(ov)}"
                    )
    assert not fails, "pole-cell vertex mismatches:\n" + "\n".join(fails[:15])
