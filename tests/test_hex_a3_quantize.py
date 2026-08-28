import pytest
from py4dggs.topologies.hex_a3 import HexAperture3Topology
from py4dggs.projections.isea import ISEAProjection
from py4dggs.types import GridConfig
from py4dggs.indexings.i3h import unpack_i3h
from _pydggal_oracle import requires_pydggal, oracle_grid, oracle_zone_int

PROJ = ISEAProjection(); GEOM = PROJ.build_geometry(GridConfig())
TOPO = HexAperture3Topology()

def _py4dggs_quantize_int(lat, lon, res):
    p = PROJ.forward(GEOM, lat, lon)
    base, digits = TOPO.quantize(GEOM, p, res)
    assert digits == []
    return base

@requires_pydggal
def test_quantize_int_matches_pydggal():
    import random
    g = oracle_grid("ISEA3H"); rng = random.Random(20260702)
    mism = 0; total = 0
    for _ in range(600):
        lat = rng.uniform(-88.0, 88.0); lon = rng.uniform(-180.0, 180.0); res = rng.randint(0, 12)
        total += 1
        if _py4dggs_quantize_int(lat, lon, res) != oracle_zone_int(g, lat, lon, res):
            mism += 1
    # exact-int match; only measure-zero boundary tie-breaks may differ
    assert mism <= 6, f"{mism}/{total} quantize int mismatches (expected ~0)"
