from py4dggs.types import GridConfig
from py4dggs.projections.isea import ISEAProjection
from py4dggs.topologies.hex_a7 import HexAperture7Topology

PROJ = ISEAProjection(); GEOM = PROJ.build_geometry(GridConfig())
TOPO = HexAperture7Topology()


def test_aperture_is_7():
    assert TOPO.aperture == 7


def test_quantize_anchor_lisbon():
    """Self-contained regression anchor (no oracle needed): Lisbon res 5 must
    quantize to base 0, digits [6,4,1,5,6] -- the canonical "0064156" cell
    exercised throughout this repo's docs/tests."""
    p = PROJ.forward(GEOM, 38.7223, -9.1393)
    base, digits = TOPO.quantize(GEOM, p, 5)
    assert (base, digits) == (0, [6, 4, 1, 5, 6])


def test_planar_vertices_count():
    hexv = TOPO.planar_vertices(GEOM, 0, [6, 4, 1, 5, 6])    # hexagon
    penta = TOPO.planar_vertices(GEOM, 0, [0, 0, 0, 0, 0])   # centre path -> pentagon
    assert len(hexv) == 6 and len(penta) == 5
