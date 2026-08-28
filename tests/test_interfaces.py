# tests/test_interfaces.py
from py4dggs.types import GeoPoint, PlanarPoint, GridConfig, InvalidZoneError
from py4dggs.interfaces import Projection, Topology, Indexing

def test_geopoint_is_namedtuple():
    p = GeoPoint(38.7, -9.1)
    assert p.lat == 38.7 and p.lon == -9.1

def test_planarpoint_carries_face():
    q = PlanarPoint(face=3, x=2.5, y=2.25)
    assert q.face == 3 and q.x == 2.5 and q.y == 2.25

def test_gridconfig_defaults_to_igeo7():
    c = GridConfig()
    assert c.orientation_lon_deg == -11.20 and c.authalic is True

def test_protocols_define_methods():
    # structural: the protocols must declare the methods later tasks rely on
    assert hasattr(Projection, "forward") and hasattr(Projection, "inverse")
    assert hasattr(Topology, "quantize") and hasattr(Topology, "planar_centroid")
    assert hasattr(Indexing, "encode") and hasattr(Indexing, "decode")
