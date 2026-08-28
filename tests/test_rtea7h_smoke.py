"""Task smoke: RTEA7H is registered, constructs cells, has geometry distinct
from IGEO7 and IVEA7H, and leaves IGEO7/IVEA7H unaffected. Full correctness vs
pydggal (forward / centroid / vertices / neighbours) is verified in
test_rtea7h_fuzz.py."""
from py4dggs import IGEO7, IVEA7H, RTEA7H, get_grid, Zone


def test_rtea7h_registered_and_constructs():
    assert get_grid("RTEA7H") is RTEA7H
    z = RTEA7H.zone_from_geo(38.7223, -9.1393, 5)
    assert isinstance(z, Zone) and z.resolution == 5 and len(z.text_id) >= 2
    assert z.centroid is not None            # geometry runs without error


def test_rtea7h_geometry_differs_from_igeo7_and_ivea7h():
    # Same cell id, different projection -> a different geographic centroid.
    zi = IGEO7.zone_from_text("0064156")
    zv = IVEA7H.zone_from_text("0064156")
    zr = RTEA7H.zone_from_text("0064156")
    assert (zi.centroid.lat, zi.centroid.lon) != (zr.centroid.lat, zr.centroid.lon)
    assert (zv.centroid.lat, zv.centroid.lon) != (zr.centroid.lat, zr.centroid.lon)


def test_igeo7_and_ivea7h_unaffected_by_new_grid():
    assert IGEO7.zone_from_geo(38.7223, -9.1393, 5).text_id == "0064156"
    zv = IVEA7H.zone_from_geo(38.7223, -9.1393, 5)
    assert isinstance(zv, Zone) and zv.resolution == 5
