"""Task 2 smoke: IVEA7H is registered, constructs cells, has geometry distinct
from IGEO7, and leaves IGEO7 unaffected. Full correctness vs pydggal (forward /
centroid / vertices / neighbours) is verified in test_ivea7h_fuzz.py (Task 3)."""
from py4dggs import IGEO7, IVEA7H, get_grid, Zone


def test_ivea7h_registered_and_constructs():
    assert get_grid("IVEA7H") is IVEA7H
    z = IVEA7H.zone_from_geo(38.7223, -9.1393, 5)
    assert isinstance(z, Zone) and z.resolution == 5 and len(z.text_id) >= 2
    assert z.centroid is not None            # geometry runs without error


def test_ivea7h_geometry_differs_from_igeo7():
    # Same cell id, different projection -> a different geographic centroid.
    zi = IGEO7.zone_from_text("0064156")
    zv = IVEA7H.zone_from_text("0064156")
    assert (zi.centroid.lat, zi.centroid.lon) != (zv.centroid.lat, zv.centroid.lon)


def test_igeo7_unaffected_by_new_grid():
    assert IGEO7.zone_from_geo(38.7223, -9.1393, 5).text_id == "0064156"
