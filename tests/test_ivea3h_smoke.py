"""Task: IVEA3H smoke — the IVEA projection dropped into the same aperture-3
I3H topology/indexing as ISEA3H (a pure combination, zero new implementation
code). Mirrors test_isea3h_smoke.py; full correctness vs pydggal (quantize /
centroid / vertices / neighbours) is verified in test_ivea3h_fuzz.py."""
from _pydggal_oracle import (requires_pydggal, oracle_grid, oracle_zone_int,
                             centroid_of_int, vertices_of_int, neighbour_centroids_of_int)

@requires_pydggal
def test_pydggal_ivea3h_int_lever():
    g = oracle_grid("IVEA3H")
    v = oracle_zone_int(g, 38.7223, -9.1393, 5)          # Lisbon res 5
    assert v == 306244774661193870                        # same packed I3H cell as ISEA3H at this point, verified 2026-07-03
    clat, clon = centroid_of_int(g, v)
    assert abs(clat - 37.68431) < 1e-4 and abs(clon - (-9.28996)) < 1e-4
    assert len(vertices_of_int(g, v)) == 6
    assert len(neighbour_centroids_of_int(g, v)) == 6


def test_ivea3h_registered_and_constructs():
    from py4dggs import IVEA3H, get_grid, Zone
    assert get_grid("IVEA3H") is IVEA3H
    z = IVEA3H.zone_from_geo(38.7223, -9.1393, 5)
    assert z.resolution == 5 and len(z.vertices) == 6
    assert repr(z).startswith("Zone('")   # canonical rhombic text-id, no crash


def test_ivea3h_geometry_differs_from_isea3h():
    # Same packed I3H cell id (identical uint64 value: same rhombus root/index/
    # sub-hex), different projection -> a different geographic centroid.
    from py4dggs import ISEA3H, IVEA3H, Zone
    val = 306244774661193870
    zi = Zone(ISEA3H, val)
    zv = Zone(IVEA3H, val)
    assert (zi.centroid.lat, zi.centroid.lon) != (zv.centroid.lat, zv.centroid.lon)


def test_isea3h_unaffected_by_new_grid():
    from py4dggs import ISEA3H
    z = ISEA3H.zone_from_geo(38.7223, -9.1393, 5)
    assert z.value == 306244774661193870
