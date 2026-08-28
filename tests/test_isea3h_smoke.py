from _pydggal_oracle import (requires_pydggal, oracle_grid, oracle_zone_int,
                             centroid_of_int, vertices_of_int, neighbour_centroids_of_int)

@requires_pydggal
def test_pydggal_isea3h_int_lever():
    g = oracle_grid("ISEA3H")
    v = oracle_zone_int(g, 38.7223, -9.1393, 5)          # Lisbon res 5
    assert v == 306244774661193870                        # C2-23-C, verified 2026-07-02
    clat, clon = centroid_of_int(g, v)
    assert abs(clat - 37.73156) < 1e-4 and abs(clon - (-9.22933)) < 1e-4
    assert len(vertices_of_int(g, v)) == 6
    assert len(neighbour_centroids_of_int(g, v)) == 6


def test_isea3h_registered_and_constructs():
    from py4dggs import ISEA3H, get_grid, Zone
    assert get_grid("ISEA3H") is ISEA3H
    z = ISEA3H.zone_from_geo(38.7223, -9.1393, 5)
    assert z.resolution == 5 and len(z.vertices) == 6
    assert repr(z).startswith("Zone('")   # canonical rhombic text-id, no crash
