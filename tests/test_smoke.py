def test_py4dggs_imports():
    import py4dggs            # the new package
    assert py4dggs is not None


def test_igeo7_anchor():
    """Self-contained sanity anchor -- no oracle dependency."""
    from py4dggs import IGEO7
    assert IGEO7.zone_from_geo(38.7223, -9.1393, 5).text_id == "0064156"
