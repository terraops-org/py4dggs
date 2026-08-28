"""Task 0 gate: pydggal is wired and exposes IVEA7H_Z7."""
from _pydggal_oracle import requires_pydggal, oracle_grid, forward_textid, centroid_of

@requires_pydggal
def test_ivea7h_z7_is_exposed_and_works():
    g = oracle_grid("IVEA7H_Z7")                      # must not AttributeError
    tid = forward_textid(g, 38.7223, -9.1393, 5)
    assert isinstance(tid, str) and len(tid) >= 2     # a real IVEA7H text id
    clat, clon = centroid_of(g, tid)                  # Degrees -> float works
    assert -90.0 <= clat <= 90.0 and -180.0 <= clon <= 180.0

@requires_pydggal
def test_isea7h_z7_oracle_matches_known_anchor():
    g = oracle_grid("ISEA7H_Z7")
    assert forward_textid(g, 38.7223, -9.1393, 5) == "0064156"   # pydggal ISEA7H anchor
