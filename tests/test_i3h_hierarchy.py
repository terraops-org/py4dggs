"""A2 I3H hierarchy (parents/children/centroid) vs pydggal (exact int-sets)."""
import random
import pytest
from py4dggs import ISEA3H, IVEA3H, RTEA3H, Zone
from py4dggs.indexings.i3h import unpack_i3h
from py4dggs.topologies.hex_a3 import _i3h_centroid_child, _i3h_is_centroid_child
from _pydggal_oracle import (requires_pydggal, oracle_grid, oracle_zone_int,
                             centroid_child_of_int, is_centroid_child_of_int)

GRIDS = [("ISEA3H", ISEA3H), ("IVEA3H", IVEA3H), ("RTEA3H", RTEA3H)]


@requires_pydggal
@pytest.mark.parametrize("name,grid", GRIDS)
def test_centroid_child_and_flag_vs_pydggal(name, grid):
    g = oracle_grid(name); rng = random.Random(11)
    for _ in range(500):
        val = oracle_zone_int(g, rng.uniform(-89, 89), rng.uniform(-180, 180), rng.randint(1, 20))
        assert _i3h_centroid_child(*unpack_i3h(val)) == centroid_child_of_int(g, val)
        assert _i3h_is_centroid_child(*unpack_i3h(val)) == is_centroid_child_of_int(g, val)


@requires_pydggal
@pytest.mark.parametrize("name,grid", GRIDS)
def test_parents_children_centroidparent_vs_pydggal(name, grid):
    from py4dggs.topologies.hex_a3 import _i3h_get_parents, _i3h_get_children, _i3h_centroid_parent
    from _pydggal_oracle import parents_of_int, children_of_int, centroid_parent_of_int
    g = oracle_grid(name); rng = random.Random(22)
    pts = [(rng.uniform(-89, 89), rng.uniform(-180, 180), rng.randint(1, 28)) for _ in range(800)]
    pts += [(-58.4, -168.8, 12), (89.99, 0, 5), (-89.99, 0, 8), (0, 11.2, 7), (58.28, 11.2, 6)]
    for lat, lon, res in pts:
        val = oracle_zone_int(g, lat, lon, res)
        assert set(_i3h_get_parents(val)) == parents_of_int(g, val)
        assert set(_i3h_get_children(val)) == children_of_int(g, val)
        assert _i3h_centroid_parent(val) == centroid_parent_of_int(g, val)


@requires_pydggal
@pytest.mark.parametrize("name,grid", GRIDS)
def test_zone_hierarchy_vs_pydggal(name, grid):
    from _pydggal_oracle import parents_of_int, children_of_int, centroid_parent_of_int
    g = oracle_grid(name); rng = random.Random(33)
    for _ in range(400):
        val = oracle_zone_int(g, rng.uniform(-89, 89), rng.uniform(-180, 180), rng.randint(1, 22))
        z = Zone(grid, val)
        assert {p.value for p in z.parents} == parents_of_int(g, val)
        assert {c.value for c in z.children} == children_of_int(g, val)
        ocp = centroid_parent_of_int(g, val)
        assert (z.centroid_parent.value if z.centroid_parent else None) == ocp
        assert z.is_centroid_child == is_centroid_child_of_int(g, val)
        # invariants
        assert z.parent == (z.parents[0] if z.parents else None)
        for c in z.children:
            assert any(p == z for p in c.parents)


def test_z7_children_and_parent_unchanged():
    from py4dggs import IGEO7
    z = IGEO7.zone_from_geo(48, 11.2, 4)
    assert len(z.children) == 7 and all(c.text_id.startswith(z.text_id) for c in z.children)  # digit-path
    assert z.parent is not None and z.parent.text_id == z.text_id[:-1]
    assert z.is_centroid_child is True and len(z.parents) == 1


@requires_pydggal
@pytest.mark.parametrize("name,grid", GRIDS)
def test_ancestor_sibling_helpers(name, grid):
    # built on the (pydggal-verified) parents/children; checked for self-consistency
    # + transitivity rather than DGGAL's finicky isZoneAncestorOf(maxDepth) signature.
    g = oracle_grid(name); rng = random.Random(44)
    for _ in range(200):
        val = oracle_zone_int(g, rng.uniform(-89, 89), rng.uniform(-180, 180), rng.randint(2, 16))
        z = Zone(grid, val)
        for p in z.parents:
            assert z.is_immediate_child_of(p) and p.is_ancestor_of(z)
            for gp in p.parents:                 # transitivity: grandparent is an ancestor
                assert gp.is_ancestor_of(z)
        for c in z.children:
            assert c.is_immediate_child_of(z) and z.is_ancestor_of(c)
        # two children that both count z among their parents are siblings
        shared = [c for c in z.children if any(p == z for p in c.parents)]
        for i in range(len(shared)):
            for j in range(i + 1, len(shared)):
                assert shared[i].is_sibling_of(shared[j])


@requires_pydggal
@pytest.mark.parametrize("name,grid", GRIDS)
def test_hierarchy_exact_any_seed(name, grid):
    """Different seed than the runs above (res 1-30) — exactness is not seed luck."""
    from _pydggal_oracle import parents_of_int, children_of_int, centroid_parent_of_int
    g = oracle_grid(name); rng = random.Random(20260703)
    for _ in range(1500):
        val = oracle_zone_int(g, rng.uniform(-89.5, 89.5), rng.uniform(-180, 180), rng.randint(1, 30))
        z = Zone(grid, val)
        assert {p.value for p in z.parents} == parents_of_int(g, val)
        assert {c.value for c in z.children} == children_of_int(g, val)
        assert (z.centroid_parent.value if z.centroid_parent else None) == centroid_parent_of_int(g, val)


def test_hierarchy_invariants_no_oracle():
    """Oracle-independent structural invariants (runs without dggal)."""
    from py4dggs.topologies.hex_a3 import _i3h_get_parents_raw
    for grid in (ISEA3H, IVEA3H, RTEA3H):
        rng = random.Random(7)
        for _ in range(1500):
            z = grid.zone_from_geo(rng.uniform(-89.5, 89.5), rng.uniform(-180, 180), rng.randint(1, 28))
            npoints = 5 if z.is_pentagon else 6
            ch = z.children
            assert len(ch) == npoints + 1                        # centroid child + one per vertex
            for c in ch:
                assert any(p == z for p in c.parents)            # z is a parent of each child
            assert len(_i3h_get_parents_raw(z.value)) in (1, 3)  # pre-null-filter parent count
            cp = z.centroid_parent
            assert cp is None or any(p == cp for p in z.parents)  # centroid_parent is among parents
