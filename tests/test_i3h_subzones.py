"""A3 sub-zones: reference-vector replay, per-slice differential fuzz vs
pydggal, cross-invariants. Oracle = the high-level pydggal API (I3H's packing
cannot reach 2^63 -- see documentation/superpowers/specs/2026-07-06-a3-isea3h-subzones-design.md)."""
import pytest

from py4dggs.registry import get_grid
from py4dggs.types import InvalidZoneError
from py4dggs.topologies.hex_a3 import (
    _i3h_vertices,
    _i3h_sub_zones,
    _i3h_count_sub_zones,
    move5x6_vertex_v3,
    _crosses5x6_interruption,
    _intersects5x6_interruption,
    _i3h_is_edge_hex,
)
from py4dggs.indexings.i3h import unpack_i3h, pack_i3h, I3HIndexing

from _pydggal_oracle import (
    requires_pydggal,
    oracle_grid,
    geopoint,
    subzone_ints_of_int,
    count_subzones_of_int,
)


@requires_pydggal
def test_vertex_order_matches_pydggal_mod_wrap():
    """_i3h_vertices' index order must match eC getVertices, up to the
    documented +/-5 toroidal-wrap representation ambiguity."""
    g = oracle_grid("ISEA3H")
    cases = [
        (30.0, -70.0, 4),   # interior hexagon
        (10.0, 36.0, 2),    # near a pentagon-dense region
        (90.0, 0.0, 2),     # north pole
        (-90.0, 0.0, 2),    # south pole
    ]
    for lat, lon, level in cases:
        z = g.getZoneFromWGS84Centroid(level, geopoint(lat, lon))
        value = int(z)
        l9r, root, rix, sh = unpack_i3h(value)
        ours = _i3h_vertices(l9r, root, rix, sh)
        theirs_arr = g.getZoneCRSVertices(z, 0)
        n = theirs_arr.count
        for i in range(n):
            ox, oy = ours[i]
            tx, ty = theirs_arr[i].x, theirs_arr[i].y
            # allow the +/-5 wrap representation of the same 5x6 point
            close = (abs(ox - tx) < 1e-6 and abs(oy - ty) < 1e-6) or \
                    (abs(ox + 5 - tx) < 1e-6 and abs(oy + 5 - ty) < 1e-6) or \
                    (abs(ox - 5 - tx) < 1e-6 and abs(oy - 5 - ty) < 1e-6)
            assert close, f"vertex {i} mismatch at ({lat},{lon},{level}): ours={ (ox,oy) } theirs={ (tx,ty) }"


def test_count_sub_zones_closed_form():
    # hexagon: 7/13/37 sub-zones at depth 1/2/3 (spot-verified during A2/A3 design)
    assert _i3h_count_sub_zones(6, 1) == 7
    assert _i3h_count_sub_zones(6, 2) == 13
    assert _i3h_count_sub_zones(6, 3) == 37
    assert _i3h_count_sub_zones(6, 0) == 1


@requires_pydggal
def test_count_sub_zones_matches_pydggal():
    g = oracle_grid("ISEA3H")
    z = g.getZoneFromWGS84Centroid(4, geopoint(30.0, -70.0))
    for depth in (1, 2, 3, 4):
        assert _i3h_count_sub_zones(6, depth) == count_subzones_of_int(g, z, depth)


def test_move5x6_vertex_v3_no_crossing():
    """Test move5x6_vertex_v3 with a step that stays within one rhombus
    (no interruption crossing)."""
    # Start at (2.5, 2.5) and step by small (0.1, 0.1)
    # This stays well within a single rhombus, no crossing
    cx, cy = 2.5, 2.5
    dx, dy = 0.1, 0.1
    result = move5x6_vertex_v3(cx, cy, dx, dy)
    # Should equal the naive (cx + dx, cy + dy) when no crossing
    expected_x = cx + dx
    expected_y = cy + dy
    assert len(result) == 2
    assert abs(result[0] - expected_x) < 1e-10
    assert abs(result[1] - expected_y) < 1e-10


def test_move5x6_vertex_v3_crosses_north_interruption():
    """Test move5x6_vertex_v3 crossing a north interruption edge.

    This test crosses the segment [(0,0), (1,0)] from above.
    Start at (0.4, 0.2) and step by (0.3, -0.5) crosses this edge.
    The crossing point is near (0.6, 0), which is on the first north
    interruption edge. When crossed, the result should reflect via
    _cross5x6_interruption rather than naively land at (0.7, -0.3).

    Hand-traced derivation (verified against actual code output):
    - _crosses5x6_interruption(0.4, 0.2, 0.3, -0.5) returns:
      * src (crossing point) = (0.52, 0.0)
      * dst (reflected crossing) = [5.0, 4.48]
      * north = True (north interruption)
    - Using north branch formula: vx = i2x - 2*(dy - (i1y - cy))
      * dy - (i1y - cy) = -0.5 - (0.0 - 0.2) = -0.3
      * vx = 5.0 - 2*(-0.3) = 5.0 + 0.6 = 5.6
    - Using north branch formula: vy = i2y + dx - (i1x - cx)
      * dx - (i1x - cx) = 0.3 - (0.52 - 0.4) = 0.18
      * vy = 4.48 + 0.18 = 4.66
    - Wrap conditions: vx=5.6 not in (>5, <0) and vy=4.66 not in (>5, <1)
    - Final result: [5.6, 4.66]
    """
    cx, cy = 0.4, 0.2
    dx, dy = 0.3, -0.5

    # First verify that _crosses5x6_interruption detects this crossing
    crossed, src, dst, north = _crosses5x6_interruption(cx, cy, dx, dy)
    assert crossed, "Expected this step to cross an interruption"
    assert north is not None, "Expected north/south hemisphere info"
    assert src is not None and dst is not None, "Expected crossing and reflection points"

    # Now test move5x6_vertex_v3 returns the reflected point, not the naive point
    result = move5x6_vertex_v3(cx, cy, dx, dy)
    naive_x = cx + dx
    naive_y = cy + dy

    # Result should be the reflected point, not the naive endpoint
    # (for this case, the naive endpoint would wrap or reflect via interruption)
    assert len(result) == 2
    # Tight assertion on exact expected value (derived and verified above)
    assert abs(result[0] - 5.6) < 1e-12, f"Expected vx≈5.6, got {result[0]}"
    assert abs(result[1] - 4.66) < 1e-12, f"Expected vy≈4.66, got {result[1]}"


def test_intersects5x6_interruption_basic():
    """Test _intersects5x6_interruption with a segment that crosses an edge."""
    # Segment from (0.3, 0.2) to (0.8, -0.1) crosses horizontal edge [(0,0), (1,0)]
    a0x, a0y = 0.3, 0.2
    a1x, a1y = 0.8, -0.1
    # Horizontal edge from (0,0) to (1,0)
    b0x, b0y = 0.0, 0.0
    b1x, b1y = 1.0, 0.0

    found, ix, iy, t = _intersects5x6_interruption(a0x, a0y, a1x, a1y, b0x, b0y, b1x, b1y)

    assert found, "Expected intersection to be found"
    assert 0 <= t <= 1, f"Expected parametric t in [0,1], got {t}"
    assert ix is not None and iy is not None, "Expected intersection coordinates"
    # Intersection should be on the horizontal line y=0
    assert abs(iy - 0.0) < 1e-10, f"Expected y≈0 at crossing, got {iy}"
    # Intersection should be between x=0 and x=1
    assert -1e-10 <= ix <= 1.0 + 1e-10, f"Expected x in [0,1], got {ix}"


def test_intersects5x6_interruption_no_crossing():
    """Test _intersects5x6_interruption with parallel segments (no crossing)."""
    # Horizontal segment from (0.3, 0.5) to (0.8, 0.5)
    a0x, a0y = 0.3, 0.5
    a1x, a1y = 0.8, 0.5
    # Another horizontal edge from (0,0) to (1,0)
    b0x, b0y = 0.0, 0.0
    b1x, b1y = 1.0, 0.0

    found, ix, iy, t = _intersects5x6_interruption(a0x, a0y, a1x, a1y, b0x, b0y, b1x, b1y)

    assert not found, "Expected no intersection for parallel segments"
    assert ix is None and iy is None and t is None


def _is_edge_hex(l9r, root, rix, sh):
    """Test-side wrapper around the shared ``hex_a3._i3h_is_edge_hex`` core
    boundary predicate (final-review DRY extraction, Finding 2 -- this used to
    be an independent third copy of the production ``edge_hex`` logic). The
    explicit ``root >= 10`` / ``(rix == 0 and sh <= 1)`` pentagon/pole
    exclusion below is REDUNDANT with ``_i3h_is_edge_hex``'s own ``rix != 0``
    check (every pentagon-like cell reachable in this codebase -- non-polar
    pentagons and both poles alike -- has ``rix == 0``, confirmed during this
    review), but is kept here as explicit, readable documentation of WHY
    pentagons/poles never count as edge-hex, for this test file's own
    ``_is_interior_hex``/enumeration helpers that read it directly."""
    if root >= 10 or (rix == 0 and sh <= 1):
        return False
    return _i3h_is_edge_hex(l9r, root, rix, sh)


def _is_interior_hex(l9r, root, rix, sh):
    if root >= 10:
        return False
    if rix == 0 and sh <= 1:
        return False
    return not _is_edge_hex(l9r, root, rix, sh)


@requires_pydggal
@pytest.mark.parametrize("lat,lon,level", [
    (30.0, -70.0, 4), (10.0, -100.0, 3), (45.0, 10.0, 3), (-20.0, 60.0, 4),
])
def test_first_sub_zone_centroid_interior_hex(lat, lon, level):
    from py4dggs.topologies.hex_a3 import _i3h_first_sub_zone_centroid, _i3h_from_centroid
    g = oracle_grid("ISEA3H")
    z = g.getZoneFromWGS84Centroid(level, geopoint(lat, lon))
    value = int(z)
    l9r, root, rix, sh = unpack_i3h(value)
    for depth in (1, 2, 3, 4):
        cx, cy = _i3h_first_sub_zone_centroid(l9r, root, rix, sh, depth)
        ours = pack_i3h(*_i3h_from_centroid(level + depth, cx, cy))
        theirs = int(g.getFirstSubZone(z, depth))
        assert ours == theirs


@requires_pydggal
@pytest.mark.parametrize("grid_name", ["ISEA3H", "IVEA3H", "RTEA3H"])
def test_sub_zones_interior_hex_exact(grid_name):
    import random
    random.seed(20260706)
    g = oracle_grid(grid_name)
    matched = total = 0
    for _ in range(60):
        lat = random.uniform(-85, 85)
        lon = random.uniform(-179, 179)
        level = random.choice([2, 3, 4, 5])
        z = g.getZoneFromWGS84Centroid(level, geopoint(lat, lon))
        value = int(z)
        l9r, root, rix, sh = unpack_i3h(value)
        if not _is_interior_hex(l9r, root, rix, sh):
            continue
        depth = random.choice([1, 2, 3, 4])
        total += 1
        theirs = subzone_ints_of_int(g, z, depth)
        ours = _i3h_sub_zones(value, depth)
        if ours == theirs:
            matched += 1
    assert total >= 30, "fixture drew too few interior-hex samples -- widen the sweep"
    assert matched == total, f"{matched}/{total} interior-hex sub-zone lists matched pydggal exactly"


@requires_pydggal
@pytest.mark.parametrize("grid_name", ["ISEA3H", "IVEA3H", "RTEA3H"])
def test_sub_zones_edge_hex_exact(grid_name):
    import random
    random.seed(20260706)
    g = oracle_grid(grid_name)
    matched = total = 0
    attempts = 0
    while total < 20 and attempts < 2000:
        attempts += 1
        lat = random.uniform(-85, 85)
        lon = random.uniform(-179, 179)
        level = random.choice([2, 3, 4, 5])
        z = g.getZoneFromWGS84Centroid(level, geopoint(lat, lon))
        value = int(z)
        l9r, root, rix, sh = unpack_i3h(value)
        if not _is_edge_hex(l9r, root, rix, sh):
            continue
        depth = random.choice([1, 2, 3, 4])
        total += 1
        theirs = subzone_ints_of_int(g, z, depth)
        ours = _i3h_sub_zones(value, depth)
        if ours == theirs:
            matched += 1
    assert total >= 20, "could not find enough edge-hex samples -- widen the sweep"
    assert matched == total, f"{matched}/{total} edge-hex sub-zone lists matched pydggal exactly"


@requires_pydggal
@pytest.mark.parametrize("textid,depths", [
    ("B6-1-A", [1, 2, 3, 4, 5, 6]),  # North edge hex, even parent (A) -- doc case #7/#9
    ("B5-3-A", [1, 2, 3, 4, 5, 6]),  # South edge hex, even parent (A) -- doc case #8/#10
    ("B6-1-D", [1, 2, 3, 4, 5]),     # Same rix as B6-1-A but sub_hex=D (interior, not edge)
    ("B5-3-D", [1, 2, 3, 4, 5]),     # Same rix as B5-3-A but sub_hex=D (interior, not edge)
])
def test_sub_zones_documented_edge_vectors(textid, depths):
    """Replays parent zones cited by name in I3HSubZones.ec's own 28-case
    documentation comment (lines 11-98) as interruption-spanning ("Edge")
    examples for the odd-depth generators. ``B6-1-D``/``B5-3-D`` sit at the
    SAME rhombus position as their ``-A`` siblings but carry sub_hex D (odd,
    off-center) -- included here specifically because they are NOT edge-hex
    (see :func:`_is_edge_hex`'s docstring) and must keep matching via the
    plain interior path, guarding against a regression that widens edge-hex
    detection to ignore the sub_hex filter."""
    g = oracle_grid("ISEA3H")
    z = g.getZoneFromTextID(textid)
    value = int(z)
    for depth in depths:
        theirs = subzone_ints_of_int(g, z, depth)
        ours = _i3h_sub_zones(value, depth)
        assert ours == theirs, f"{textid} depth={depth}: mismatch"


def _zone_id_of(value):
    """Round-trip a packed I3H int through our own ``to_text`` -- avoids the
    ``dggal.I3HZone(int)`` re-wrap ``OverflowError`` noted in Task 1 Step 9."""
    return I3HIndexing().to_text(value)


@requires_pydggal
@pytest.mark.parametrize("grid_name", ["ISEA3H", "IVEA3H", "RTEA3H"])
def test_sub_zones_nonpolar_pentagon_exact(grid_name):
    """Enumerates all 10 non-polar pentagons directly (``root in range(10),
    rix == 0, sub_hex in (0, 1)``, both ``level_i9r`` parities) rather than
    random-sampling them, since they are sparse and specific -- exactly the
    brief's Step 3 test."""
    g = oracle_grid(grid_name)
    matched = total = 0
    for level_i9r in (1, 2):
        for root in range(10):
            for sub_hex in (0, 1):
                value = pack_i3h(level_i9r, root, 0, sub_hex)
                for depth in (1, 2, 3):
                    total += 1
                    z = g.getZoneFromTextID(_zone_id_of(value))
                    theirs = subzone_ints_of_int(g, z, depth)
                    ours = _i3h_sub_zones(value, depth)
                    if ours == theirs:
                        matched += 1
    assert matched == total, f"{matched}/{total} non-polar-pentagon sub-zone lists matched"


@requires_pydggal
@pytest.mark.parametrize("textid,depths", [
    ("A6-0-B", [1, 2, 3, 4]),  # North non-polar pentagon, odd parent -- doc cases #13/#19 region
    ("A5-0-B", [1, 2, 3, 4]),  # South non-polar pentagon, odd parent -- doc cases #14/#20 region
    ("A2-0-A", [1, 2, 3, 4, 5, 6]),  # North non-polar pentagon, even parent -- doc cases #15/#17
    ("A3-0-A", [1, 2, 3, 4, 5, 6]),  # South non-polar pentagon, even parent -- doc cases #16/#18
])
def test_sub_zones_documented_pentagon_vectors(textid, depths):
    """Replays named non-polar-pentagon zones (``nv == 5``, confirmed via
    ``getZoneWGS84Vertices`` count) as a permanent regression fixture,
    independent of the enumeration sweep above.

    NOTE: I3HSubZones.ec's own 28-case documentation comment (lines 11-98)
    cites parent zone IDs like ``A6-0-D``/``A5-0-D``/``A6-0-E`` for these
    "Non-polar pentagons" cases (#13-20) -- but those strings use a stale or
    different zone-ID letter scheme (``A6-0-E`` would decode ``subHex=4``,
    out of our valid A-D/0-3 range) and do NOT round-trip through this
    codebase's ``to_text``/``getZoneFromTextID`` to an actual ``nv==5``
    pentagon: querying pydggal directly shows ``A6-0-D``/``A5-0-D`` are
    ordinary ``nv==6`` hexagons (``sub_hex=D=3``), not the pentagons the
    comment's section header claims. Verified independently and replaced
    with the genuine non-polar-pentagon zones ``A6-0-B``/``A5-0-B`` (same
    root/rix, ``sub_hex=B=1``, confirmed ``nv==5``) instead of trusting the
    doc string literally."""
    g = oracle_grid("ISEA3H")
    z = g.getZoneFromTextID(textid)
    value = int(z)
    for depth in depths:
        theirs = subzone_ints_of_int(g, z, depth)
        ours = _i3h_sub_zones(value, depth)
        assert ours == theirs, f"{textid} depth={depth}: mismatch"


@requires_pydggal
@pytest.mark.parametrize("grid_name", ["ISEA3H", "IVEA3H", "RTEA3H"])
def test_sub_zones_polar_pentagon_exact(grid_name):
    """Task 4 Step 3: enumerates BOTH polar pentagons directly -- the North
    (``root == 10``) and South (``root == 11``) icosahedron vertices, ``rix ==
    0`` always, both ``sub_hex`` values (poles are pentagons, ``sub_hex <= 1``),
    several ``level_i9r`` and depths -- rather than random-sampling (there are
    only 2 such cells per level, and forward-quantizing a lat/lon almost never
    lands exactly on a pole). Depths run to 6 so the multi-stage
    ``cross5x6Interruption`` chaining and ``oddRow`` parity corrections that only
    fire at larger ``r`` are actually exercised (a shallow green is false
    confidence for this cell class)."""
    g = oracle_grid(grid_name)
    matched = total = 0
    for level_i9r in (0, 1, 2, 3):
        for root in (10, 11):
            for sub_hex in (0, 1):
                value = pack_i3h(level_i9r, root, 0, sub_hex)
                z = g.getZoneFromTextID(_zone_id_of(value))
                assert int(z) == value, f"pole {_zone_id_of(value)} did not round-trip"
                for depth in (1, 2, 3, 4, 5, 6):
                    total += 1
                    theirs = subzone_ints_of_int(g, z, depth)
                    ours = _i3h_sub_zones(value, depth)
                    if ours == theirs:
                        matched += 1
    assert matched == total, f"{matched}/{total} polar-pentagon sub-zone lists matched"


# I3HSubZones.ec:72-94 documents FIRST sub-zone text-IDs for the 8 polar-pentagon
# cases (#21-28), spanning all 4 (oddParent, oddDepth) parity dispatches x N/S.
# Each tuple: (label, root, level_i9r, sub_hex, depth, comment_FIRST_string).
# Parents are constructed directly from (pole root, level_i9r, sub_hex) derived
# from the comment's P-LEVEL/SZ-LEVEL/DEPTH columns -- NOT by decoding the
# comment's own parent IDs (``A0-0-G``/``A9-0-H``/``A0-0-B``/``A9-0-C``), whose
# ``G``/``H``/``B``/``C`` sub-hex letters use a STALE lettering scheme neither
# our ``to_text`` nor pydggal's native ``getZoneTextID`` reproduces (same
# staleness Task 3 hit for the non-polar pentagons).
_EC_POLAR_REFERENCE_VECTORS = [
    ("#21 N oddParent/oddDepth", 10, 0, 1, 3, "C2-6-A"),
    ("#21 N oddParent/oddDepth", 10, 0, 1, 5, "D2-12-A"),
    ("#22 S oddParent/oddDepth", 11, 0, 1, 3, "C9-36-A"),
    ("#22 S oddParent/oddDepth", 11, 0, 1, 5, "D9-1E6-A"),
    ("#23 N evenParent/oddDepth", 10, 0, 0, 3, "B0-5-D"),
    ("#23 N evenParent/oddDepth", 10, 0, 0, 5, "C0-21-D"),
    ("#24 S evenParent/oddDepth", 11, 0, 0, 3, "B9-7-D"),
    ("#24 S evenParent/oddDepth", 11, 0, 0, 5, "C9-39-D"),
    ("#25 N evenParent/evenDepth", 10, 0, 0, 4, "C0-21-A"),
    ("#26 S evenParent/evenDepth", 11, 0, 0, 4, "C9-39-A"),
    ("#27 N oddParent/evenDepth", 10, 0, 1, 4, "C0-6-D"),
    ("#27 N oddParent/evenDepth", 10, 0, 1, 6, "D0-12-D"),
    ("#28 S oddParent/evenDepth", 11, 0, 1, 4, "C1-36-D"),
    ("#28 S oddParent/evenDepth", 11, 0, 1, 6, "D1-1E6-D"),
]


@requires_pydggal
@pytest.mark.parametrize("label,root,level_i9r,sub_hex,depth,comment_first", _EC_POLAR_REFERENCE_VECTORS)
def test_reference_vectors_from_ec_comment_table_polar(label, root, level_i9r, sub_hex, depth, comment_first):
    """Task 4 Step 7: replay the eC comment table's documented FIRST sub-zone for
    each polar-pentagon parity case (I3HSubZones.ec cases #21-28), belt-and-
    suspenders against BOTH a live pydggal call AND the literal comment string.

    Two independent failure modes are checked:
      1. **Port bug** -- our first sub-zone (indeed the whole list) must equal
         pydggal's, asserted as packed ints (the exact-match signal).
      2. **Comment transcription** -- the comment's FIRST text-ID must match,
         WITH a documented exception: for the 6 odd-level-sub-zone cases the
         comment writes the sub-hex letter ``D`` where the current DGGAL scheme
         (both our ``to_text`` AND pydggal's own ``getZoneTextID``) writes ``B``
         for the identical zone. That is a stale-scheme artifact in the eC source
         comment itself (verified: pydggal's binary returns ``B0-5-B`` while its
         own doc comment says ``B0-5-D``), so we assert the ``root-rix`` prefix
         literally and normalize only the trailing sub-hex letter via the
         verified D->B mapping. The even-level (``-A``) cases match verbatim."""
    g = oracle_grid("ISEA3H")
    ix = I3HIndexing()
    value = pack_i3h(level_i9r, root, 0, sub_hex)
    z = g.getZoneFromTextID(ix.to_text(value))

    theirs = subzone_ints_of_int(g, z, depth)
    ours = _i3h_sub_zones(value, depth)

    assert ours == theirs, f"{label} depth={depth}: sub-zone list mismatch vs pydggal"

    our_first_tid = ix.to_text(ours[0])
    py_first_tid = g.getZoneTextID(ours[0])
    assert our_first_tid == py_first_tid, f"{label}: our to_text disagrees with pydggal native"

    exp_prefix, exp_sub = comment_first.rsplit("-", 1)
    got_prefix, got_sub = our_first_tid.rsplit("-", 1)
    assert got_prefix == exp_prefix, (
        f"{label} depth={depth}: first sub-zone prefix {got_prefix!r} != comment {exp_prefix!r}"
    )
    normalized = {"D": "B"}.get(exp_sub, exp_sub)
    assert got_sub == normalized, (
        f"{label} depth={depth}: sub-hex letter {got_sub!r} != comment {exp_sub!r} "
        f"(normalized {normalized!r}); comment uses a stale sub-hex scheme"
    )


# --------------------------------------------------------------------------- #
# Task 5: Grid/Zone dispatch + generic index/at-index + cross-invariants
# --------------------------------------------------------------------------- #
def test_grid_sub_zones_dispatch():
    grid = get_grid("ISEA3H")
    z = grid.zone_from_geo(30.0, -70.0, 4)
    subs = grid.sub_zones(z.value, 2)
    assert isinstance(subs, tuple)
    assert len(subs) == grid.count_sub_zones(z.value, 2)
    assert subs[0] == grid.first_sub_zone(z.value, 2)


def test_subzone_cross_invariants_no_oracle():
    grid = get_grid("ISEA3H")
    z = grid.zone_from_geo(30.0, -70.0, 3)
    for depth in (1, 2, 3):
        subs = grid.sub_zones(z.value, depth)
        assert grid.count_sub_zones(z.value, depth) == len(subs)
        assert subs[0] == grid.first_sub_zone(z.value, depth)
        for i in (0, len(subs) // 2, len(subs) - 1):
            assert grid.sub_zone_at_index(z.value, depth, i) == subs[i]
            assert grid.sub_zone_index(z.value, subs[i]) == i


def test_zone_sub_zones_wrapper():
    grid = get_grid("ISEA3H")
    z = grid.zone_from_geo(30.0, -70.0, 3)
    subs = z.sub_zones(2)
    assert len(subs) == z.count_sub_zones(2)
    assert subs[0] == z.first_sub_zone(2)
    assert z.sub_zone_index(subs[0]) == 0


# --------------------------------------------------------------------------- #
# Final-review polish: relative_depth<=0 edge cases + Z7 NotImplementedError
# --------------------------------------------------------------------------- #
def test_sub_zone_index_shallower_returns_minus_one_but_the_zone_itself_is_zero():
    """`Grid.sub_zone_index` for candidates that are not strictly deeper.

    An ANCESTOR is not a sub-zone -> -1. The zone ITSELF, however, is its own
    sub-zone at index 0: `sub_zones(v, 0) == (v,)`, `count_sub_zones(v, 0) == 1`
    and `sub_zone_at_index(v, 0, 0) == v` all already said so.

    Corrected 2026-08-21: this test previously asserted -1 for the zone itself,
    folding depth 0 into the "shallower" rejection. That was never checked
    against the oracle -- pydggal's own `getSubZoneIndex(v, v)` returns **0**
    (and -1 for the parent), so the old assertion diverged from DGGAL."""
    grid = get_grid("ISEA3H")
    z = grid.zone_from_geo(30.0, -70.0, 4)
    v = z.value
    assert grid.sub_zone_index(v, v) == 0  # relative_depth == 0: the zone itself
    parents = grid.parents(v)
    assert grid.sub_zone_index(v, parents[0]) == -1  # relative_depth == -1: shallower


def test_subzone_relative_depth_zero_is_the_zone_itself():
    """relative_depth == 0 means "the zone itself" -- verified against pydggal's
    getSubZones(zone, 0), which returns exactly [zone] for an interior hexagon,
    a non-polar pentagon AND a polar pentagon alike (all three checked directly
    against the oracle during this review; see the fix in Grid.count_sub_zones/
    first_sub_zone/sub_zones). This is a `Grid`-level generic short-circuit
    (mirroring dggrs.ec's generic getSubZones), NOT delegated to the I3H
    geometric generators, which are not designed for depth 0 -- calling into
    them directly for relative_depth == 0 raises/misbehaves (`3 ** ((0 - 2) //
    2)` produces a float exponent); confirmed by reading
    `_i3h_sub_zone_centroids`/`_gen_even_parent_even_depth` before writing this
    test, not assumed."""
    grid = get_grid("ISEA3H")
    z = grid.zone_from_geo(30.0, -70.0, 4)
    v = z.value
    assert grid.count_sub_zones(v, 0) == 1
    assert grid.first_sub_zone(v, 0) == v
    assert grid.sub_zones(v, 0) == (v,)


def test_z7_grid_sub_zones_raises_not_implemented():
    """The Z7 (congruent-digit) grids have no `Topology.sub_zones` override --
    only the I3H family's `HexAperture3Topology` implements it (confirmed by
    reading `src/py4dggs/topologies/hex_a7.py`: no count_sub_zones/first_sub_zone/
    sub_zones methods on `HexAperture7Topology`). `Grid.sub_zones` must raise
    `NotImplementedError` for a Z7 grid such as IGEO7, both for a normal depth
    and for relative_depth == 0 (the depth-0 short-circuit lives in `Grid`, but
    strictly AFTER the `getattr(...) is None` check, so it must never mask the
    missing-capability error for a grid that has no sub-zone order at all)."""
    grid = get_grid("IGEO7")
    z = grid.zone_from_geo(48, 11.2, 4)
    with pytest.raises(NotImplementedError):
        grid.sub_zones(z.value, 3)
    with pytest.raises(NotImplementedError):
        grid.sub_zones(z.value, 0)


def test_negative_relative_depth_rejected():
    """count_sub_zones/first_sub_zone/sub_zones/sub_zone_at_index previously
    only special-cased relative_depth == 0, so a negative depth fell through
    to the geometric generators: count_sub_zones silently returned 1 (as if
    depth were 0), first_sub_zone silently returned an unrelated (shallower)
    zone, sub_zones raised an unhandled TypeError from `3 ** ((depth-1)//2)`
    producing a float exponent, and sub_zone_at_index silently returned the
    same wrong zone (trusting count_sub_zones' bogus answer for its bounds
    check). All four must now reject relative_depth < 0 the same way the
    sibling `sub_zone_index` already did (found via /code-review, fixed in
    `Grid._sub_zone_fn`)."""
    grid = get_grid("ISEA3H")
    z = grid.zone_from_geo(30.0, -70.0, 4)
    v = z.value
    for relative_depth in (-1, -2):
        with pytest.raises(ValueError):
            grid.count_sub_zones(v, relative_depth)
        with pytest.raises(ValueError):
            grid.first_sub_zone(v, relative_depth)
        with pytest.raises(ValueError):
            grid.sub_zones(v, relative_depth)
        with pytest.raises(ValueError):
            grid.sub_zone_at_index(v, relative_depth, 0)


def test_sub_zone_past_max_resolution_rejected():
    """count_sub_zones/first_sub_zone/sub_zones previously never checked the
    resulting resolution against `Indexing.max_resolution` (unlike
    `zone_from_geo`, which enforces it strictly) -- a zone already at the max
    resolution (33 for I3H) could silently produce sub-zones at resolution 35
    with no error. Now rejected the same way `zone_from_geo` rejects an
    out-of-range resolution directly (found via /code-review, fixed in
    `Grid._sub_zone_fn`)."""
    grid = get_grid("ISEA3H")
    assert grid.indexing.max_resolution == 33
    z = grid.zone_from_geo(30.0, -70.0, 33)
    v = z.value
    with pytest.raises(InvalidZoneError):
        grid.count_sub_zones(v, 2)
    with pytest.raises(InvalidZoneError):
        grid.first_sub_zone(v, 2)
    with pytest.raises(InvalidZoneError):
        grid.sub_zones(v, 2)
    # exactly at the boundary (resolution 33) must still work
    assert grid.count_sub_zones(v, 0) == 1
