"""Aperture-3 hexagonal topology on the 5x6 ISEA planar grid.

This is the topology half of ISEA3H/IVEA3H (the I3H family). Unlike the Z7
aperture-7 topology (:mod:`py4dggs.topologies.hex_a7`), I3H is *not* addressed by a
base cell + variable-length digit path: :func:`py4dggs.indexings.i3h.pack_i3h` packs
the whole cell (level, rhombus root, linear rhombus index, sub-hex selector) into
one integer, which the ``Topology`` protocol's ``base`` parameter carries whole
(``digits`` stays ``[]`` — see ``I3HIndexing`` in ``py4dggs/indexings/i3h.py``).

This module implements the full ``Topology`` protocol for I3H: ``quantize``
(``I3HZone::fromCentroid``), ``planar_centroid`` (``I3HZone::centroid``) and
``planar_vertices`` (``I3HZone::getVertices``), plus the 5x6 coordinate helpers
they share. Hierarchy (parent/children) stays out of scope, mirroring how
``I3HIndexing`` defers its own out-of-scope hierarchy methods.

Provenance (the *why* of each block is cited inline against these sources):
  - ``dggrs/RI3H.ec:2138-2168``   — the ``centroid`` property: the non-polar
                                    rhombus anchor (``rowOP``/``colOP``/``ixOP``
                                    decomposition of ``rhombusIX``) plus the
                                    even/odd sub-hex (A/B/C/D) offset, and the two
                                    polar shortcuts (root 10/11).
  - ``projections/ri5x6.ec:1388-1434`` — ``move5x6Vertex`` (v1): offset a vertex
                                    from a centroid, re-expressing the offset
                                    across a rhombus interruption ("dent") if it
                                    is crossed. NOT ``move5x6Vertex2`` (the
                                    aperture-7 topology's ``move5x6_vertex`` is a
                                    port of that sibling, v2 — a different,
                                    "cross-early" variant); this is the plain v1.
  - ``projections/ri5x6.ec:1562-1605`` — the exported ``canonicalize5x6``: fold a
                                    5x6 point into the fundamental domain,
                                    detecting the polar diagonals and rhombus
                                    interruption seams via the ``(int)floor``-cast
                                    row/col, then reflecting across the seam via
                                    the private ``cross5x6Interruption`` (not the
                                    "V2" the aperture-7 topology's own inline
                                    ``canonicalize5x6`` actually descends from —
                                    that one is a different, RI7H.ec-local inline
                                    block that happens to share the name; this is
                                    a fresh port of the real exported function,
                                    verified against the eC's own embedded
                                    ``cross5x6InterruptionTest`` vectors).

``planar_centroid`` does NOT call ``canonicalize5x6``/``move5x6_vertex_v1`` — the
eC ``centroid`` property itself returns the raw anchor + sub-hex offset with no
wrap step, and its callers (``getZoneWGS84Centroid``) pass that straight to the
projection's ``inverse`` uncanonicalized. Both helpers are ported here (module
helpers, per the Task 2 interface) because later I3H tasks (vertices, quantize)
need them; :func:`_i3h_anchor` is the one ``planar_centroid`` actually uses.

The *arithmetic* is held byte-for-byte faithful to the eC (same expressions,
groupings, epsilons, evaluation order): ``math.floor(...)`` where the eC writes
an explicit ``(int)floor(...)`` cast, ``math.trunc(...)`` where it writes a bare
``(int)(...)`` cast (truncation toward zero, which differs from floor on
negatives).
"""
from __future__ import annotations

import math
from typing import Any

from py4dggs.indexings.i3h import pack_i3h, unpack_i3h, _ilrc_from_lrti, _from_i9r, is_pentagon_cell, absolute_level
from py4dggs.types import NULL_ZONE, PlanarPoint


# --------------------------------------------------------------------------- #
# 5x6 vertex offset across a rhombus interruption (ri5x6.ec move5x6Vertex, v1)
# --------------------------------------------------------------------------- #
def move5x6_vertex_v1(cx, cy, dx, dy):
    """Offset a vertex from centroid ``(cx,cy)`` by ``(dx,dy)``, re-expressing the
    offset across a rhombus interruption ("dent") when the naive offset point
    crosses one (``move5x6Vertex`` v1, ``ri5x6.ec:1388-1434``). Unlike its ``v2``
    sibling (ported in ``hex_a7.move5x6_vertex``), this version never pre-crosses
    the *source* point — it only reflects the naive ``v = c + d`` result. The four
    ``vx</vy`` branches mirror the four crossing directions in the eC and
    intentionally share bodies (two identical pairs) — left as-is for
    faithfulness. ``icx``/``icy`` are truncating casts (bare ``(int)(...)`` in the
    eC), distinct from the ``ivx``/``ivy`` floor casts (explicit ``(int)floor``)."""
    icx = math.trunc(cx + 1e-11)
    icy = math.trunc(cy + 1e-11)

    vx = cx + dx
    vy = cy + dy

    sgn_dx = 0 if dx == 0 else (1 if dx > 0 else -1)   # eC Sgn(): Sgn(0)==0 (ri5x6.ec:1394)
    sgn_dy = 0 if dy == 0 else (1 if dy > 0 else -1)
    ivx = math.floor(cx + dx - sgn_dx * 1e-11)
    ivy = math.floor(cy + dy - sgn_dy * 1e-11)

    if (((ivx != icx and abs(vy - ivy) > 1e-11) or
         (ivy != icy and abs(vx - ivx) > 1e-11)) and
            (ivy - ivx > 1 or ivy < ivx)):
        if ivx < icx:
            # Stepping over bottom dent to the left
            vx = icx - (cy - icy) + dx - dy
            vy = icy + dx
        elif ivx > icx:
            # Stepping over top dent to the right
            vx = icx - (cy - icy) + dx - dy
            vy = icy + dx
        elif ivy < icy:
            # Stepping over top dent to the left
            vx = icx + dy
            vy = icy - (cx - icx) - dx + dy
        elif ivy > icy:
            # Stepping over bottom dent to the right
            vx = icx + dy
            vy = icy - (cx - icx) - dx + dy
    return [vx, vy]


# --------------------------------------------------------------------------- #
# 5x6 vertex offset, crossEarly variant (ri5x6.ec move5x6Vertex2, v2)
# --------------------------------------------------------------------------- #
def move5x6_vertex_v2(cx, cy, dx, dy, cross_early):
    """Offset a point from centroid ``(cx,cy)`` by ``(dx,dy)``, re-expressing the
    offset across a rhombus interruption when crossed — the ``crossEarly`` sibling
    of :func:`move5x6_vertex_v1` (``move5x6Vertex2``, ``ri5x6.ec:1470-1560``). This
    is the mover ``I3HZone::getNeighbor`` actually calls (``RI3H.ec:1146``), and is
    NOT the same as ``move5x6_vertex_v1`` (nor ``hex_a7.move5x6_vertex``, which is
    *structurally* the v1 naive-then-reflect body — eC ``move5x6Vertex``,
    ``ri5x6.ec:1388`` — though not byte-for-byte identical to this file's v1: those
    two differ in ``floor`` vs ``trunc`` and in ``Sgn(0)`` on integer-boundary ties):
    v2 additionally (a) pre-crosses the *centroid* itself when ``cross_early`` and
    the naive step already sits past a dent, (b) computes the crossing at the
    half-distance mid-edge via the ``pi1``/``pi2`` interruption points rather than
    the direct re-expression v1 uses, and (c) snaps a result landing on the polar
    diagonals to the "North" ``(1,0)`` / "South" ``(4,6)`` pole. ``icx/icy`` are
    ``(int)floor`` casts; ``icx2/icy2`` are the ``-e`` floor used to detect a
    centroid sitting exactly on an integer rhombus boundary."""
    e = 1e-11
    c_x, c_y = cx, cy
    icx = math.floor(c_x + e); icy = math.floor(c_y + e)
    icx2 = math.floor(c_x - e); icy2 = math.floor(c_y - e)
    nx = c_x + dx; ny = c_y + dy
    if nx < 0: nx += 5
    elif nx > 5: nx -= 5
    if ny < 0: ny += 5
    elif ny > 5 and c_y > 6 - e: ny -= 5
    n_top_right = (nx > c_x + e and nx - c_x < 3) or c_x - nx > 3
    n_top_left = (nx < c_x - e and c_x - nx < 3) or nx - c_x > 3
    n_bottom_right = (ny > c_y and ny - c_y < 3) or c_y - ny > 3  # eC: no epsilon on ny > c.y
    n_bottom_left = (ny < c_y - e and c_y - ny < 3) or ny - c_y > 3
    at_top_dent_r = icx2 != icx and c_x > c_y + e and n_top_right
    at_top_dent_l = icy2 != icy and c_x > c_y + e and n_top_left
    at_bottom_dent_l = icx2 != icx and c_y > c_x + 1 + e and n_bottom_left
    at_bottom_dent_r = icy2 != icy and c_y > c_x + 1 + e and n_bottom_right

    # Cross already for cases where crossing does not happen mid-edge
    if cross_early:
        if at_top_dent_r:
            c_x, c_y = icx + 1.0 - (c_y - icy), icy + 1
        elif at_top_dent_l:
            c_x, c_y = icx, icy - (c_x - icx)
        elif at_bottom_dent_l:
            c_x, c_y = icx - (c_y - icy), icy
        elif at_bottom_dent_r:
            c_x, c_y = icx + 1, icy + (icx + 1 - c_x)
        if c_x > 5 or c_y > 6 + e:
            c_x -= 5; c_y -= 5
        elif c_x < 0 or c_y < -e:
            c_x += 5; c_y += 5
        icx = math.floor(c_x + e); icy = math.floor(c_y + e)

    vx = c_x + dx; vy = c_y + dy
    ivx = math.floor(c_x + dx + 1e-11); ivy = math.floor(c_y + dy + 1e-11)
    if (((ivx != icx and abs(vy - ivy) > 1e-11) or (ivy != icy and abs(vx - ivx) > 1e-11)) and
            (ivy - ivx > 1 or ivy < ivx)):
        # Assuming the crossing point is at half the distance
        if abs(vx - vy - 1) < 1e-10:
            vx, vy = 1, 0  # "North" pole
        elif abs(vy - vx - 2) < 1e-10:
            vx, vy = 4, 6  # "South" pole
        elif ivx < icx and vx - ivx < 1 - e:
            # Stepping over bottom dent to the left
            pi1x, pi1y = icx, c_y + 0.5 * (icx - c_x)
            pi2x, pi2y = icx - (pi1y - icy), icy
            vx = pi2x + pi1y - c_y; vy = pi2y + pi1x - c_x
        elif ivx > icx and vx - ivx > e:
            # Stepping over top dent to the right
            pi1x, pi1y = icx + 1, c_y + 0.5 * (icx + 1 - c_x)
            pi2x, pi2y = icx + 1 + (icy + 1 - pi1y), icy + 1
            vx = pi2x + pi1y - c_y; vy = pi2y + pi1x - c_x
        elif ivy < icy and vy - ivy < 1 - e:
            # Stepping over top dent to the left
            pi1x, pi1y = c_x + 0.5 * (icy - c_y), icy
            pi2x, pi2y = icx, icy - (pi1x - icx)
            vx = pi2x + pi1y - c_y; vy = pi2y + pi1x - c_x
        elif ivy > icy and vy - ivy > e:
            # Stepping over bottom dent to the right
            pi1x, pi1y = vx + 0.5 * (icy + 1 - c_y), icy + 1
            pi2x, pi2y = icx + 1, icy + 1 + (icx + 1 - pi1x)
            vx = pi2x + pi1y - c_y; vy = pi2y + pi1x - c_x
    if vx > 5:
        vx -= 5; vy -= 5
    elif vx < 0:
        vx += 5; vy += 5
    return [vx, vy]


# --------------------------------------------------------------------------- #
# 5x6 interruption-crossing test + reflecting mover (ri5x6.ec:1128-1155, 1304-1374, 1438-1463)
# --------------------------------------------------------------------------- #
_INTERRUPTIONS_5X6 = [
    [  # North (h=0), r=0..4, each [left, right]
        [(0, 0), (1, 0)], [(1, 0), (1, 1)],
        [(1, 1), (2, 1)], [(2, 1), (2, 2)],
        [(2, 2), (3, 2)], [(3, 2), (3, 3)],
        [(3, 3), (4, 3)], [(4, 3), (4, 4)],
        [(4, 4), (5, 4)], [(5, 4), (5, 5)],
    ],
    [  # South (h=1)
        [(0, 1), (0, 2)], [(0, 2), (1, 2)],
        [(1, 2), (1, 3)], [(1, 3), (2, 3)],
        [(2, 3), (2, 4)], [(2, 4), (3, 4)],
        [(3, 4), (3, 5)], [(3, 5), (4, 5)],
        [(4, 5), (4, 6)], [(4, 6), (5, 6)],
    ],
]
_INTERRUPTIONS_5X6 = [
    [[_INTERRUPTIONS_5X6[h][2 * r], _INTERRUPTIONS_5X6[h][2 * r + 1]] for r in range(5)]
    for h in range(2)
]


def _intersects5x6_interruption(a0x, a0y, a1x, a1y, b0x, b0y, b1x, b1y):
    """``intersects5x6Interruption`` (``ri5x6.ec:1128-1155``): segment-segment
    intersection of ``a0->a1`` against a fixed interruption edge ``b0->b1``,
    returning the intersection point and the parametric distance ``t`` along
    ``a0->a1`` (used to pick the *nearest* crossing)."""
    e = 1e-12
    s1x = a1x - a0x
    s1y = a1y - a0y
    s2x = b1x - b0x
    s2y = b1y - b0y
    dx = a0x - b0x
    dy = a0y - b0y
    d = s1x * s2y - s2x * s1y
    if abs(d) > 1e-13:
        factor = 1.0 / d
        s = (s1x * dy - s1y * dx) * factor
        if s - e >= 0 and s + e <= 1:
            t = (s2x * dy - s2y * dx) * factor
            if t - e >= 0 and t + e <= 1:
                return True, a0x + t * s1x, a0y + t * s1y, t
    return False, None, None, None


def _crosses5x6_interruption(cx, cy, dx, dy):
    """``crosses5x6Interruption`` (``ri5x6.ec:1304-1374``): does the segment
    ``(cx,cy) -> (cx+dx,cy+dy)`` cross one of the 10 fixed interruption edges?
    Returns ``(crossed, (src_x,src_y), (dst_x,dst_y), north)`` where ``src`` is
    the crossing point and ``dst`` is its mirrored far-side representation
    (via the already-ported :func:`_cross5x6_interruption`)."""
    def snap(v):
        for k in range(7):
            if abs(v - k) < 1e-12:
                return float(k)
        return v

    c_x, c_y = snap(cx), snap(cy)
    if math.trunc(c_x + 1) == math.trunc(c_x + dx + 1) and math.trunc(c_y + 1) == math.trunc(c_y + dy + 1):
        return False, None, None, None

    min_t = float("inf")
    src = None
    cross_h = cross_s = None
    found = False
    for h in range(2):
        for r in range(5):
            for s in range(2):
                (bx0, by0), (bx1, by1) = _INTERRUPTIONS_5X6[h][r][s]
                ok, ix, iy, t = _intersects5x6_interruption(c_x, c_y, c_x + dx, c_y + dy, bx0, by0, bx1, by1)
                if ok and t < min_t:
                    src = (ix, iy)
                    min_t = t
                    cross_h, cross_s = h, s
                    found = True
    if not found:
        return False, None, None, None
    dst = _cross5x6_interruption(src[0], src[1], cross_h == 1, cross_s == 0)
    north = cross_h == 0
    return True, src, dst, north


def move5x6_vertex_v3(cx, cy, dx, dy):
    """``move5x6Vertex3`` (``ri5x6.ec:1438-1463``): offset a point from
    ``(cx,cy)`` by ``(dx,dy)``, reflecting across a rhombus interruption via
    the *nearest-crossing* test (:func:`_crosses5x6_interruption`), distinct
    from both :func:`move5x6_vertex_v1` (naive-then-reflect) and
    ``move5x6_vertex_v2`` (crossEarly, ``ri5x6.ec:1470-1560``). This is the
    mover ``I3HSubZones.ec``'s sub-zone generators call directly (never v2)."""
    crossed, i1, i2, north = _crosses5x6_interruption(cx, cy, dx, dy)
    if crossed:
        i1x, i1y = i1
        i2x, i2y = i2
        if north:
            vx = i2x - 2 * (dy - (i1y - cy))
            vy = i2y + dx - (i1x - cx)
        else:
            vx = i2x + dy - (i1y - cy)
            vy = i2y + 2 * (dx - (i1x - cx))
    else:
        vx = cx + dx
        vy = cy + dy
    if vx > 5 and vy > 5:
        vx -= 5
        vy -= 5
    elif vx < 0 and vy < 1:
        vx += 5
        vy += 5
    return [vx, vy]


# --------------------------------------------------------------------------- #
# Fundamental-domain folding (ri5x6.ec canonicalize5x6 + private cross5x6Interruption)
# --------------------------------------------------------------------------- #
def _cross5x6_interruption(src_x, src_y, south, left):
    """Reflect a point sitting exactly on a rhombus interruption seam to its
    other-side representation (the private, non-"V2" ``cross5x6Interruption``,
    ``ri5x6.ec:1247-1299``). ``south``/``left`` select which of the four seam
    crossings; :func:`canonicalize5x6` only ever calls this with ``left=False``,
    but both are ported for faithfulness. Verified against the eC's own embedded
    ``cross5x6InterruptionTest`` vectors (``ri5x6.ec:1166-1229``, ``TEST_CROSSING``)."""
    if not south:
        if left:
            ix = math.trunc(src_y + 1e-11)
            dst_x, dst_y = src_y, ix - (src_x - ix)
        else:
            iy = math.trunc(src_x - 1 + 1e-11)
            dst_x, dst_y = iy + 2 - (src_y - iy), src_x
    else:
        if left:
            iy = math.trunc(src_x + 1 + 1e-11)
            dst_x, dst_y = iy - 1 - (src_y - iy), src_x + 1
        else:
            ix = math.trunc(src_y - 2 + 1e-11)
            dst_x, dst_y = src_y - 1, ix + 3 - (src_x - ix)

    if dst_x > 5 - 1e-11 and dst_y > 5 - 1e-11:
        dst_x -= 5; dst_y -= 5
    elif dst_x < 1e-11 and dst_y < 1 - 1e-11:
        dst_x += 5; dst_y += 5
    return [dst_x, dst_y]


def canonicalize5x6(cx, cy):
    """Fold a 5x6 planar point into the canonical fundamental domain: wrap the
    two out-of-range corners, detect the polar diagonal points and the rhombus
    interruption seams from the ``(int)floor``-cast row/col, then reflect a
    seam point via :func:`_cross5x6_interruption` (``canonicalize5x6``,
    ``ri5x6.ec:1562-1605``). This is a fresh port of the eC's *exported*
    function — distinct from ``hex_a7.canonicalize5x6``, which (despite the
    shared name) is actually an inlined block from ``RI7H.ec``'s
    ``I7HZone::fromCentroid`` and uses a different algorithm entirely; the two
    were compared line-by-line and are not the same port."""
    x, y = cx, cy
    if x > 5 - 1e-11 and y > 5 - 1e-11:
        x -= 5; y -= 5
    if x < -1e-11 or y < -1e-11:
        x += 5; y += 5

    icx = math.floor(x + 1e-11)
    icy = math.floor(y + 1e-11)

    cross = False
    south = False
    np_ = False
    sp = False

    if icy == 0:
        cross = abs(x - 1) < 1e-11
        np_ = cross and abs(y - 0) < 1e-11
    elif icy == 1:
        cross = abs(x - 2) < 1e-11
        np_ = cross and abs(y - 1) < 1e-11
    elif icy == 2:
        cross = abs(x - 3) < 1e-11
        np_ = cross and abs(y - 2) < 1e-11
    elif icy == 3:
        cross = abs(x - 4) < 1e-11
        np_ = cross and abs(y - 3) < 1e-11
    elif icy == 4:
        cross = abs(x - 5) < 1e-11
        np_ = cross and abs(y - 4) < 1e-11

    if icx == 0:
        if abs(y - 2) < 1e-11:
            cross = True; south = True; sp = abs(x - 0) < 1e-11
    elif icx == 1:
        if abs(y - 3) < 1e-11:
            cross = True; south = True; sp = abs(x - 1) < 1e-11
    elif icx == 2:
        if abs(y - 4) < 1e-11:
            cross = True; south = True; sp = abs(x - 2) < 1e-11
    elif icx == 3:
        if abs(y - 5) < 1e-11:
            cross = True; south = True; sp = abs(x - 3) < 1e-11
    elif icx == 4:
        if abs(y - 6) < 1e-11:
            cross = True; south = True; sp = abs(x - 4) < 1e-11

    if sp:
        return [4, 6]
    if np_:
        return [1, 0]
    if cross:
        return _cross5x6_interruption(x, y, south, False)
    if abs(x - 5) < 1e-11:
        return [0, y - 5]
    return [x, y]


# --------------------------------------------------------------------------- #
# I3H rhombus anchor (RI3H.ec centroid property, non-polar branch)
# --------------------------------------------------------------------------- #
def _i3h_anchor(level_i9r, root, rix):
    """The rhombus anchor (top-left corner, ``tl``) in 5x6-fractional units for a
    non-polar I3H cell (``I3HZone::centroid``, ``RI3H.ec:2150-2158``). ``rix``
    (``rhombusIX``) is decomposed into a row/col pair via ``rowOP``/``colOP``
    (root's row/col in units of the aperture-3 level's cell count ``p``) plus
    ``ixOP = rix // p`` (the eC's comment: "distributivity on: rix - ixOP*p for
    rix % p" — i.e. ``col = colOP*p + (rix mod p)``, computed here the same way
    the eC does, without the intermediate ``uint64`` wraparound the eC relies on
    being mathematically inert)."""
    p = 3 ** level_i9r
    row_op = (root + 1) >> 1
    col_op = root >> 1
    ix_op = rix // p
    row = row_op * p + ix_op
    col = (col_op - ix_op) * p + rix
    d = 1.0 / p
    return col * d, row * d


def _i3h_centroid_xy(level_i9r, root, rix, sub_hex):
    """The cell centroid in raw 5x6-fractional units (``I3HZone::centroid``,
    ``RI3H.ec:2138-2168``): the two polar shortcuts (root 10/11) then the
    :func:`_i3h_anchor` top-left corner plus the C/D sub-hex offsets. Shared by
    :meth:`HexAperture3Topology.planar_centroid` and the neighbour port
    (:func:`_i3h_get_neighbor`, which needs the same raw ``(x, y)`` the eC's
    ``getNeighbor`` reads from ``this.centroid``)."""
    if root == 10:
        return 1.0, 0.0  # "North" pole
    if root == 11:
        return 4.0, 6.0  # "South" pole
    tlx, tly = _i3h_anchor(level_i9r, root, rix)
    d = 1.0 / (3 ** level_i9r)
    if sub_hex == 0 or sub_hex == 1:
        return tlx, tly  # Even level A or Odd level B hex
    if sub_hex == 2:
        return tlx + 2 * d / 3, tly + d / 3  # Odd level C hex
    return tlx + d / 3, tly + 2 * d / 3  # Odd level D hex


# --------------------------------------------------------------------------- #
# fromCentroid: quantize a 5x6 point to an I3H zone (RI3H.ec I3HZone::fromCentroid,
# RI3H.ec:1330-1672) — closed-form (NOT RI7H's round-to-nearest + candidate-parent
# search).
# --------------------------------------------------------------------------- #
# `NULL_ZONE` (dggrs.ec:9 `nullZone`, imported from py4dggs.types) is what the
# out-of-range guard below returns. NOTE: it IS reachable from ordinary valid
# lat/lon — within ~0.1 deg of a pole at the coarsest resolutions
# (level_i9r == 0), on the order of 1 in 10^4 near-pole points. DGGAL itself
# agrees there (same sentinel, same meaningless centroid), so the arithmetic is
# faithful; what this port must not do is dress the sentinel up as an ordinary
# zone in the text id. See `I3HIndexing.to_text`.


def _i3h_from_centroid(level, cx, cy):
    """Quantize a 5x6 planar point to the ``(levelI9R, rootRhombus, rhombusIX,
    subHex)`` tuple for an I3H zone at ``level`` (``I3HZone::fromCentroid``,
    ``RI3H.ec:1330-1672``). Returns ``None`` for the eC's ``nullZone`` guard
    (malformed input; not expected for valid WGS84 points).

    Three parts, ported byte-for-byte (same expressions/groupings/epsilons/
    evaluation order as the eC):

    1. Pole/dent/wrap preamble (``RI3H.ec:1337-1363``) — structurally the same
       preamble as ``hex_a7.from_centroid``'s call to its own ``canonicalize5x6``
       (``RI7H.ec:1439-1462``), EXCEPT the final out-of-range branch calls
       ``move5x6Vertex`` (our :func:`move5x6_vertex_v1`, invoked here as
       ``move5x6_vertex_v1(5, 5, c_x, c_y)`` — the eC's
       ``move5x6Vertex(c, {5,5}, c.x, c.y)`` aliases its output param ``v`` to
       ``c``, using ``c``'s own (pre-call) ``x, y`` as the offset ``dx, dy``
       from the fixed center ``(5, 5)``) rather than a plain ``+= 5``. This is
       NOT the module's own :func:`canonicalize5x6` (a different, later-stage
       wrap used by ``planar_vertices``) — it tracks north/south-pole-diagonal
       booleans instead of returning early, because the pole cases still need a
       ``root``/``rix`` (10/0 or 11/0) and a level-parity ``sub_hex`` below.
       The trailing bottom-right wrap (``RI3H.ec:1361-1363``) mirrors the extra
       fix-up ``hex_a7.from_centroid`` applies just after its own
       ``canonicalize5x6`` call.

    2. Floor into the aperture-9 (3^levelI9R x 3^levelI9R) cell within the root
       rhombus (``RI3H.ec:1365-1417``): truncating (bare ``(int)``/``(uint64)``
       casts, not ``floor``) into ``cx, cy`` (root-rhombus row/col) and
       ``x, y`` (within-rhombus row/col), with three edge fix-ups for the rare
       ``x == p`` / ``y == p`` / dent-corner (``cy - cx > 1``, ``cy < cx``)
       cases the straightforward floor can land on. ``rix``/``xd``/``yd`` are
       computed from the FINAL (post-fix-up) ``x, y, cx, cy``.

    3. Fractional-thirds sub-hex classification (``RI3H.ec:1419-1668``): the
       pole booleans short-circuit to ``root=10/11, rix=0``; otherwise odd
       ``level`` (``RI3H.ec:1434-1550``) classifies B (centre third) vs the
       rhombus-crossing/non-polar-pentagon B reassignments vs C/D, and even
       ``level`` (``RI3H.ec:1552-1667``) classifies A with the same
       rhombus-crossing/pentagon reassignments (topRight/bottomLeft/bottomRight
       portions of the hex, per the ``xd/yd`` half-plane tests)."""
    l9r = level // 2
    p = 3 ** l9r
    d = 1.0 / p

    c_x, c_y = cx, cy
    is_north_pole = False
    is_south_pole = False

    if abs(c_x - c_y - 1) < 1e-10:
        is_north_pole = True
    elif abs(c_y - c_x - 2) < 1e-10:
        is_south_pole = True
    elif c_y < -1e-11 and c_x > -1e-11:
        c_x -= c_y; c_y = 0
    elif math.floor(c_x + 1e-11) > math.floor(c_y + 1e-11):
        # Over top dent to the right
        iy = min(5, math.floor(c_y + 1e-11))
        c_x += (iy + 1 - c_y); c_y = iy + 1
    elif math.floor(c_y + 1e-11) - math.floor(c_x + 1e-11) > 1:
        # Over bottom dent to the right
        ix = min(4, math.floor(c_x + 1e-11))
        c_y += (ix + 1 - c_x); c_x = ix + 1
    elif c_x < -1e-11 or c_y < -1e-11:
        c_x, c_y = move5x6_vertex_v1(5, 5, c_x, c_y)

    if (c_x > 5 - 1e-11 and c_y > 5 - 1e-11 and  # bottom-right wrap, e.g. A9-0E/A9-0F
            c_x + c_y > 5.0 + 5.0 - d - 1e-11):
        c_x -= 5; c_y -= 5

    cx = min(4, math.trunc(c_x + 1e-11))   # Coordinate of root rhombus
    cy = min(5, math.trunc(c_y + 1e-11))
    root = cx + cy
    x = math.trunc((c_x - cx) * p + 1e-6)  # Row and column within root rhombus
    y = math.trunc((c_y - cy) * p + 1e-6)

    # REVIEW (eC): x == p or y == p are currently possible, yet the code below
    # assumed it was not.
    if y == p:  # IVEA3H B1-2-A
        cy += 1
        root += 1
        y -= p
        c_y = cy + y / p
    if x == p:  # IVEA3H B4-6-A
        cx += 1
        root += 1
        if root == 10:
            cx -= 5
            cy -= 5
            root = 0
            c_y = cy + y / p
        x -= p
        c_x = cx + x / p
    if cy - cx > 1 and not y:  # IVEA3H B9-3-A, C9-12-A
        cx += 1
        y = p - x
        x = 0
        root += 1
        c_y = cy + y / p
        c_x = cx
    elif cy < cx and not x:  # RTEA3H B4-1-A
        cy += 1
        x = p - y
        y = 0
        root += 1
        c_x = cx + x / p
        c_y = cy

    rix = y * p + x  # Index within root rhombus
    xd = (c_x - cx) * p - x
    yd = (c_y - cy) * p - y

    if is_north_pole:
        sh = 1 if (level & 1) else 0
        root, rix = 10, 0
    elif is_south_pole:
        sh = 1 if (level & 1) else 0
        root, rix = 11, 0
    else:
        right_sr = (x == p - 1); top_sr = (y == 0)
        left_sr = (x == 0); bottom_sr = (y == p - 1)
        np_sub_rhombus = right_sr and top_sr and not (root & 1)
        sp_sub_rhombus = bottom_sr and left_sr and (root & 1)

        if cy < cx or xd < -1 or yd < -1 or x >= p or y >= p or rix >= p * p:
            return None  # nullZone -- y cannot be smaller than x

        if level & 1:  # Odd level
            left_third = 3 * xd < 1; top_third = 3 * yd < 1
            if left_third and top_third:
                sh = 1  # B
            else:
                right_third = 3 * xd > 2; bottom_third = 3 * yd > 2
                if right_third and bottom_third:
                    if bottom_sr and right_sr:
                        # Non-polar pentagon
                        root = (root + 2) % 10; rix = 0
                    elif bottom_sr:
                        # Indexed to another root rhombus
                        if root & 1:
                            # Crossing South interruption to the right
                            root = (root + 2) % 10; rix = (p - 1 - x) * p
                        else:
                            root += 1; rix = x + 1
                    elif right_sr:
                        # Indexed to another root rhombus
                        if not (root & 1):
                            # Crossing North interruption to the right
                            root = (root + 2) % 10; rix = p - 1 - y
                        else:
                            root = (root + 1) % 10; rix = p * (y + 1)
                    else:
                        rix += p + 1
                    sh = 1  # B
                elif bottom_third:
                    if 3 * (yd - xd) > 2:
                        if sp_sub_rhombus:
                            root, rix, sh = 11, 0, 1  # "South" pole B
                        else:
                            if bottom_sr:
                                # Indexed to another root rhombus
                                if root & 1:
                                    # Crossing South interruption to the right
                                    root = (root + 2) % 10; rix = (p - x) * p
                                else:
                                    rix = x; root += 1
                            else:
                                rix += p
                            sh = 1  # B
                    else:
                        sh = 3  # D
                elif right_third:
                    if 3 * (xd - yd) > 2:
                        if np_sub_rhombus:
                            root, rix, sh = 10, 0, 1  # "North" pole B
                        else:
                            if right_sr:
                                if not (root & 1):
                                    # Crossing North interruption to the right
                                    root = (root + 2) % 10; rix = p - y
                                else:
                                    root = (root + 1) % 10; rix = p * y
                            else:
                                rix += 1
                            sh = 1  # B
                    else:
                        sh = 2  # C
                elif xd > yd:
                    sh = 2  # C
                else:
                    sh = 3  # D
        else:  # Even level
            top_right = bottom_left = bottom_right = False
            if xd - 1 > -yd:  # Bottom-Right portion
                if xd > yd * 2:
                    top_right = True  # Top-right hexagon
                elif 2 * xd < yd:
                    bottom_left = True  # Bottom-left hexagon
                else:
                    bottom_right = True  # Bottom-right hexagon
            else:  # Top-Left portion
                if 2 * xd > yd + 1:
                    top_right = True  # Top-right hexagon
                elif xd + 1 < 2 * yd:
                    bottom_left = True  # Bottom-left hexagon

            sh = 0  # A
            if top_right:
                if np_sub_rhombus:
                    root, rix, sh = 10, 0, 0  # "North" pole A
                else:
                    if right_sr:
                        if not (root & 1):
                            # Crossing North interruption to the right
                            root = (root + 2) % 10; rix = p - y
                        else:
                            root = (root + 1) % 10; rix = p * y
                    else:
                        rix += 1
            elif bottom_left:
                if sp_sub_rhombus:
                    root, rix, sh = 11, 0, 0  # "South" pole A
                else:
                    if bottom_sr:
                        # Indexed to another root rhombus
                        if root & 1:
                            # Crossing South interruption to the right
                            root = (root + 2) % 10; rix = (p - x) * p
                        else:
                            rix = x; root += 1
                    else:
                        rix += p
            elif bottom_right:
                if bottom_sr and right_sr:
                    # Non-polar pentagon
                    root = (root + 2) % 10; rix = 0
                elif bottom_sr:
                    # Indexed to another root rhombus
                    if root & 1:
                        # Crossing South interruption to the right
                        root = (root + 2) % 10; rix = (p - 1 - x) * p
                    else:
                        root += 1; rix = x + 1
                elif right_sr:
                    # Indexed to another root rhombus
                    if not (root & 1):
                        # Crossing North interruption to the right
                        root = (root + 2) % 10; rix = p - 1 - y
                    else:
                        root = (root + 1) % 10; rix = p * (y + 1)
                else:
                    rix += p + 1

    return l9r, root, rix, sh


# --------------------------------------------------------------------------- #
# I3H boundary vertices (RI3H.ec I3HZone::getVertices, corner-anchored)
# --------------------------------------------------------------------------- #
def _i3h_vertices(level_i9r, root, rix, sub_hex):
    """Raw (uncanonicalized) 5x6 boundary vertices for an I3H cell
    (``I3HZone::getVertices``, ``RI3H.ec:1674-1790``). Corner-anchored: recomputes
    the same ``tl`` anchor as :func:`_i3h_anchor` for non-polar roots (the
    ``else`` branch below reuses it verbatim), but the pole roots (10/11) get
    their OWN ``tl`` — the eC's ``row``/``col`` ternary at the top of
    ``getVertices`` gives ``(col, row) = (p-1, 0)`` for the "North" pole and
    ``(4p, 6p-1)`` for the "South" pole. This is NOT the same point
    :func:`_i3h_anchor` would produce if handed ``root=10``/``11`` (that helper's
    ``rowOP``/``colOP`` formula is only valid for roots 0-9), nor the same as
    ``planar_centroid``'s hardcoded pole shortcuts ``(1,0)``/``(4,6)`` (those are
    *centroid* shortcuts; the vertex anchor here is offset from them by ``d``).

    Poles (root 10/11) emit an explicit 5-vertex fan: one ``move5x6_vertex_v1``
    call seeds a single vertex ``v``, then four more vertices are formed by
    walking the 5x6 diagonal in whole-unit steps from ``v`` (``v.x+k, v.y+k`` —
    NOT further ``move5x6Vertex`` calls), mirroring the five rhombi that meet at
    each pole. The ``sub_hex==0``/"South" branch's last step is ``+1`` rather
    than the ``-4`` its ``-0,-1,-2,-3`` siblings would suggest (RI3H.ec:1713-
    1717) — preserved verbatim (``v.x-4``/``v.x+1`` are the same point once
    5x6-wrapped, but the eC literally writes ``+1``, so this port does too).

    Regular (non-polar) cells use the fixed ``d/3``-unit offset tables A
    (``sub_hex==0``, even level), B (``sub_hex==1``, odd level), C
    (``sub_hex==2``), D (``sub_hex==3``) from ``tl``. Only A and B ever touch a
    rhombus corner, so only they suppress a vertex for the ``rix==0`` pentagon
    cells — dropping exactly one of their two corner-adjacent vertices depending
    on ``south = root & 1`` (RI3H.ec:1723-1727, 1762-1766). C and D hexes are
    interior to the rhombus and always emit all 6 vertices.

    Returns a list of ``[x, y]`` pairs; the caller canonicalizes each one
    (mirrors ``getZoneCRSVertices``/``getZoneWGS84Vertices``, RI3H.ec:337-386,
    which call ``canonicalize5x6`` on every raw vertex from ``getVertices``
    before use — ``getVertices`` itself never wraps)."""
    p = 3 ** level_i9r
    d = 1.0 / p

    if root == 10:
        tlx, tly = (p - 1) * d, 0.0
    elif root == 11:
        tlx, tly = 4 * p * d, (6 * p - 1) * d
    else:
        tlx, tly = _i3h_anchor(level_i9r, root, rix)

    south = bool(root & 1)
    verts: list[list[float]] = []

    if sub_hex == 0:  # Even level
        if root == 10:  # "North" pole
            vx, vy = move5x6_vertex_v1(tlx, tly, d / 3, -d / 3)
            verts = [[vx + k, vy + k] for k in (1, 2, 3, 4, 5)]
        elif root == 11:  # "South" pole
            vx, vy = move5x6_vertex_v1(tlx, tly, -d / 3, d / 3)
            verts = [[vx - 0, vy - 0], [vx - 1, vy - 1], [vx - 2, vy - 2], [vx - 3, vy - 3], [vx + 1, vy + 1]]
        else:  # Regular case
            verts.append(move5x6_vertex_v1(tlx, tly, 2 * d / 3, d / 3))
            verts.append(move5x6_vertex_v1(tlx, tly, d / 3, 2 * d / 3))
            if not south or rix:  # 0 rhombusIndex are pentagons
                verts.append(move5x6_vertex_v1(tlx, tly, -d / 3, d / 3))
            verts.append(move5x6_vertex_v1(tlx, tly, -2 * d / 3, -d / 3))
            if south or rix:  # 0 rhombusIndex are pentagons
                verts.append(move5x6_vertex_v1(tlx, tly, -d / 3, -2 * d / 3))
            verts.append(move5x6_vertex_v1(tlx, tly, d / 3, -d / 3))
    elif sub_hex == 1:  # Odd level -- type B
        if root == 10:  # "North" pole
            vx, vy = move5x6_vertex_v1(tlx, tly, 2 * d / 3, 0)
            verts = [[vx + k, vy + k] for k in (1, 2, 3, 4, 5)]
        elif root == 11:  # "South" pole
            vx, vy = move5x6_vertex_v1(tlx, tly, d / 3, d)
            verts = [[vx - k, vy - k] for k in (0, 1, 2, 3, 4)]
        else:
            if south or rix:  # 0 rhombusIndex are pentagons
                verts.append(move5x6_vertex_v1(tlx, tly, d / 3, 0))
            verts.append(move5x6_vertex_v1(tlx, tly, d / 3, d / 3))
            verts.append(move5x6_vertex_v1(tlx, tly, 0, d / 3))
            if not south or rix:  # 0 rhombusIndex are pentagons
                verts.append(move5x6_vertex_v1(tlx, tly, -d / 3, 0))
            verts.append(move5x6_vertex_v1(tlx, tly, -d / 3, -d / 3))
            verts.append(move5x6_vertex_v1(tlx, tly, 0, -d / 3))
    elif sub_hex == 2:  # Odd level -- type C
        verts = [
            move5x6_vertex_v1(tlx, tly, d / 3, 0),
            move5x6_vertex_v1(tlx, tly, 2 * d / 3, 0),
            move5x6_vertex_v1(tlx, tly, d, d / 3),
            move5x6_vertex_v1(tlx, tly, d, 2 * d / 3),
            move5x6_vertex_v1(tlx, tly, 2 * d / 3, 2 * d / 3),
            move5x6_vertex_v1(tlx, tly, d / 3, d / 3),
        ]
    else:  # sub_hex == 3, Odd level -- type D
        verts = [
            move5x6_vertex_v1(tlx, tly, 0, d / 3),
            move5x6_vertex_v1(tlx, tly, d / 3, d / 3),
            move5x6_vertex_v1(tlx, tly, 2 * d / 3, 2 * d / 3),
            move5x6_vertex_v1(tlx, tly, 2 * d / 3, d),
            move5x6_vertex_v1(tlx, tly, d / 3, d),
            move5x6_vertex_v1(tlx, tly, 0, 2 * d / 3),
        ]

    return verts


# --------------------------------------------------------------------------- #
# Exact aperture-3 neighbours (RI3H.ec I3HZone::getNeighbor / getNeighbors)
# --------------------------------------------------------------------------- #
# I3HNeighbor enum order (RI3H.ec:841-852) — getNeighbors iterates it in order,
# and its dedup/relabel logic depends on these exact values.
_TOP, _BOTTOM, _LEFT, _RIGHT, _TOP_LEFT, _TOP_RIGHT, _BOTTOM_LEFT, _BOTTOM_RIGHT = range(8)


def _i3h_get_neighbor(level_i9r, root, rix, sub_hex, which):
    """The I3H neighbour in one ``I3HNeighbor`` direction, or ``None`` for
    ``nullZone`` (``I3HZone::getNeighbor``, ``RI3H.ec:1005-1161``). A geometric
    planar method (not row/col arithmetic): pick a 5x6 offset ``(x, y)`` in units
    of ``d = 1/3^levelI9R`` from the cell centroid — with the even-level
    (``RI3H.ec:1021-1078``) and odd-level (``RI3H.ec:1080-1140``) offset tables
    and their pole / interruption special-cases — then apply it via
    :func:`move5x6_vertex_v2` (honouring ``crossEarly``) and re-quantize with
    :func:`_i3h_from_centroid` at the same resolution ``2*l9r + (sub_hex>0)``. The
    special-cases (poles, ``south``/``north`` interruptions) are precisely what
    the grid-agnostic edge k-ring lacks, hence its ~0.01% aperture-3 overshoot."""
    cxf, cyf = _i3h_centroid_xy(level_i9r, root, rix, sub_hex)  # this.centroid
    cx = math.floor(cxf + 1e-11)
    cy = math.floor(cyf + 1e-11)
    south = (cyf - cxf - 1e-11) > 1  # Not counting pentagons as south or north
    north = (cxf - cyf - 1e-11) > 0
    north_pole = north and abs(cxf - cyf - 1.0) < 1e-11
    south_pole = south and abs(cyf - cxf - 2.0) < 1e-11
    l9r = level_i9r
    p = 3 ** l9r
    d = 1.0 / p
    x = 0.0
    y = 0.0
    sh = sub_hex
    cross_early = True

    if sh == 0:  # Even level
        if which == _TOP:
            if south and cxf - cx < 1e-11:
                cross_early = False
                if south_pole:
                    x, y = -3, -3 - d
                else:  # Extra top neighbor at south interruptions
                    y = -d
        elif which == _BOTTOM:
            if north and cyf - cy < 1e-11:
                cross_early = False
                if north_pole:
                    x, y = 2 - d, 2
                else:  # Extra bottom neighbor at north interruptions
                    x = -d
        elif which == _LEFT:
            x, y = -d, -d
        elif which == _RIGHT:
            x, y = d, d
        elif which == _TOP_LEFT:
            if north_pole:
                cross_early = False; x, y = 3 - d, 3
            elif south_pole:
                cross_early = False; y = -d
            else:
                y = -d
        elif which == _BOTTOM_LEFT:
            if south_pole:
                cross_early = False; x, y = -2, -2 - d
            else:
                x = -d
        elif which == _TOP_RIGHT:
            if north_pole:
                cross_early = False; x, y = 4 - d, 4
            elif south_pole:
                cross_early = False; x, y = -4, -d - 4
            else:
                x = d
        elif which == _BOTTOM_RIGHT:
            if south_pole:
                cross_early = False; x, y = -1, -1 - d
            else:
                y = d
    else:  # Odd level
        do3 = d / 3
        if which == _TOP:
            if south_pole:
                x, y, cross_early = do3 - 5, -do3 - 5, False
            elif not north_pole:
                x, y = do3, -do3
        elif which == _BOTTOM:
            if north_pole:
                x, y, cross_early = 1 - do3, 1 + do3, False
            elif not south_pole:
                x, y = -do3, do3
        elif which == _TOP_LEFT:
            if north_pole:
                x, y, cross_early = 2 - do3, 2 + do3, False
            elif south_pole:
                x, y = do3, -do3
            else:
                x, y = -do3, -2 * do3
        elif which == _BOTTOM_LEFT:
            if north_pole:
                x, y, cross_early = 4 - do3, 4 + do3, False
            elif south_pole:
                x, y = do3 - 2, -do3 - 2
            else:
                x, y = -2 * do3, -do3
        elif which == _TOP_RIGHT:
            if north_pole:
                x, y, cross_early = 3 - do3, 3 + do3, False
            elif south_pole:
                x, y = do3 - 4, -do3 - 4
            else:
                x, y = 2 * do3, do3
        elif which == _BOTTOM_RIGHT:
            if north_pole:
                x, y, cross_early = 5 - do3, 5 + do3, False
            elif south_pole:
                x, y = do3 - 1, -do3 - 1
            else:
                x, y = do3, 2 * do3
        elif which == _RIGHT:  # Stand-in for second bottom / top neighbor at interruptions
            if north and not north_pole and cyf - cy < 1e-11:
                cross_early = False; y = do3; x = -do3
            elif south and not south_pole and cxf - cx < 1e-11:
                cross_early = False; x = do3; y = -do3

    if x or y:
        vx, vy = move5x6_vertex_v2(cxf, cyf, x, y, cross_early)
        res = _i3h_from_centroid(2 * l9r + (1 if sh > 0 else 0), vx, vy)
        if res is None:
            return None
        packed = pack_i3h(*res)
        if packed == pack_i3h(level_i9r, root, rix, sub_hex):
            return None  # result == this; should not happen
        return packed
    return None


def _i3h_get_neighbors(level_i9r, root, rix, sub_hex):
    """All I3H neighbours (5 for a pentagon, 6 for a hexagon), as packed ints
    (``I3HZone::getNeighbors``, ``RI3H.ec:1163-1201``). Iterates the eight
    ``I3HNeighbor`` directions in enum order, dropping ``nullZone`` results, and
    applies the eC's relabel/dedup: a ``topRight`` that equals the preceding
    ``topLeft`` collapses to a single ``top`` (``bottomRight``/``bottomLeft`` ->
    ``bottom`` likewise); a ``topRight``/``bottomRight`` following a different
    direction relabels to ``top``/``bottom``. The direction labels found here are
    what A1's ``getParents`` will consume; this function returns only the zone
    list, which is what :meth:`HexAperture3Topology.neighbors` needs."""
    neighbors = []
    directions = []
    for n in range(8):  # I3HNeighbor::enumSize
        nb = _i3h_get_neighbor(level_i9r, root, rix, sub_hex, n)
        if nb is None:
            continue
        which = n
        if neighbors:
            prev = directions[-1]
            if n == _TOP_RIGHT and prev == _TOP_LEFT and neighbors[-1] == nb:
                directions[-1] = _TOP
                continue
            elif n == _BOTTOM_RIGHT and prev == _BOTTOM_LEFT and neighbors[-1] == nb:
                directions[-1] = _BOTTOM
                continue
            elif n == _TOP_RIGHT and prev != _TOP_LEFT:
                which = _TOP
            elif n == _BOTTOM_RIGHT and prev != _BOTTOM_LEFT:
                which = _BOTTOM
        directions.append(which)
        neighbors.append(nb)
    return neighbors


# --------------------------------------------------------------------------- #
# Non-congruent hierarchy (RI3H.ec: parent0 / getParents / getChildren /
# centroidChild / isCentroidChild / centroidParent) — A2, geometric, exact vs pydggal.
# --------------------------------------------------------------------------- #
def _i3h_parent0(l9r, root, rix, sh):
    """parent0 (RI3H.ec:977-1003): the primary parent (analytic, one resolution
    up). Odd (sh>0) -> reset subHex to 0; even -> walk up one I9R level via
    iLRCFromLRtI + row/col mod-3, picking the parent's odd sub-hex (B/C/D). None
    at the res-0 root."""
    if l9r == 0 and sh == 0:
        return None
    if sh > 0:
        return pack_i3h(l9r, root, rix, 0)
    if root < 10:
        ir = _ilrc_from_lrti(l9r, root, rix)
        level, row, col = ir if ir is not None else (-1, -1, -1)
    else:
        level = l9r if (root <= 12 and rix == 0) else -1
        row, col = 0, 0
    p = 3 ** level if level >= 0 else 1
    r = rix // p
    c = rix % p
    rm3 = r % 3
    cm3 = c % 3
    sub = 1 if root > 9 else (2 if cm3 > 1 else (3 if rm3 > 1 else 1))
    res = _from_i9r(level - 1, row // 3, col // 3, sub, root if root > 9 else 0)
    return pack_i3h(*res) if res is not None else None


def _i3h_is_centroid_child(l9r, root, rix, sh):
    """isCentroidChild (RI3H.ec:2170-2196): does this cell have a single parent
    (odd B; polar even A; even A with row/col or row+col a multiple of 3)."""
    if sh > 0:
        return sh == 1
    if root == 10 or root == 11:
        return True
    p = 3 ** l9r
    r = rix // p
    c = rix % p
    return (r % 3 == 0 and c % 3 == 0) or ((r + c) % 3 == 0)


def _i3h_centroid_child(l9r, root, rix, sh):
    """centroidChild (RI3H.ec:2012-2042): the child sharing this cell's centroid."""
    if sh == 0:
        return pack_i3h(l9r, root, rix, 1)
    if root > 9:
        return pack_i3h(l9r + 1, root, 0, 0)
    p = 3 ** l9r
    row_op = (root + 1) >> 1
    col_op = root >> 1
    ix_op = rix // p
    row = row_op * p + ix_op
    col = (col_op - ix_op) * p + rix
    r = row * 3 + (2 if sh == 3 else (1 if sh == 2 else 0))
    c = col * 3 + (2 if sh == 2 else (1 if sh == 3 else 0))
    res = _from_i9r(l9r + 1, r, c, 0, 0)
    return pack_i3h(*res) if res is not None else None


def _nb(packed, direction):
    """A neighbour of an already-packed cell in one I3HNeighbor direction, or None."""
    if packed is None:
        return None
    return _i3h_get_neighbor(*unpack_i3h(packed), direction)


# --------------------------------------------------------------------------- #
# nullZone hierarchy: mirrored from pydggal verbatim
# --------------------------------------------------------------------------- #
# Policy (2026-07-29, extended to hierarchy 2026-08-21): mirror pydggal exactly
# for nullZone rather than inventing our own signalling. `centroid`/`text_id`
# were already mirrored; `parents`/`children` were not, and the generic path
# manufactured ordinary-looking zones from the "no such zone" sentinel --
# including a resolution-0 zone out of `_i3h_centroid_child`'s pack_i3h(31+1, ...).
#
# These are LOOKUPS, not arithmetic, and cannot be otherwise: NULL_ZONE is
# all-ones, and bits 62-63 lie outside every I3H bitfield (levelI9R:5:57,
# rootRhombus:4:53, rhombusIX:51:2, subHex:2:0), so any value rebuilt from an
# unpacked tuple loses them. DGGAL keeps them by manipulating the raw uint64.
#
# Captured from pydggal 0.0.6 (ISEA3H, DGGRS_getZoneParents / getZoneChildren on
# nullZone) and pinned against the LIVE oracle by tests/test_null_zone_hierarchy.py,
# so a dggal upgrade that changes them fails loudly instead of drifting.
# The duplicates and the self-reference in the children list are DGGAL's own --
# "mirror exactly" includes its degenerate output.
_NULL_ZONE_PARENTS = [
    18446744073709551612,
    4515077520210372780,
    4515077520210372776,
]
_NULL_ZONE_CHILDREN = [
    18446744073709551615,
    16478188292383308288,
    16496202690892790272,
    16478188292383308288,
    16496202690892790272,
    16550245886421236224,
]


def _i3h_get_parents_raw(value):
    """getParents (RI3H.ec:1247-1329): [parent0] if a centroid child, else
    [parent0, p1, p2] (entries MAY be None near pentagons/poles). The non-primary
    parents are neighbours of parent0 chosen by an odd-level right/topRight|
    bottomRight rule or an even-level centroid-delta dispatch."""
    l9r, root, rix, sh = unpack_i3h(value)
    p0 = _i3h_parent0(l9r, root, rix, sh)
    if _i3h_is_centroid_child(l9r, root, rix, sh):
        return [] if p0 is None else [p0]
    if sh > 0:  # Odd level
        p1 = _nb(p0, _RIGHT)
        p2 = _nb(p0, _TOP_RIGHT if sh == 2 else _BOTTOM_RIGHT)
    else:  # Even level -- centroid-delta dispatch
        cx, cy = _i3h_centroid_xy(l9r, root, rix, sh)
        p0cx, p0cy = _i3h_centroid_xy(*unpack_i3h(p0))
        dx = cx - p0cx
        dy = cy - p0cy
        p0cxi = math.floor(p0cx + 1e-11)
        p0cyi = math.floor(p0cy + 1e-11)
        on_bottom_left = (p0cy - p0cx + 1e-11 > 1) and (p0cx - p0cxi < 1e-11)
        on_top_right = (p0cx - p0cy + 1e-11 > 0) and (p0cy - p0cyi < 1e-11)
        if abs(dx) < 1e-11:
            if dy > 0:
                on_top_right_neg = (p0cx - p0cy - 1e-11 > 0) and (p0cy - p0cyi < 1e-11)
                p1 = _nb(p0, _BOTTOM_RIGHT)
                p2 = _nb(p0, _BOTTOM_LEFT if on_bottom_left else (_RIGHT if on_top_right_neg else _BOTTOM))
            else:
                p1 = _nb(p0, _TOP_LEFT)
                p2 = _nb(p0, _TOP)
        elif abs(dy) < 1e-11:
            if dx > 0:
                on_bottom_left_neg = (p0cy - p0cx - 1e-11 > 1) and (p0cx - p0cxi < 1e-11)
                p1 = _nb(p0, _TOP_RIGHT)
                p2 = _nb(p0, _TOP_LEFT if on_top_right else (_RIGHT if on_bottom_left_neg else _TOP))
            else:
                p1 = _nb(p0, _BOTTOM_LEFT)
                p2 = _nb(p0, _BOTTOM)
        else:
            if dx > 0:
                p1 = _nb(p0, _TOP_RIGHT)
                p2 = _nb(p0, _BOTTOM_RIGHT)
            else:
                p1 = _nb(p0, _TOP_LEFT)
                p2 = _nb(p0, _BOTTOM_LEFT)
    return [p0, p1, p2]


def _i3h_get_parents(value):
    """The cell's parents (1 or 3), nulls filtered, as packed ints."""
    if value == NULL_ZONE:
        return list(_NULL_ZONE_PARENTS)  # see the mirror table above
    return [p for p in _i3h_get_parents_raw(value) if p is not None]


def _i3h_centroid_parent(value):
    """centroidParent (RI3H.ec:1203-1224): parent0 if it is itself a centroid
    child, else the first centroid-child among the other parents; None if none."""
    l9r, root, rix, sh = unpack_i3h(value)
    cp = _i3h_parent0(l9r, root, rix, sh)
    if cp is not None and _i3h_is_centroid_child(*unpack_i3h(cp)):
        return cp
    for p in _i3h_get_parents_raw(value)[1:]:
        if p is not None and _i3h_is_centroid_child(*unpack_i3h(p)):
            return p
    return None


def _i3h_get_children(value):
    """getChildren (RI3H.ec:2044-2085): centroidChild + one child per boundary
    vertex (quantize the cell's own RAW planar vertices at the next resolution) —
    nPoints+1 total (6 pentagon / 7 hexagon). Poles use explicit fan formulas."""
    if value == NULL_ZONE:
        return list(_NULL_ZONE_CHILDREN)  # see the mirror table above
    l9r, root, rix, sh = unpack_i3h(value)
    out = [_i3h_centroid_child(l9r, root, rix, sh)]
    if root > 9:
        p = 3 ** l9r
        if sh == 0:
            if root == 10:
                out += [pack_i3h(l9r, (i - 1) * 2, p - 1, 2) for i in range(1, 6)]
            else:
                out += [pack_i3h(l9r, (i - 1) * 2 + 1, p * (p - 1), 3) for i in range(1, 6)]
        else:
            if root == 10:
                out += [pack_i3h(l9r + 1, (i - 1) * 2, 3 * p - 1, 0) for i in range(1, 6)]
            else:
                out += [pack_i3h(l9r + 1, (i - 1) * 2 + 1, 3 * p * (3 * p - 1), 0) for i in range(1, 6)]
    else:
        next_level = 2 * l9r + 1 + (1 if sh > 0 else 0)
        for vx, vy in _i3h_vertices(l9r, root, rix, sh):
            res = _i3h_from_centroid(next_level, vx, vy)
            out.append(pack_i3h(*res) if res is not None else None)
    return [c for c in out if c is not None]


# --------------------------------------------------------------------------- #
# Sub-zone closed-form count (RI3H.ec I3HZone::getSubZonesCount)
# --------------------------------------------------------------------------- #
def _i3h_count_sub_zones(nv: int, depth: int) -> int:
    """``I3HZone::getSubZonesCount`` (``RI3H.ec:2198-2202``): the closed-form
    count of sub-zones at relative ``depth`` for a cell with ``nv`` boundary
    points (6 = hexagon, 5 = pentagon)."""
    n_hex_sub_zones = 3 ** depth + 3 ** ((depth + 1) // 2) + 1 if depth > 0 else 1
    return (n_hex_sub_zones * nv + 5) // 6


def _i3h_sz_level(level_i9r: int, sub_hex: int, relative_depth: int) -> int:
    """The absolute I3H level (``level_i9r`` units, i.e. rows-of-9-cells) a
    sub-zone lives at, ``relative_depth`` below a parent cell described by
    ``(level_i9r, sub_hex)``: the parent's own absolute level is
    ``absolute_level(level_i9r, sub_hex)`` (``I3HIndexing.resolution``'s own
    formula -- shared via :func:`py4dggs.indexings.i3h.absolute_level` so the two
    never drift apart), plus ``relative_depth``. Shared by
    :func:`_i3h_sub_zone_centroids`, :func:`_i3h_sub_zones` and
    ``HexAperture3Topology.first_sub_zone`` -- final-review DRY extraction
    (Finding 3), pure refactor, no behavior change."""
    return absolute_level(level_i9r, sub_hex) + relative_depth


def _i3h_is_edge_hex(level_i9r: int, root: int, rix: int, sub_hex: int) -> bool:
    """Whether ``(level_i9r, root, rix, sub_hex)`` is an "edge hexagon" -- a
    hexagonal sub-hex sitting exactly on an interruption-adjacent rhombus
    boundary ``rix`` (``I3HSubZones.ec:1817``: ``rhombusIX && (subHex==0||
    subHex==1) && ...``). Requires ``sub_hex in (0, 1)`` -- a C/D sub-hex (2/3)
    at the same boundary ``rix`` is geometrically interior (see
    :func:`_i3h_centroid_xy`'s C/D offset), NOT edge-hex, even though its
    underlying grid cell touches the interruption (empirically confirmed
    against pydggal: forcing edge-hex on such a cell, e.g. ``B6-1-D``,
    produces a mismatch -- see :func:`_i3h_sub_zone_centroids`'s docstring).

    Some call sites historically also gated this on ``nv == 6`` (the cell's
    boundary-point count); that term is REDUNDANT and intentionally omitted
    here: ``rix != 0`` already excludes every pentagon-like cell reachable in
    this codebase (non-polar pentagons always have ``rix == 0``; the 2 polar
    pentagons, ``root >= 10``, are ALSO only ever produced with ``rix == 0`` --
    confirmed both algebraically and by a 20000-sample near-pole quantization
    sweep finding zero ``root >= 10 and rix != 0`` cells -- during this
    review's DRY extraction), so whenever ``sub_hex in (0, 1) and rix != 0``
    hold, ``nv`` is provably 6.

    Final-review DRY extraction (Finding 2): unifies the predicate previously
    duplicated at :func:`_i3h_first_sub_zone_centroid` (with the redundant
    ``nv == 6`` term), :func:`_i3h_sub_zone_centroids` (without it) and the
    ``tests/test_i3h_subzones.py`` test helper (a third copy). Pure refactor,
    no behavior change -- verified by an exhaustive sweep over reachable
    ``(level_i9r, root, rix, sub_hex)`` combinations before and after."""
    divs = 3 ** level_i9r
    south = bool(root & 1)
    return sub_hex in (0, 1) and rix != 0 and ((rix % divs == 0) if south else (rix // divs == 0))


# --------------------------------------------------------------------------- #
# I3H sub-zones (Task 1: interior hexagon; Task 2: + edge hexagon)
# --------------------------------------------------------------------------- #
def _i3h_first_sub_zone_centroid(level_i9r, root, rix, sub_hex, depth, edge_hex=None):
    """The starting centroid for sub-zone iteration (``getI3HFirstSubZoneCentroid``,
    ``I3HSubZones.ec:102-247``) -- a specific vertex of the parent cell's own
    polygon, picked by a tie-break search over :func:`_i3h_vertices`' output.

    Three of the four start-corner branches (odd-parent/odd-depth, odd-parent/
    even-depth, even-parent/odd-depth) run an unconditional vertex-search loop
    in the eC (guarded by ``#if 1``, with the index-shortcut ``#else`` arm
    actually dead code that still executes unconditionally alongside it,
    overriding only the ``root==10``/pentagon cases -- verified by reading the
    literal source, not a paraphrase) -- so these three need no edge-hex
    override: the search is already geometry-driven and edge-agnostic.

    The fourth (even-parent/even-depth, ``I3HSubZones.ec:200-238``) is
    index-based, not search-based (its loop is ``#if 0``-disabled): ``ix =
    (root==10 and subHex==0) ? 0 : (nv==6) ? 5 : 4``, then ``if(edgeHex &&
    southRhombus) ix = 4`` (``I3HSubZones.ec:231-235``) -- a south edge-hex
    parent starts from ``vertices[4]`` instead of ``vertices[5]``; north
    edge-hex is unaffected (still ``vertices[5]``, same as interior).

    Non-polar pentagon (``nv==5``, ``root<10``) overrides (Task 3) and polar
    (``root in (10,11)``, Task 4 -- the ``root==10``/``root==11`` sub-terms of
    each eC ternary), read directly off the literal eC source:
      - odd-parent/odd-depth (``:132-135``): ``tl = southRhombus ? 4 : 3``
        (the ``root==11`` sub-term of the eC's ternary never fires for
        ``root<10``).
      - odd-parent/even-depth (``:159-162``): ``tl = 3`` unconditionally
        (again the ``root==11`` sub-term is polar-only).
      - even-parent/odd-depth (``:184-188``): only fires when
        ``southRhombus`` (``left = 2``); north non-polar pentagon keeps the
        vertex-search loop's result unchanged (no override in the eC).
      - even-parent/even-depth (``:200-238``): both north and south non-polar
        pentagon land on ``vertices[4]`` -- north explicitly (the active
        ``if(nv==5 && (root<10||subHex!=0) && !southRhombus) { top=4; ... }``
        branch, ``:200-203``; the commented-out ``/* ... */`` block at
        ``:204-213`` is dead prose, not live code, per the brief), south via
        the SAME index formula the edge-hex/interior path already uses
        (``ix = (nv==6) ? 5 : 4`` defaults to 4 for ``nv==5`` since the
        south-edge-hex override at ``:234-235`` requires ``nv==6``) -- so
        this port needs no separate pentagon branch there, just deriving
        ``nv`` from ``len(verts)`` instead of assuming 6."""
    verts = _i3h_vertices(level_i9r, root, rix, sub_hex)
    nv = len(verts)
    odd_parent = sub_hex > 0
    odd_depth = depth & 1
    south = bool(root & 1)

    if odd_parent:
        if odd_depth:
            tl = 0
            for i in range(1, len(verts)):
                if verts[i][1] < verts[tl][1] or (
                    abs(verts[i][1] - verts[tl][1]) < 1e-11 and verts[i][0] > verts[tl][0]
                ):
                    tl = i
            # eC I3HSubZones.ec:132-135 -- the #else override arm (runs unconditionally
            # alongside the search, overriding only root==10/pentagon). Poles (Task 4):
            # North pole (root 10, subHex 1) -> tl=0; South pole (root 11, subHex 1) -> tl=1.
            if root == 10 and sub_hex == 1:
                tl = 0
            elif nv == 5:
                tl = (1 if (root == 11 and sub_hex == 1) else 4) if south else 3
            c = verts[tl]
        else:
            tl = 0
            for i in range(1, len(verts)):
                if verts[i][1] < verts[tl][1] or (
                    abs(verts[i][1] - verts[tl][1]) < 1e-11 and verts[i][0] < verts[tl][0]
                ):
                    tl = i
            # eC I3HSubZones.ec:159-162. Poles (Task 4): North pole (root 10, subHex 1)
            # -> tl=4; South pole (root 11, subHex 1) -> tl=0.
            if root == 10 and sub_hex == 1:
                tl = 4
            elif nv == 5:
                tl = 0 if (root == 11 and sub_hex == 1) else 3
            c = verts[tl]
    else:
        if odd_depth:
            left = 0
            for i in range(1, len(verts)):
                if verts[i][0] < verts[left][0]:
                    left = i
            # eC I3HSubZones.ec:184-188. Poles (Task 4): North pole (root 10, subHex 0)
            # -> left=0; South pole (root 11, subHex 0) -> left=4.
            if root == 10 and sub_hex == 0:
                left = 0
            elif nv == 5 and south:
                left = 4 if (root == 11 and sub_hex == 0) else 2
            c = verts[left]
        else:
            # (root==10 and sub_hex==0) -> 0; else nv==6 -> 5 (south edge-hex ->
            # 4); nv==5 (non-polar pentagon, either N or S) -> 4 as well, via
            # the same default (see docstring).
            if edge_hex is None:  # caller may already have it (avoids a redundant call)
                edge_hex = _i3h_is_edge_hex(level_i9r, root, rix, sub_hex)
            if root == 10 and sub_hex == 0:
                c = verts[0]
            else:
                ix = 5 if nv == 6 else 4
                if edge_hex and south:
                    ix = 4
                c = verts[ix]

    cx, cy = c
    if cx > 5 and cy > 5:
        cx -= 5
        cy -= 5
    elif cx < 0 and cy < 1:
        cx += 5
        cy += 5
    return cx, cy


def _gen_odd_parent_odd_depth(first_cx, first_cy, depth, u, south, edge_hex, is_pentagon=False, polar_pentagon=False):
    """``generateOddParentOddDepth`` (``I3HSubZones.ec:250-432``), interior +
    edge-hex + non-polar pentagon + polar pentagon: row-start offsets
    ``I3HSubZones.ec:292-293`` (first half) / ``373-374`` (second half); in-row
    scanline ``I3HSubZones.ec:337-351`` (first half) / ``383-421`` (second half).

    Polar pentagon (Task 4, ``polar_pentagon = root >= 10``; both poles are
    ``nv==5`` so ``is_pentagon`` is also ``True``, and exactly one of
    ``north_pentagon``/``south_pentagon`` holds). The eC tests ``if(polarPentagon)``
    FIRST in every scanline dispatch, so the polar branch is prepended ahead of
    the non-polar ones. First half (``:290-336``): diagonal row-start
    ``move5x6Vertex(rc, first, ±r*u, ±r*u)``; the in-row is a two-stage
    ``move5x6Vertex3`` + ``cross5x6Interruption`` chain with a ``maxB``/``c``
    overflow into a SECOND cross (``:317-332``). Second half (``:362-421``):
    ``nCols--`` per row (like north pentagon, ``:365``), and the merged
    ``polarPentagon || edgeHex || northPentagon`` in-row branch (``:383-401``)
    where polar uses ``move5x6Vertex``+cross then a final offset with dx forced
    to 0 and dy = ``±b*u`` by hemisphere.

    North non-polar pentagon (``northPentagon = nv==5 and not south``) reuses
    almost all of the north edge-hex code path verbatim -- the eC merges the
    two conditions (``edgeHex || northPentagon``) everywhere except one spot:
    the second-half in-row ``b``-truthy (crossing) sub-case, where edge-hex
    explicitly crosses via ``move5x6Vertex`` + ``cross5x6Interruption`` but
    pure north-pentagon (not edge-hex) crosses implicitly via
    ``move5x6Vertex3`` directly (``:389-397``) -- ported as the ``edge_hex``
    if/else below. North pentagon ALSO shrinks the second-half row's column
    count by one extra per row (``nCols--`` at ``:365-366``, compounding via
    the loop's own per-row decrement -- ported as the extra ``n_cols -= 1``
    inside the loop body, not just the trailing one).

    South non-polar pentagon (``southPentagon``) does NOT share code with
    south edge-hex here at all -- the first half is untouched (plain interior
    default, matching the eC which only tests ``northPentagon``/``edgeHex``,
    never ``southPentagon``, in the first half); the second half only adds a
    column-skip mechanism (``I3HSubZones.ec:362-421``): ``n = (nCols-r)/2``
    computed at each row, then ``if(col==n) col = nCols-1-n`` (skip an
    interruption-straddling run of columns) -- ported as the ``while``-loop
    ``col`` mutation below (needed because Task 1/2 never required this
    mechanism, only introduced for pentagons)."""
    north_pentagon = is_pentagon and not south
    south_pentagon = is_pentagon and south
    sgn = -1 if south_pentagon else 1
    out = []
    n_half_rows = 3 ** ((depth - 1) // 2)
    max_cols = 2 * n_half_rows + 1
    n_cols = max_cols - (max_cols - 1) // 2

    r = 0
    while r <= n_half_rows:
        if polar_pentagon:
            rc = move5x6_vertex_v1(first_cx, first_cy, sgn * r * u, sgn * r * u)
        elif north_pentagon or (edge_hex and not south):
            rc = move5x6_vertex_v1(first_cx, first_cy, 0, r * u)
        else:
            rc = move5x6_vertex_v1(first_cx, first_cy, -r * u, 0)
        for col in range(n_cols):
            if polar_pentagon:
                a = min(col, r); b = col - a
                t = move5x6_vertex_v3(rc[0], rc[1], 0, (a if south_pentagon else -a) * u)
                if b:
                    max_b = n_cols - 1 - 2 * r; c = 0
                    if b > max_b:
                        c = b - max_b; b = max_b
                    left = north_pentagon or r >= n_half_rows
                    i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, left)
                    mover = move5x6_vertex_v1 if c else move5x6_vertex_v3
                    cen = mover(i2[0], i2[1], (b if south_pentagon else -b) * u, (b if south_pentagon else -b) * u)
                    if c:
                        i2 = _cross5x6_interruption(cen[0], cen[1], south_pentagon, left)
                        cen = move5x6_vertex_v3(
                            i2[0], i2[1],
                            (0 if r >= n_half_rows else c * u) if south_pentagon else -c * u,
                            -c * u if (south_pentagon and r >= n_half_rows) else 0,
                        )
                    out.append(cen)
                else:
                    out.append(t)
            elif edge_hex or north_pentagon:
                a = min(col, n_half_rows)
                b = col - a
                mover = move5x6_vertex_v1 if b else move5x6_vertex_v3
                t = mover(rc[0], rc[1], a * u, a * u if south else 0)
                if b:
                    i2 = _cross5x6_interruption(t[0], t[1], south, False)
                    out.append(move5x6_vertex_v3(i2[0], i2[1], b * u, 0 if south else b * u))
                else:
                    out.append(t)
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, col * u))
        r += 1
        n_cols += 1

    r, n_cols = 1, max_cols - 1
    while r <= n_half_rows:
        n = -1
        if north_pentagon or polar_pentagon:
            n_cols -= 1
        elif south_pentagon:
            n = (n_cols - r) // 2

        if polar_pentagon:
            v = (-1 - r * u) if south_pentagon else (1 + r * u)
            rc = move5x6_vertex_v1(first_cx, first_cy, v, v)
        elif north_pentagon or (edge_hex and not south):
            rc = move5x6_vertex_v1(first_cx, first_cy, r * u, (r + n_half_rows) * u)
        else:
            rc = move5x6_vertex_v1(first_cx, first_cy, -n_half_rows * u, r * u)

        col = 0
        while col < n_cols:
            if polar_pentagon or edge_hex or north_pentagon:
                a = min(col, n_half_rows - r)
                b = col - a
                dx_a = (-a if (polar_pentagon and south_pentagon) else a) * u
                dy_a = a * u if (south and not polar_pentagon) else 0
                if b:
                    if polar_pentagon or edge_hex:
                        t = move5x6_vertex_v1(rc[0], rc[1], dx_a, dy_a)
                        t = _cross5x6_interruption(t[0], t[1], south, south_pentagon and polar_pentagon)
                    else:  # north_pentagon only (not edge_hex): direct v3 crossing
                        t = move5x6_vertex_v3(rc[0], rc[1], dx_a, dy_a)
                    dx_b = 0 if polar_pentagon else b * u
                    dy_b = (-b * u if polar_pentagon else 0) if south else b * u
                    out.append(move5x6_vertex_v1(t[0], t[1], dx_b, dy_b))
                else:
                    out.append(move5x6_vertex_v3(rc[0], rc[1], dx_a, dy_a))
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, col * u))

            if col == n:
                col = n_cols - 1 - n  # Skip interruption (south pentagon; n=-1 for poles)
            col += 1
        r += 1
        n_cols -= 1
    return out


def _gen_even_parent_odd_depth(first_cx, first_cy, depth, u, south, edge_hex, is_pentagon=False, polar_pentagon=False):
    """``generateEvenParentOddDepth`` (``I3HSubZones.ec:435-821``).

    Polar pentagon (Task 4, prepended ``if(polarPentagon)`` branches, both
    halves). First half row-start (``:488-501``) has a ``r > nHalfRows/2``
    ``cross5x6Interruption`` reflection; the in-row (``:539-604``) is the file's
    most elaborate chain: a ``crossingLeft`` parity flag, an ``oddRow``
    correction that decrements ``b`` and nudges ``t``/``i2`` by ``sgn*u/3`` on
    even rows (``:576-584``), and a ``b > n`` overflow into a second cross. The
    second-half row-start (``:664-678``) double-crosses for ``r > nHalfRows/2``;
    its in-row (``:720-770``) is a three-stage ``b > n`` cascade. All offsets
    carry ``sgn = southPentagon ? -1 : 1`` and the loop bound already shrinks by
    ``r`` for ``nv==5`` (the shared ``n_cols_eff`` below). South
    edge-hex needs a distinct row-start/scanline for ``r > nHalfRows/2`` in
    BOTH halves (``I3HSubZones.ec:503-526`` row-start / ``607-627`` in-row,
    first half; second-half row-start/in-row reuse the same shape at
    ``680-689``/``774-794``, ported inline below since the eC recomputes the
    same ``i0``/``i1`` reflection each half rather than sharing it). North
    edge-hex only diverges in the SECOND half, and only for ``r >
    nHalfRows/2`` (``I3HSubZones.ec:690-703`` row-start / ``772-773`` in-row,
    the ``crosses`` flag); for ``r <= nHalfRows/2`` and the entire first half,
    north edge-hex is arithmetically identical to the interior/basic path.

    Non-polar pentagon (Task 3, ``polarPentagon`` hardcoded ``False``):
      - South pentagon (``southPentagon = nv==5 and south``) shares the SAME
        row-start as south edge-hex for ``r > nHalfRows/2`` (``:503,
        :680``: the eC tests ``edgeHex || southPentagon`` throughout this
        generator). The in-row computation also mostly matches, except right
        at the boundary row ``r == nHalfRows`` in the first half (``:616``:
        ``!southPentagon || r < nHalfRows``), where south pentagon skips the
        explicit ``cross5x6Interruption`` call and crosses directly via
        ``move5x6Vertex3`` instead (ported as the ``r < n_half_rows`` guard
        below); and the second half's ``n`` (crossing threshold) differs
        entirely for south pentagon (``nHalfRows - r``, ``:776``) vs. south
        edge-hex (``nCols - (nHalfRows-2r) - 1``), with south pentagon also
        skipping the explicit cross (direct ``move5x6Vertex3``, ``:783-784``).
      - North pentagon (``northPentagon = nv==5 and not south``) has NO
        row-start override anywhere in this generator (only ``southRhombus &&
        (edgeHex||southPentagon)`` and the north-edge-hex-only ``crosses``
        branch test row-start; north pentagon falls to the plain interior
        default both halves). The in-row computation gets its own dedicated
        branch only in the first half, for ``r > nHalfRows/2``
        (``:628-643``); the second half's north-pentagon in-row branch
        (``:795-806``) is unconditional on ``r``.
      - ``nCols`` is reduced by ``r`` in the second half for EITHER pentagon
        orientation (``int skip = nv==5 ? r : 0;``, ``:654``, ``:712``) --
        ported as the ``n_cols_eff`` bound below."""
    north_pentagon = is_pentagon and not south
    south_pentagon = is_pentagon and south
    sgn = -1 if south_pentagon else 1
    out = []
    n_half_rows = 3 ** ((depth - 1) // 2)
    max_cols = 2 * n_half_rows + 1
    n_cols = max_cols - (max_cols - 1) // 2

    # South edge-hex/south-pentagon reflection anchor: depends only on
    # first_cx/first_cy, which are fixed for the whole call, so it's
    # loop-invariant across both scanline loops below (each used to
    # recompute it independently, once per qualifying row).
    if south and (edge_hex or south_pentagon):
        ix = math.trunc(first_cx - 1e-11)
        iy = math.trunc(first_cy - 1e-11)
        i0x = first_cx + ((iy + 1) - first_cy)
        i0y = iy + 1
        i1x = i0y - 1
        i1y = iy + 1 + (ix + 1 - i0x)

    r = 0
    while r <= n_half_rows:
        if polar_pentagon:
            if r > n_half_rows // 2:
                a = n_half_rows // 2; b = r - (a + 1)
                rc = move5x6_vertex_v1(first_cx, first_cy, (a * 2 + 1) * sgn * u / 3, (a * 1 + 1) * sgn * u / 3)
                i1 = _cross5x6_interruption(rc[0], rc[1], south_pentagon, south_pentagon)
                rc = move5x6_vertex_v3(i1[0], i1[1], (1 + b * 1) * sgn * u / 3, (1 + b * 2) * sgn * u / 3)
            else:
                rc = move5x6_vertex_v1(first_cx, first_cy, r * (-2 if south_pentagon else 2) * u / 3, r * sgn * u / 3)
        elif south and (edge_hex or south_pentagon) and r > n_half_rows // 2:
            b = r - n_half_rows // 2 - 1
            rc = move5x6_vertex_v3(i1x, i1y, (b * 2 + 1) * u / 3, ((n_half_rows + 1) // 2 + b) * u / 3)
        else:
            rc = move5x6_vertex_v1(first_cx, first_cy, r * u / 3, r * 2 * u / 3)

        for col in range(n_cols):
            if polar_pentagon:
                col_rem = col
                crossing_left = (r >= n_half_rows) if south_pentagon else (r < n_half_rows)
                if r > n_half_rows // 2:
                    n = (r - n_half_rows // 2) * 2 - 1
                    a = min(col, n)
                    col_rem = col - a
                    mover = move5x6_vertex_v1 if col_rem else move5x6_vertex_v3
                    start = mover(rc[0], rc[1], a * sgn * u / 3, -1 * a * sgn * u / 3)
                    if col_rem:
                        n = (n_half_rows - r) // 2 + (n_half_rows - r)
                        start = _cross5x6_interruption(start[0], start[1], south_pentagon, crossing_left)
                else:
                    start = rc
                    n = (n_half_rows + r) // 2
                if col_rem:
                    a = min(col_rem, n); b = col_rem - a
                    mover = move5x6_vertex_v1 if b else move5x6_vertex_v3
                    t = mover(start[0], start[1], -a * sgn * u / 3, -2 * a * sgn * u / 3)
                    if b:
                        odd_row = r & 1
                        if not odd_row:
                            b -= 1
                            t = [t[0] - sgn * u / 3, t[1] - sgn * u / 3]
                        i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, crossing_left)
                        if not odd_row:
                            i2 = [i2[0] - sgn * u / 3, i2[1] - sgn * u / 3]
                        if b > n:
                            a = min(b, n); b -= a
                            mover2 = move5x6_vertex_v1 if a else move5x6_vertex_v3
                            t = mover2(i2[0], i2[1], a * sgn * -2 * u / 3, -a * sgn * u / 3)
                            if a:
                                i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, crossing_left)
                            else:
                                i2 = t
                            out.append(move5x6_vertex_v3(i2[0], i2[1], -b * sgn * u / 3, b * sgn * u / 3))
                        else:
                            out.append(move5x6_vertex_v3(i2[0], i2[1], b * -2 * sgn * u / 3, -b * sgn * u / 3))
                    else:
                        out.append(t)
                else:
                    out.append(start)
            elif south and (edge_hex or south_pentagon) and r > n_half_rows // 2:
                n = (r - n_half_rows // 2) * 2 - 1
                a = min(col, n)
                b = col - a
                if b:
                    if not south_pentagon or r < n_half_rows:
                        t = move5x6_vertex_v1(rc[0], rc[1], -a * u / 3, -2 * a * u / 3)
                        i2 = _cross5x6_interruption(t[0], t[1], True, True)
                    else:
                        i2 = move5x6_vertex_v3(rc[0], rc[1], -a * u / 3, -2 * a * u / 3)
                    out.append(move5x6_vertex_v3(i2[0], i2[1], b * u / 3, -b * u / 3))
                else:
                    out.append(move5x6_vertex_v3(rc[0], rc[1], -a * u / 3, -2 * a * u / 3))
            elif north_pentagon and r > n_half_rows // 2:
                n = 2 * n_half_rows - r
                a = min(col, n)
                b = col - a
                mover = move5x6_vertex_v1 if b else move5x6_vertex_v3
                t = mover(rc[0], rc[1], a * u / 3, -a * u / 3)
                if b:
                    i2 = _cross5x6_interruption(t[0], t[1], False, False)
                    out.append(move5x6_vertex_v3(i2[0], i2[1], b * 2 * u / 3, b * u / 3))
                else:
                    out.append(t)
            else:
                out.append(move5x6_vertex_v3(rc[0], rc[1], col * u / 3, -col * u / 3))
        r += 1
        n_cols += 1

    r, n_cols = 1, max_cols - 1
    while r <= n_half_rows:
        crosses = False
        if polar_pentagon:
            a = n_half_rows // 2; b = n_half_rows - (a + 1); aa = min(r, n_half_rows // 2)
            rc = move5x6_vertex_v1(first_cx, first_cy, (a * 2 + 1) * sgn * u / 3, (a + 1) * sgn * u / 3)
            i1 = _cross5x6_interruption(rc[0], rc[1], south_pentagon, south_pentagon)
            rc = move5x6_vertex_v3(i1[0], i1[1], (1 + b + aa * 2) * sgn * u / 3, (1 + b * 2 + aa) * sgn * u / 3)
            if r > n_half_rows // 2:
                bb = r - (aa + 1)
                rc = [rc[0] + sgn * u / 3, rc[1] + sgn * u / 3]
                i1 = _cross5x6_interruption(rc[0], rc[1], south_pentagon, south_pentagon)
                rc = move5x6_vertex_v3(i1[0], i1[1], (1 + bb) * sgn * u / 3, (1 + bb * 2) * sgn * u / 3)
        elif south and (edge_hex or south_pentagon):
            b = n_half_rows - n_half_rows // 2 - 1
            rc = move5x6_vertex_v3(
                i1x, i1y,
                u / 3 + b * 2 * u / 3 + r * u / 3,
                (n_half_rows + 1) / 2.0 * u / 3 + b * u / 3 - r * u / 3,
            )
        elif not south and edge_hex and r > n_half_rows // 2:
            i_ = r - n_half_rows // 2 - 1
            iix = first_cx + n_half_rows * 2 * u / 3
            iiy = first_cy + n_half_rows * 2 * u / 3
            iy = math.trunc(iiy - 1e-11)
            tx = iy + 2 - (iiy - iy)
            ty = iix
            rc = move5x6_vertex_v1(tx, ty, -(n_half_rows - 1) / 2.0 * u / 3 + i_ * u / 3, u / 3 + i_ * 2 * u / 3)
            crosses = True
        else:
            dx = n_half_rows * u / 3 + r * 2 * u / 3
            dy = n_half_rows * 2 * u / 3 + r * u / 3
            rc = move5x6_vertex_v1(first_cx, first_cy, dx, dy)

        n_cols_eff = n_cols - r if is_pentagon else n_cols
        for col in range(n_cols_eff):
            if polar_pentagon:
                b = col
                n = n_half_rows - 2 * r
                if r <= n_half_rows // 2:
                    a = min(col, n); b = col - a
                    mover = move5x6_vertex_v1 if b else move5x6_vertex_v3
                    t = mover(rc[0], rc[1], sgn * a * u / 3, sgn * -a * u / 3)
                    if b:
                        i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, south_pentagon)
                        t = i2
                        n = r
                else:
                    n = n_half_rows - r
                    t = rc
                if b:
                    if b > n:
                        a = min(b, n); b -= a
                        i2 = move5x6_vertex_v3(t[0], t[1], sgn * a * 2 * u / 3, sgn * a * u / 3)
                        if b > n:
                            a = min(b, n); b -= a
                            t = move5x6_vertex_v1(i2[0], i2[1], sgn * a * u / 3, sgn * a * 2 * u / 3)
                            i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, south_pentagon)
                            out.append(move5x6_vertex_v3(i2[0], i2[1], sgn * -b * u / 3, sgn * b * u / 3))
                        else:
                            out.append(move5x6_vertex_v3(i2[0], i2[1], sgn * b * u / 3, sgn * b * 2 * u / 3))
                    else:
                        out.append(move5x6_vertex_v3(t[0], t[1], sgn * b * 2 * u / 3, sgn * b * u / 3))
                else:
                    out.append(t)
            elif crosses:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * 2 * u / 3, col * u / 3))
            elif south and (edge_hex or south_pentagon):
                n = (n_half_rows - r) if south_pentagon else (n_cols - (n_half_rows - 2 * r) - 1)
                a = min(col, n)
                b = max(0, col - n)
                if b:
                    if south_pentagon:
                        i2 = move5x6_vertex_v3(rc[0], rc[1], -a * u / 3, -2 * a * u / 3)
                    else:
                        t = move5x6_vertex_v1(rc[0], rc[1], -a * u / 3, -2 * a * u / 3)
                        i2 = _cross5x6_interruption(t[0], t[1], True, True)
                    out.append(move5x6_vertex_v3(i2[0], i2[1], b * u / 3, -b * u / 3))
                else:
                    out.append(move5x6_vertex_v3(rc[0], rc[1], -a * u / 3, -2 * a * u / 3))
            elif north_pentagon:
                n = n_half_rows - r
                a = min(col, n)
                b = col - a
                t = move5x6_vertex_v3(rc[0], rc[1], a * u / 3, -a * u / 3)
                if b:
                    out.append(move5x6_vertex_v3(t[0], t[1], b * 2 * u / 3, b * u / 3))
                else:
                    out.append(t)
            else:
                out.append(move5x6_vertex_v3(rc[0], rc[1], col * u / 3, -col * u / 3))
        r += 1
        n_cols -= 1
    return out


def _gen_even_parent_even_depth(first_cx, first_cy, depth, u, south, edge_hex, is_pentagon=False, polar_pentagon=False):
    """Polar pentagon (Task 4): this generator is STRUCTURALLY ASYMMETRIC in the
    eC and is replicated as-is, NOT unified. In the row-start the SOUTH pole is a
    top-level ``if(polarPentagon && southRhombus)`` (``:862``) while the NORTH
    pole is NESTED inside ``(edgeHex||nv==5) && !southRhombus -> if(polarPentagon)``
    (``:874-887``); in the in-row BOTH poles are nested inside the north/south
    ``(edgeHex||nv==5)`` branches (``:903-985``, ``:1092-1186``, ``:1259-1333``).
    Each polar sub-case runs its own ``nb``/``c`` overflow-into-second-cross
    chain, and the main-portion south pole even reroutes into the
    ``colSkip && southPentagon`` branch (``:1058-1090``). The second-cap
    row-start uses fixed reflection points ``{2+nCapRows*u, 4-nCapRows*u}`` /
    ``{3-nCapRows*u, 2+nCapRows*u}`` (``:1206,:1211``). Do not refactor.

    ``generateEvenParentEvenDepth`` (``I3HSubZones.ec:824-1363``), the
    cap/main/cap stateful-``nCols`` scanline. Edge-hex branches: first cap
    row-start/in-row ``I3HSubZones.ec:874-890``/``903-985``; main row-start/
    in-row ``1021-1047``/``1092-1188``; second cap row-start/in-row
    ``1197-1256``/``1263-1335`` -- each filtered to the ``edgeHex`` (not
    ``nv==5``) sub-case, with north (``!southRhombus``) and south branches
    kept separate exactly as the eC does (they are not mirror images: north
    uses ``move5x6Vertex``/south uses ``move5x6Vertex3`` as the b==0 "no
    crossing yet" mover in the first cap, for instance).

    Non-polar pentagon (Task 3, ``polarPentagon`` hardcoded ``False``; this
    generator's own local booleans are ``southPentagon = nv==5 and
    southRhombus`` -- there is no separate ``northPentagon`` local, the eC
    tests ``nv==5 && !southRhombus`` directly):
      - First cap: north/south pentagon share the SAME row-start/in-row
        branches as north/south edge-hex (the eC's condition is literally
        ``edgeHex || nv==5``) -- only one south-specific line differs
        (``:978``: ``move5x6Vertex3(centroid, nv==5 ? t : i2, b*u, nv==5 ?
        b*u : 0)`` -- pentagon uses the PRE-cross ``t`` with a different
        offset than edge-hex's post-cross ``i2``).
      - Main portion: ``colSkip = (nv==5 && r>nMidRows/2) ? r-nMidRows/2 : 0``
        shrinks the per-row column count for BOTH orientations once past the
        midpoint (``:996``). For south pentagon past the midpoint, this
        additionally reroutes into an entirely separate ``if(colSkip &&
        southPentagon)`` in-row branch (``:1058-1090``) that bypasses the
        edge-hex-shared logic altogether. North pentagon past the midpoint
        stays in the shared ``(edgeHex||nv==5) && !southRhombus`` branch, but
        forces the ``v3``-then-plain-offset path unconditionally
        (``!b || (nv==5 && r>nMidRows/2)``, ``:1099``) instead of the
        b-dependent mover choice edge-hex uses, and skips the explicit
        cross (``:1116-1117``, direct offset from ``t``). South pentagon at
        or before the midpoint reuses the south edge-hex branch's shape but
        skips the cross too (``:1151-1155``, direct offset from the
        PRE-cross ``t``, not ``i2``).
      - Second cap: row-loop bound shrinks by ``endCapSkip = nv==5 ?
        (nCapRows+1)/2 : 0`` (``:846``, ``:1195``) for EITHER orientation.
        In-row: pentagon reuses edge-hex's mover choices but skips the
        explicit cross when ``a`` and ``b`` are both truthy (``i2`` computed
        directly via ``move5x6Vertex3`` at ``:1283``/``move5x6Vertex`` at
        ``:1309`` instead of cross-then-v3), rescales ``b`` by a shift
        term (``:1276``, ``:1319``), and an unconditional ``if(nv==5 &&
        col==n) col = nCols-n-1`` column-skip (``:1351``, mirrors the
        odd-parent/odd-depth generator's south-pentagon skip mechanism)."""
    north_pentagon = is_pentagon and not south
    south_pentagon = is_pentagon and south
    out = []
    n_cap_rows = 3 ** ((depth - 2) // 2)
    n_mid_rows = 2 * n_cap_rows + 1
    end_cap_skip = (n_cap_rows + 1) // 2 if is_pentagon else 0
    n_cols = 1

    # First cap
    r = 0
    while r < n_cap_rows:
        # row-start (eC :862-892) -- south pole is top-level; north pole nested
        if polar_pentagon and south:
            if r > n_cap_rows // 2:
                t = move5x6_vertex_v1(first_cx, first_cy, -n_cap_rows * u, -(n_cap_rows - r) * u)
                i2 = _cross5x6_interruption(t[0], t[1], True, True)
                rc = move5x6_vertex_v1(i2[0], i2[1], 0, -(2 * (r - n_cap_rows // 2) - 1) * u)
            else:
                rc = move5x6_vertex_v3(first_cx, first_cy, -2 * r * u, -r * u)
        elif (edge_hex or north_pentagon) and not south:
            if polar_pentagon:
                if r > n_cap_rows // 2:
                    t = move5x6_vertex_v1(first_cx, first_cy, n_cap_rows * u, (n_cap_rows - r) * u)
                    i2 = _cross5x6_interruption(t[0], t[1], False, False)
                    rc = move5x6_vertex_v1(i2[0], i2[1], 0, (2 * (r - n_cap_rows // 2) - 1) * u)
                else:
                    rc = move5x6_vertex_v1(first_cx, first_cy, 2 * r * u, r * u)
            else:
                rc = move5x6_vertex_v1(first_cx, first_cy, -r * u, r * u)
        else:
            rc = move5x6_vertex_v1(first_cx, first_cy, r * -2 * u, r * -1 * u)

        for col in range(n_cols):
            if (edge_hex or north_pentagon) and not south:
                if r > n_cap_rows // 2:
                    n = (2 * r - n_cap_rows) if polar_pentagon else (n_cap_rows + r)
                    a = min(col, n)
                    b = col - a
                    mover = move5x6_vertex_v1 if b else move5x6_vertex_v3
                    t = mover(rc[0], rc[1], 0 if polar_pentagon else a * u, -a * u if polar_pentagon else 0)
                    if b:
                        i2 = _cross5x6_interruption(t[0], t[1], False, polar_pentagon)
                        if polar_pentagon:
                            nb = 2 * n_cap_rows - r
                            b = min(b, nb); c = col - (a + b)
                            mover2 = move5x6_vertex_v1 if c else move5x6_vertex_v3
                            t2 = mover2(i2[0], i2[1], -b * u, -b * u)
                            if c:
                                i2 = _cross5x6_interruption(t2[0], t2[1], False, True)
                                out.append(move5x6_vertex_v3(i2[0], i2[1], -c * u, 0))
                            else:
                                out.append(t2)
                        else:
                            out.append(move5x6_vertex_v3(i2[0], i2[1], b * u, b * u))
                    else:
                        out.append(t)
                elif polar_pentagon:
                    out.append(move5x6_vertex_v1(rc[0], rc[1], -col * u, -col * u))
                else:
                    out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, 0))
            elif (edge_hex or south_pentagon) and south:
                if r > n_cap_rows // 2:
                    n = (2 * r - n_cap_rows) if polar_pentagon else (n_cap_rows + r)
                    a = min(col, n)
                    b = col - a
                    mover = move5x6_vertex_v1 if b else move5x6_vertex_v3
                    t = mover(rc[0], rc[1], 0 if polar_pentagon else a * u, a * u)
                    if b:
                        i2 = _cross5x6_interruption(t[0], t[1], True, False)
                        if polar_pentagon:
                            nb = 2 * n_cap_rows - r
                            b = min(b, nb); c = col - (a + b)
                            mover2 = move5x6_vertex_v1 if c else move5x6_vertex_v3
                            t2 = mover2(i2[0], i2[1], b * u, b * u)
                            if c:
                                i2 = _cross5x6_interruption(t2[0], t2[1], True, False)
                                out.append(move5x6_vertex_v3(i2[0], i2[1], c * u, 0))
                            else:
                                out.append(t2)
                        elif south_pentagon:
                            out.append(move5x6_vertex_v3(t[0], t[1], b * u, b * u))
                        else:
                            out.append(move5x6_vertex_v3(i2[0], i2[1], b * u, 0))
                    else:
                        out.append(t)
                else:
                    out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, col * u))
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, col * u))
        r += 1
        n_cols += 3

    # Main portion
    r = 0
    while r < n_mid_rows:
        col_skip = (r - n_mid_rows // 2) if (is_pentagon and r > n_mid_rows // 2) else 0

        # row-start (eC :1005-1047) -- south pole top-level; north pole nested
        if polar_pentagon and south:
            a = min(r, n_cap_rows); b = r - a
            t = move5x6_vertex_v1(first_cx, first_cy, -n_cap_rows * u, 0)
            i2 = _cross5x6_interruption(t[0], t[1], True, True)
            t2 = move5x6_vertex_v1(i2[0], i2[1], -a * u, -(a // 2 + n_cap_rows) * u)
            if b:
                i2 = _cross5x6_interruption(t2[0], t2[1], True, True)
                rc = move5x6_vertex_v1(i2[0], i2[1], -(b // 2) * u, -b * u)
            else:
                rc = t2
        elif (edge_hex or north_pentagon) and not south:
            if polar_pentagon:
                a = min(r, n_cap_rows); b = r - a
                t = move5x6_vertex_v1(first_cx, first_cy, n_cap_rows * u, 0)
                i2 = _cross5x6_interruption(t[0], t[1], False, False)
                t2 = move5x6_vertex_v1(i2[0], i2[1], a * u, (a // 2 + n_cap_rows) * u)
                if b:
                    i2 = _cross5x6_interruption(t2[0], t2[1], False, False)
                    rc = move5x6_vertex_v1(i2[0], i2[1], (b // 2) * u, b * u)
                else:
                    rc = t2
            else:
                rc = move5x6_vertex_v1(
                    first_cx, first_cy, n_cap_rows * -u + ((r + 1) >> 1) * u, n_cap_rows * u + r * u
                )
        else:
            rc = move5x6_vertex_v1(
                first_cx, first_cy,
                n_cap_rows * -2 * u + (r >> 1) * -u,
                n_cap_rows * -1 * u + (r >> 1) * u + (r & 1) * u,
            )

        for col in range(n_cols - col_skip):
            if col_skip and south_pentagon:
                jump_col = (n_cols - col_skip) >> 1
                if col <= jump_col:
                    if polar_pentagon:
                        out.append(move5x6_vertex_v1(rc[0], rc[1], -col * u, 0))
                    else:
                        out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, col * u))
                else:
                    if polar_pentagon:
                        t = move5x6_vertex_v1(rc[0], rc[1], -jump_col * u, 0)
                        i2 = _cross5x6_interruption(t[0], t[1], True, True)
                        out.append(move5x6_vertex_v1(i2[0], i2[1], 0, -(col - jump_col) * u))
                    else:
                        start_shift = col_skip + jump_col + 1
                        extra_cols = col - jump_col - 1
                        t = move5x6_vertex_v1(rc[0], rc[1], start_shift * u, start_shift * u)
                        out.append(move5x6_vertex_v1(t[0], t[1], extra_cols * u, extra_cols * u))
            elif (edge_hex or north_pentagon) and not south:
                n = (n_cap_rows + r // 2) if polar_pentagon else (n_cols - 1 - n_cap_rows - r // 2)
                a = min(col, n)
                b = col - a
                pentagon_tail = is_pentagon and r > n_mid_rows // 2
                if not b or pentagon_tail:
                    if polar_pentagon and r <= n_cap_rows:
                        t = move5x6_vertex_v3(rc[0], rc[1], 0, -a * u)
                    else:
                        t = move5x6_vertex_v3(rc[0], rc[1], a * u, 0)
                else:
                    if polar_pentagon and r <= n_cap_rows:
                        t = move5x6_vertex_v1(rc[0], rc[1], 0, -a * u)
                    else:
                        t = move5x6_vertex_v1(rc[0], rc[1], a * u, 0)
                if b:
                    if pentagon_tail:
                        out.append(move5x6_vertex_v1(t[0], t[1], 0 if polar_pentagon else b * u, b * u))
                    else:
                        i2 = _cross5x6_interruption(t[0], t[1], False, polar_pentagon)
                        if polar_pentagon:
                            nb = n_cap_rows - r
                            b = min(b, nb); c = col - (a + b)
                            t2 = move5x6_vertex_v1(i2[0], i2[1], -b * u, -b * u)
                            if c:
                                i2 = _cross5x6_interruption(t2[0], t2[1], False, True)
                                out.append(move5x6_vertex_v3(i2[0], i2[1], -c * u, 0))
                            else:
                                out.append(t2)
                        else:
                            out.append(move5x6_vertex_v1(i2[0], i2[1], b * u, b * u))
                else:
                    out.append(t)
            elif (edge_hex or south_pentagon) and south:
                n = (n_cap_rows + r // 2) if polar_pentagon else (n_cols - 1 - n_cap_rows - r // 2)
                a = min(col, n)
                b = col - a
                if b:
                    if is_pentagon and r <= n_mid_rows // 2 and not polar_pentagon:
                        t = move5x6_vertex_v3(rc[0], rc[1], a * u, a * u)
                        out.append(move5x6_vertex_v1(t[0], t[1], b * u, b * u))
                    else:
                        t = move5x6_vertex_v3(rc[0], rc[1], 0 if (polar_pentagon and r <= n_cap_rows) else a * u, a * u)
                        i2 = _cross5x6_interruption(t[0], t[1], True, polar_pentagon and r >= n_mid_rows // 2)
                        if polar_pentagon:
                            nb = n_cap_rows - r
                            b = min(b, nb); c = col - (a + b)
                            t2 = move5x6_vertex_v1(i2[0], i2[1], b * u, b * u)
                            if c:
                                i2 = _cross5x6_interruption(t2[0], t2[1], True, r >= n_mid_rows // 2)
                                if r >= n_mid_rows // 2:
                                    out.append(move5x6_vertex_v3(i2[0], i2[1], 0, -c * u))
                                else:
                                    out.append(move5x6_vertex_v3(i2[0], i2[1], c * u, 0))
                            else:
                                out.append(t2)
                        else:
                            out.append(move5x6_vertex_v1(i2[0], i2[1], b * u, b * u if is_pentagon else 0))
                else:
                    out.append(move5x6_vertex_v3(rc[0], rc[1], 0 if (polar_pentagon and r <= n_cap_rows) else a * u, a * u))
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, col * u))
        r += 1
        n_cols += 1 if (n_cols & 1) else -1

    # Second cap
    n_cols -= 2
    r = 0
    while r < n_cap_rows - end_cap_skip:
        n = n_cap_rows - 2 * (r + 1) if ((edge_hex or is_pentagon) and r < n_cap_rows // 2) else 0

        # row-start (eC :1202-1256) -- poles use fixed reflection points
        if polar_pentagon:
            if south:
                i2 = [2 + n_cap_rows * u, 4 - n_cap_rows * u]
                rc = move5x6_vertex_v1(i2[0], i2[1], -(r + 1) * 2 * u, -(r + 1) * u)
            else:
                i2 = [3 - n_cap_rows * u, 2 + n_cap_rows * u]
                rc = move5x6_vertex_v1(i2[0], i2[1], (r + 1) * 2 * u, (r + 1) * u)
        elif (edge_hex or north_pentagon) and not south:
            if r < n_cap_rows // 2:
                a = min(r + 1, n_cap_rows // 2)
                b = r + 1 - a
                t = move5x6_vertex_v1(
                    first_cx, first_cy,
                    n_cap_rows * -u + (n_mid_rows >> 1) * u + a * 2 * u,
                    n_cap_rows * u + (n_mid_rows - 1) * u + a * u,
                )
                if b:
                    i2 = _cross5x6_interruption(t[0], t[1], False, False)
                    rc = move5x6_vertex_v1(i2[0], i2[1], b * 2 * u, b * u)
                else:
                    rc = t
            else:
                a = n_cap_rows
                ay = 3 * n_cap_rows - 1 - r
                b = 1 + 2 * (r - n_cap_rows // 2)
                t = move5x6_vertex_v1(first_cx, first_cy, a * u, a * u)
                i2 = _cross5x6_interruption(t[0], t[1], False, False)
                rc = move5x6_vertex_v3(i2[0], i2[1], -ay * u, b * u)
        elif (edge_hex or south_pentagon) and south and r >= n_cap_rows // 2:
            a = n_cap_rows
            ay = 2 * n_cap_rows + r + 1
            b = 1 + 2 * (r - n_cap_rows // 2)
            t = move5x6_vertex_v1(first_cx, first_cy, 0, a * u)
            i2 = _cross5x6_interruption(t[0], t[1], True, False)
            rc = move5x6_vertex_v3(i2[0], i2[1], b * u, ay * u)
        else:
            rc = move5x6_vertex_v3(
                first_cx, first_cy,
                n_cap_rows * -2 * u + ((n_mid_rows - 1) >> 1) * -u + (r + 1) * u,
                n_cap_rows * -1 * u + ((n_mid_rows - 1) >> 1) * u + (r + 1) * 2 * u,
            )

        col = 0
        while col < n_cols:
            if (edge_hex or north_pentagon) and not south:
                a = min(col, n)
                b = col - a
                if a:
                    if b:
                        if is_pentagon:
                            b -= (n_cols // 2 - n // 2) + r + r // 2 + 1
                            if polar_pentagon:
                                t = move5x6_vertex_v1(rc[0], rc[1], a * u, 0)
                                i2 = _cross5x6_interruption(t[0], t[1], False, False)
                            else:
                                i2 = move5x6_vertex_v3(rc[0], rc[1], a * u, 0)
                        else:
                            t = move5x6_vertex_v1(rc[0], rc[1], a * u, 0)
                            i2 = _cross5x6_interruption(t[0], t[1], False, False)
                        if polar_pentagon:
                            out.append(move5x6_vertex_v3(i2[0], i2[1], 0, b * u))
                        else:
                            out.append(move5x6_vertex_v3(i2[0], i2[1], b * u, b * u))
                    else:
                        out.append(move5x6_vertex_v3(rc[0], rc[1], a * u, 0))
                else:
                    out.append(move5x6_vertex_v3(rc[0], rc[1], b * u, b * u))
            elif (edge_hex or south_pentagon) and south:
                a = min(col, n)
                b = col - a
                if a:
                    if polar_pentagon:
                        t = move5x6_vertex_v1(rc[0], rc[1], -a * u, 0)
                    else:
                        mover = move5x6_vertex_v1 if b else move5x6_vertex_v3
                        t = mover(rc[0], rc[1], a * u, a * u)
                    if b:
                        i2 = _cross5x6_interruption(t[0], t[1], True, polar_pentagon)
                        if is_pentagon:
                            b -= (n_cols // 2 - n // 2) + r + r // 2 + 1
                            if polar_pentagon:
                                out.append(move5x6_vertex_v3(i2[0], i2[1], 0, -b * u))
                            else:
                                out.append(move5x6_vertex_v3(i2[0], i2[1], b * u, b * u))
                        else:
                            out.append(move5x6_vertex_v3(i2[0], i2[1], b * u, 0))
                    else:
                        out.append(t)
                else:
                    out.append(move5x6_vertex_v3(rc[0], rc[1], b * u, 0))
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u, col * u))

            if is_pentagon and col == n:
                col = n_cols - n - 1  # Skip interruption
            col += 1
        r += 1
        n_cols -= 3
    return out


def _gen_odd_parent_even_depth(first_cx, first_cy, depth, u, south, edge_hex, is_pentagon=False, polar_pentagon=False):
    """Polar pentagon (Task 4, prepended ``if(polarPentagon)`` branches in all
    three sections). First cap in-row (``:1447-1467``) has an ``r&1`` decrement
    of ``b`` with a ``sgn*u/3`` nudge before the cross. The main-section polar
    in-row (``:1539-1609``) is the file's deepest chain: a triple
    ``maxB``/``c`` cross-cascade whose left/right hemisphere split turns on
    ``r < nCapRows`` vs ``r >= nCapRows`` and carries an ``oddR`` parity
    correction (``:1569-1589``). Row-starts double-cross via the pentagon
    centroid reflection; ``sgn = southPentagon ? -1 : 1`` throughout.

    ``generateOddParentEvenDepth`` (``I3HSubZones.ec:1366-1767``). The first
    cap (``I3HSubZones.ec:1423-1473``) has NO edge-hex OR pentagon branch at
    all -- north/south edge-hex AND non-polar pentagon are all identical to
    interior there (only ``polarPentagon`` is tested, never ``nv``). Main
    row-start/in-row: ``1511-1524``/``1611-1648`` (south row-start only fires
    for ``r != 0``; north has no row-start override, only an in-row one).
    Second cap row-start/in-row: ``1684-1697``/``1734-1737``.

    Non-polar pentagon (Task 3, ``polarPentagon`` hardcoded ``False``):
      - ``skip = nv==5 && r>nCapRows ? r-nCapRows : 0`` shrinks the main
        section's per-row column count for EITHER orientation once past
        ``nCapRows`` (``:1478``).
      - South pentagon gets its OWN row-start branch for ``r > nCapRows``
        (``:1503-1509``, checked BEFORE the south edge-hex/pentagon shared
        branch), offsetting from the pentagon's own centroid rather than the
        edge-hex reflection point; for ``0 < r <= nCapRows`` it shares the
        south edge-hex row-start but skips the explicit cross when ``r >
        nCapRows`` (dead within this sub-range, kept for source fidelity).
        North pentagon has no row-start override (falls to interior default,
        matching north edge-hex's absence of one here too).
      - In-row: north pentagon shares the ``northPentagon || (edgeHex &&
        !southRhombus)`` branch, forcing ``n`` to a pentagon-specific value
        and skipping the cross for ``r > nCapRows`` (``:1613, :1621-1622``);
        south pentagon likewise shares ``(edgeHex||southPentagon) &&
        southRhombus && r``, forcing its own ``n`` and skipping the cross for
        ``r >= nCapRows`` (``:1632, :1640-1641``).
      - Second cap: row-loop bound shrinks by ``endCapSkip = nv==5 ?
        (nCapRows+1)/2 : 0`` (``:1417``, ``:1657``) for EITHER orientation.
        South pentagon gets its own dedicated row-start branch (``:1677-
        1682``, checked before edge-hex); north pentagon has none (falls to
        the ``else`` default, like north edge-hex's ABSENCE of an override
        here -- edge-hex only overrides when south). In-row: pentagon (either
        orientation) computes its own ``a``/``b`` from a per-row ``n`` and a
        shift-adjusted ``b`` UNCONDITIONALLY (``:1709-1715``, ahead of the
        edge-hex dispatch), then an unconditional ``if(nv==5 && col==n) col =
        nCols-n-1`` column-skip (``:1755``)."""
    north_pentagon = is_pentagon and not south
    south_pentagon = is_pentagon and south
    sgn = -1 if south_pentagon else 1
    out = []
    n_cap_rows = 3 ** ((depth - 2) // 2)
    n_mid_rows = 2 * n_cap_rows + 1
    end_cap_skip = (n_cap_rows + 1) // 2 if is_pentagon else 0
    n_cols = 1

    # First cap -- interior/edge/non-polar-pentagon are identical; polar prepended (eC :1433-1467)
    r = 0
    while r < n_cap_rows:
        if polar_pentagon:
            rc = move5x6_vertex_v1(first_cx, first_cy, sgn * r * u, sgn * r * u)
        else:
            rc = move5x6_vertex_v1(first_cx, first_cy, 0, r * u)
        for col in range(n_cols):
            if polar_pentagon:
                n = r // 2 + r
                a = min(col, n); b = col - a
                t = move5x6_vertex_v1(rc[0], rc[1], sgn * -a * u / 3, sgn * -2 * a * u / 3)
                if b:
                    if r & 1:
                        b -= 1
                        t = move5x6_vertex_v1(t[0], t[1], sgn * -u / 3, sgn * -u / 3)
                    i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, not south_pentagon)
                    out.append(move5x6_vertex_v1(i2[0], i2[1], sgn * (-b * 2 - (r & 1)) * u / 3, sgn * (-b - (r & 1)) * u / 3))
                else:
                    out.append(t)
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u / 3, -col * u / 3))
        r += 1
        n_cols += 3

    # Main section
    r = 0
    while r < n_mid_rows:
        skip = (r - n_cap_rows) if (is_pentagon and r > n_cap_rows) else 0

        if polar_pentagon:
            t = move5x6_vertex_v1(first_cx, first_cy, sgn * n_cap_rows * u, sgn * n_cap_rows * u)
            i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, south_pentagon)
            if r:
                a = r >> 1; b = r & 1
                rc = move5x6_vertex_v1(i2[0], i2[1], sgn * a * u + sgn * b * 2 * u / 3, sgn * a * u + sgn * b * u / 3)
            else:
                rc = i2
        elif south_pentagon and r > n_cap_rows:
            a = r // 2
            b = r & 1
            t = move5x6_vertex_v1(first_cx, first_cy, n_cap_rows * u, n_cap_rows * u)  # pentagon centroid
            rc = move5x6_vertex_v1(t[0], t[1], a * u + b * u / 3, n_cap_rows * u - b * u / 3)
        elif (edge_hex or south_pentagon) and south and r:
            n = r + r // 2
            t = move5x6_vertex_v1(
                first_cx, first_cy,
                (n_cap_rows + r) * 2 * u / 3 - (2 * n_cap_rows - r) * u / 3,
                (n_cap_rows + r) * u / 3 + (2 * n_cap_rows - r) * u / 3,
            )
            if is_pentagon and r > n_cap_rows:
                i2 = t
            else:
                i2 = _cross5x6_interruption(t[0], t[1], True, False)
            rc = move5x6_vertex_v1(i2[0], i2[1], n * u / 3, n * 2 * u / 3)
        else:
            rc = move5x6_vertex_v1(
                first_cx, first_cy,
                (r >> 1) * u + (r & 1) * 2 * u / 3,
                (n_cap_rows + (r >> 1)) * u + (r & 1) * u / 3,
            )

        for col in range(n_cols - skip):
            if polar_pentagon:
                r2 = (2 * n_cap_rows - r) if r > n_cap_rows else r
                n = r2 + r2 // 2
                a = min(col, n); b = col - a
                if a:
                    t = move5x6_vertex_v1(rc[0], rc[1], sgn * a * u / 3, sgn * -a * u / 3)
                else:
                    t = rc
                if b:
                    max_b = (n_cap_rows - r + (n_cap_rows - r) // 2) if r < n_cap_rows else (r - n_cap_rows)
                    c = 0
                    crossing_left = (r >= n_cap_rows) if south_pentagon else (r < n_cap_rows)
                    if b > max_b:
                        c = b - max_b; b = max_b
                    i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, crossing_left)
                    if r >= n_cap_rows:
                        cen = move5x6_vertex_v1(i2[0], i2[1], sgn * b * 2 * u / 3, sgn * b * u / 3)
                    else:
                        cen = move5x6_vertex_v1(i2[0], i2[1], sgn * -b * u / 3, sgn * -b * 2 * u / 3)
                    if r < n_cap_rows:
                        if c:
                            odd_r = r & 1
                            if not odd_r:
                                c -= 1
                                cen = move5x6_vertex_v1(cen[0], cen[1], sgn * -u / 3, sgn * -u / 3)
                            i2 = _cross5x6_interruption(cen[0], cen[1], south_pentagon, north_pentagon)
                            if c > max_b:
                                a = min(c, max_b); b = c - a
                                t = move5x6_vertex_v1(i2[0], i2[1], sgn * (-a * 2 - (1 - odd_r)) * u / 3, sgn * (-a - (1 - odd_r)) * u / 3)
                                i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, north_pentagon)
                                cen = move5x6_vertex_v1(i2[0], i2[1], sgn * -b * u / 3, sgn * b * u / 3)
                            else:
                                cen = move5x6_vertex_v1(i2[0], i2[1], sgn * (-c * 2 - (1 - odd_r)) * u / 3, sgn * (-c - (1 - odd_r)) * u / 3)
                    else:
                        i2 = cen
                        if c > max_b:
                            a = min(c, max_b); b = c - a
                            t = move5x6_vertex_v1(i2[0], i2[1], sgn * a * u / 3, sgn * a * 2 * u / 3)
                            i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, crossing_left)
                            cen = move5x6_vertex_v1(i2[0], i2[1], sgn * -b * u / 3, sgn * b * u / 3)
                        else:
                            cen = move5x6_vertex_v1(i2[0], i2[1], sgn * c * u / 3, sgn * c * 2 * u / 3)
                    out.append(cen)
                else:
                    out.append(t)
            elif north_pentagon or (edge_hex and not south):
                if north_pentagon and r > n_cap_rows:
                    n = n_cap_rows + n_cap_rows // 2 - (r - n_cap_rows) // 2
                else:
                    n = 3 * n_cap_rows - r - (r + 1) // 2
                a = min(col, n)
                b = col - a
                t = move5x6_vertex_v1(rc[0], rc[1], a * u / 3, -a * u / 3)
                if b:
                    if north_pentagon and r > n_cap_rows:
                        i2 = t
                    else:
                        i2 = _cross5x6_interruption(t[0], t[1], False, False)
                    out.append(move5x6_vertex_v1(i2[0], i2[1], b * 2 * u / 3, b * u / 3))
                else:
                    out.append(t)
            elif (edge_hex or south_pentagon) and south and r:
                if south_pentagon and r > n_cap_rows:
                    n = n_cap_rows + n_cap_rows // 2 - (r - n_cap_rows) // 2
                else:
                    n = r + r // 2
                a = min(col, n)
                b = col - a
                t = move5x6_vertex_v1(rc[0], rc[1], -a * u / 3, -a * 2 * u / 3)
                if b:
                    if south_pentagon and r >= n_cap_rows:
                        i2 = t
                    else:
                        i2 = _cross5x6_interruption(t[0], t[1], True, True)
                    out.append(move5x6_vertex_v1(i2[0], i2[1], b * u / 3, -b * u / 3))
                else:
                    out.append(t)
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u / 3, -col * u / 3))
        r += 1
        n_cols += 1 if (n_cols & 1) else -1

    # Second cap
    n_cols -= 2
    r = 0
    while r < n_cap_rows - end_cap_skip:
        n = 0
        if is_pentagon:
            n = 0 if r >= n_cap_rows // 2 else (n_cap_rows - 2 * (r + 1))

        if polar_pentagon:
            r2 = n_mid_rows - 1; a = r2 >> 1; b = r2 & 1
            t = move5x6_vertex_v1(first_cx, first_cy, sgn * n_cap_rows * u, sgn * n_cap_rows * u)
            i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, south_pentagon)
            t = move5x6_vertex_v1(i2[0], i2[1], sgn * a * u + sgn * b * 2 * u / 3, sgn * a * u + sgn * b * u / 3)
            i2 = _cross5x6_interruption(t[0], t[1], south_pentagon, south_pentagon)
            rc = move5x6_vertex_v1(i2[0], i2[1], sgn * (r + 1) * u, sgn * (r + 1) * u)
        elif south_pentagon:
            t = move5x6_vertex_v1(first_cx, first_cy, n_cap_rows * u, n_cap_rows * u)  # pentagon centroid
            rc = move5x6_vertex_v1(t[0], t[1], (n_mid_rows // 2) * u, (n_cap_rows - r - 1) * u)
        elif edge_hex and not south:
            t = move5x6_vertex_v1(first_cx, first_cy, 3 * n_cap_rows * u / 3, 3 * n_cap_rows * 2 * u / 3)
            i2 = _cross5x6_interruption(t[0], t[1], False, False)
            rc = move5x6_vertex_v1(i2[0], i2[1], (r + 1) * u, (r + 1) * u)
        elif edge_hex and south:
            t = move5x6_vertex_v1(first_cx, first_cy, 3 * n_cap_rows * 2 * u / 3, 3 * n_cap_rows * u / 3)
            i2 = _cross5x6_interruption(t[0], t[1], True, False)
            rc = move5x6_vertex_v1(i2[0], i2[1], 3 * n_cap_rows * u / 3, 3 * n_cap_rows * 2 * u / 3 - (r + 1) * u)
        else:
            rc = move5x6_vertex_v1(
                first_cx, first_cy,
                ((n_mid_rows + 1) // 2 + r) * u,
                (n_cap_rows + ((n_mid_rows - 1) >> 1)) * u,
            )

        col = 0
        while col < n_cols:
            a = b = 0
            if is_pentagon:
                a = min(col, n)
                b = col - a
                if b:
                    b -= (n_cols // 2 - n // 2) + r + r // 2 + 1

            if polar_pentagon:
                t = move5x6_vertex_v1(rc[0], rc[1], sgn * a * 2 * u / 3, sgn * a * u / 3)
                out.append(move5x6_vertex_v1(t[0], t[1], sgn * b * u / 3, sgn * b * 2 * u / 3))
            elif north_pentagon:
                t = move5x6_vertex_v1(rc[0], rc[1], a * u / 3, -a * u / 3)
                out.append(move5x6_vertex_v1(t[0], t[1], b * 2 * u / 3, b * u / 3))
            elif south_pentagon:
                t = move5x6_vertex_v1(rc[0], rc[1], -a * u / 3, -2 * a * u / 3)
                out.append(move5x6_vertex_v1(t[0], t[1], b * u / 3, -b * u / 3))
            elif edge_hex and not south:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * 2 * u / 3, col * u / 3))
            elif edge_hex and south:
                out.append(move5x6_vertex_v1(rc[0], rc[1], -col * u / 3, -2 * col * u / 3))
            else:
                out.append(move5x6_vertex_v1(rc[0], rc[1], col * u / 3, -col * u / 3))

            if is_pentagon and col == n:
                col = n_cols - n - 1  # Skip interruption
            col += 1
        r += 1
        n_cols -= 3
    return out


def _i3h_sub_zone_centroids(level_i9r, root, rix, sub_hex, depth, sz_level=None):
    """Dispatch (``iterateI3HSubZones``, ``I3HSubZones.ec:1802-1847``) into the
    4 parity generators. All four cell classes are handled: interior,
    edge-hexagon and non-polar-pentagon (Tasks 1-3) plus the 2 polar pentagons
    (the North/South icosahedron vertices, ``root >= 10``; Task 4). The
    ``polar_pentagon`` flag is threaded to each generator and
    :func:`_i3h_first_sub_zone_centroid`, which carry the eC's
    ``if(polarPentagon)`` branches; verified 0-mismatch vs pydggal by direct
    enumeration of both poles x both ``sub_hex`` x levels 0-3 x depths 1-9 x all
    3 grids (see ``.superpowers/sdd/task-4-report.md``).

    ``edge_hex`` requires ``sub_hex in (0, 1)`` (``I3HSubZones.ec:1817``:
    ``(subHex == 0 || subHex == 1)``) as well as the rix-boundary test -- a C/D
    sub-hex (2/3) at a boundary ``rix`` is NOT edge-hex even though its
    underlying grid cell touches the interruption, because the C/D offset
    (:func:`_i3h_centroid_xy`) pushes its centroid strictly into the cell's
    interior, away from the corner the interruption passes through. Empirically
    confirmed against pydggal: ``B6-1-D`` (an odd-level ``D`` cell sitting at
    the same boundary ``rix`` as the edge-hex ``B6-1-A``) generates correct
    sub-zones via the plain interior path -- forcing ``edge_hex=True`` for it
    produces a mismatch.

    ``is_pentagon`` (``rix == 0 and sub_hex <= 1``, root<10 non-polar --
    the 10 icosahedron vertices, matching ``nv == zone.nPoints == 5`` in the
    eC's dispatcher) is passed through to all 4 generators and
    :func:`_i3h_first_sub_zone_centroid`, which each carry their own
    ``polarPentagon``-hardcoded-``False`` pentagon branches (Task 3; verified
    0-mismatch vs pydggal by direct enumeration of all 10 cells, several
    depths, all 3 grids -- see ``docs/superpowers/sdd/task-3-report.md``)."""
    is_pentagon = is_pentagon_cell(rix, sub_hex)
    polar_pentagon = root >= 10  # the 2 icosahedron poles (root 10 = North, 11 = South)
    south = bool(root & 1)
    edge_hex = _i3h_is_edge_hex(level_i9r, root, rix, sub_hex)

    odd_parent = sub_hex > 0
    odd_depth = depth & 1
    if sz_level is None:  # caller (_i3h_sub_zones) may already have it
        sz_level = _i3h_sz_level(level_i9r, sub_hex, depth)
    u = 1.0 / (3 ** (sz_level // 2))
    fcx, fcy = _i3h_first_sub_zone_centroid(level_i9r, root, rix, sub_hex, depth, edge_hex)

    if odd_parent:
        if odd_depth:
            return _gen_odd_parent_odd_depth(fcx, fcy, depth, u, south, edge_hex, is_pentagon, polar_pentagon)
        return _gen_odd_parent_even_depth(fcx, fcy, depth, u, south, edge_hex, is_pentagon, polar_pentagon)
    if odd_depth:
        return _gen_even_parent_odd_depth(fcx, fcy, depth, u, south, edge_hex, is_pentagon, polar_pentagon)
    return _gen_even_parent_even_depth(fcx, fcy, depth, u, south, edge_hex, is_pentagon, polar_pentagon)


def _i3h_sub_zones(value, depth):
    """Quantize the ordered sub-zone centroids into packed I3H zone ints
    (matches ``dggrs.ec:195-222``'s generic ``getSubZones`` composed with the
    I3H-specific centroid generator)."""
    level_i9r, root, rix, sub_hex = unpack_i3h(value)
    sz_level = _i3h_sz_level(level_i9r, sub_hex, depth)
    centroids = _i3h_sub_zone_centroids(level_i9r, root, rix, sub_hex, depth, sz_level)
    # _i3h_from_centroid returns None for its nullZone guard (see :593). Every
    # other call site handles that -- quantize -> (NULL_ZONE, []), _i3h_get_children
    # -> skip, _i3h_get_neighbor -> None -- but the sub-zone path splatted it
    # into pack_i3h(*None) and raised TypeError. Mirror the sibling behaviour.
    out = []
    for (cx, cy) in centroids:
        cell = _i3h_from_centroid(sz_level, cx, cy)
        out.append(NULL_ZONE if cell is None else pack_i3h(*cell))
    return out


# --------------------------------------------------------------------------- #
# The Topology facade
# --------------------------------------------------------------------------- #
class HexAperture3Topology:
    """Aperture-3 hexagonal topology on the 5x6 ISEA grid — a ``Topology`` impl
    for the I3H family (ISEA3H/IVEA3H).

    Unlike :class:`~py4dggs.topologies.hex_a7.HexAperture7Topology`, the whole
    quantized cell (level, root, rhombus index, sub-hex) lives in ``base``
    (packed by :func:`py4dggs.indexings.i3h.pack_i3h`); ``digits`` is always ``[]``.
    """

    aperture = 3

    def quantize(self, geom: Any, p: PlanarPoint, res: int) -> tuple[int, list[int]]:
        """Quantize planar point ``p`` to an I3H cell at resolution ``res``
        (``I3HZone::fromCentroid``, see :func:`_i3h_from_centroid`). Like
        ``hex_a7``'s ``quantize``, ``p.face`` is unused — the cell is determined
        from ``(p.x, p.y)`` alone. Returns ``(packed_i3h_int, [])``: the whole
        cell packs into ``base`` (see the module docstring), so ``digits`` is
        always empty."""
        zone = _i3h_from_centroid(res, p.x, p.y)
        if zone is None:
            return NULL_ZONE, []
        return pack_i3h(*zone), []

    def planar_centroid(self, geom: Any, base: int, digits: list[int]) -> PlanarPoint:
        """The cell's centroid as a planar ``PlanarPoint`` in 5x6 space
        (``I3HZone::centroid``, ``RI3H.ec:2138-2168``). ``base`` carries the
        whole packed I3H cell (see module docstring); ``digits`` is unused
        (always ``[]``). ``face`` is the ``-1`` sentinel ("derive from x,y");
        the Grid runs the projection's inverse to get lat/lon."""
        cx, cy = _i3h_centroid_xy(*unpack_i3h(base))
        return PlanarPoint(face=-1, x=cx, y=cy)

    def planar_vertices(self, geom: Any, base: int, digits: list[int]) -> list[PlanarPoint]:
        """The cell's boundary vertices (5 for a pentagon, 6 for a hexagon) as
        planar ``PlanarPoint``s in 5x6 space, each with the ``-1`` face sentinel
        (``I3HZone::getVertices`` plus the canonicalizing wrap its callers apply,
        ``RI3H.ec:1674-1790`` and ``RI3H.ec:337-386``). ``base`` carries the whole
        packed I3H cell (see module docstring); ``digits`` is unused (always
        ``[]``)."""
        level_i9r, root, rix, sub_hex = unpack_i3h(base)
        raw = _i3h_vertices(level_i9r, root, rix, sub_hex)
        return [PlanarPoint(face=-1, x=cx, y=cy) for cx, cy in (canonicalize5x6(vx, vy) for vx, vy in raw)]

    def neighbors(self, geom: Any, base: int, digits: list[int]) -> list[tuple[int, list[int]]]:
        """The cell's exact topological neighbours (``I3HZone::getNeighbors``, see
        :func:`_i3h_get_neighbors`). This OPTIONAL ``Topology`` method makes
        ``Grid.neighbors`` use the exact aperture-3 adjacency instead of the
        grid-agnostic edge k-ring (which overshoots for ~0.01% of aperture-3 cells
        at apices / interruption seams). ``base`` carries the whole packed I3H cell
        (``digits`` is always ``[]``); each returned neighbour is ``(packed_int,
        [])`` to match the aperture-3 empty-digits convention."""
        return [(nb, []) for nb in _i3h_get_neighbors(*unpack_i3h(base))]

    # --- exact non-congruent hierarchy (A2): optional Topology overrides ---
    def parents(self, geom: Any, base: int, digits: list[int]) -> list[tuple[int, list[int]]]:
        """The cell's 1-or-3 geometric parents (``getParents``, RI3H.ec:1247-1329)."""
        return [(p, []) for p in _i3h_get_parents(base)]

    def children(self, geom: Any, base: int, digits: list[int]) -> list[tuple[int, list[int]]]:
        """The cell's geometric children (``getChildren``, RI3H.ec:2044-2085) —
        6 for a pentagon, 7 for a hexagon."""
        return [(c, []) for c in _i3h_get_children(base)]

    def centroid_parent(self, geom: Any, base: int, digits: list[int]):
        """The centroid parent (``centroidParent``, RI3H.ec:1203-1224), or None."""
        cp = _i3h_centroid_parent(base)
        return None if cp is None else (cp, [])

    def is_centroid_child(self, geom: Any, base: int, digits: list[int]) -> bool:
        """Whether this cell has a single parent (``isCentroidChild``, RI3H.ec:2170-2196)."""
        return _i3h_is_centroid_child(*unpack_i3h(base))

    # --- sub-zones (A3): optional Topology overrides, matching the
    # (geom, base, digits) -> list[(int, digits)] shape every hierarchy
    # method above uses (parents/children/centroid_parent/is_centroid_child) ---
    def count_sub_zones(self, geom: Any, base: int, digits: list[int], relative_depth: int) -> int:
        _, root, rix, sub_hex = unpack_i3h(base)
        is_pentagon = is_pentagon_cell(rix, sub_hex)
        nv = 5 if (is_pentagon or root >= 10) else 6
        return _i3h_count_sub_zones(nv, relative_depth)

    def first_sub_zone(self, geom: Any, base: int, digits: list[int], relative_depth: int) -> tuple[int, list[int]]:
        level_i9r, root, rix, sub_hex = unpack_i3h(base)
        sz_level = _i3h_sz_level(level_i9r, sub_hex, relative_depth)
        cx, cy = _i3h_first_sub_zone_centroid(level_i9r, root, rix, sub_hex, relative_depth)
        cell = _i3h_from_centroid(sz_level, cx, cy)  # None == nullZone, see _i3h_sub_zones
        return (NULL_ZONE if cell is None else pack_i3h(*cell), [])

    def sub_zones(self, geom: Any, base: int, digits: list[int], relative_depth: int) -> list[tuple[int, list[int]]]:
        return [(v, []) for v in _i3h_sub_zones(base, relative_depth)]
