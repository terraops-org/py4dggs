"""Aperture-7 hexagonal topology on the 5x6 ISEA planar grid.

This is the topology half of IGEO7/ISEA7H. It turns a projection's planar point
``(x, y)`` in DGGAL's oblique 5x6 grid into a Z7 address ``(base cell, direction
digits)`` and reconstructs planar cell geometry (centroid + 5/6 vertices) from
that address. The geographic lift (centroid/vertices in lat/lon) and neighbours
live on the *Grid*, so this stays projection-agnostic above the shared 5x6 layout.

Provenance (the *why* of each block is cited inline against these sources):
  - ``dggrs/RI7H.ec``     — the I7H zone, ``fromCentroid`` quantization (§6),
                            ``getCentroid``/``getVertices`` geometry, the 5x6
                            wrapping/interruption helpers (§11), ``isEdgeHex``,
                            ``calcCandidateParent``, ``getPrimaryChildren`` (§9).
  - ``dggrs/RI7H_Z7.ec``  — the I7H<->Z7 conversion: ``from7H``/``to7H``,
                            ``getLevelRotationOffset``, ``getChildPosition``,
                            the pentagon child-position (de)adjustments.
Spec cross-reference: ``DGGAL-ALGORITHMS.md`` §6 (fromCentroid), §9 (children),
§11 (5x6 wrapping).

It is re-expressed from ``igeo7/hex7.py``. Readability lives in the structure — a
:class:`HexAperture7Topology` facade implementing the ``Topology`` protocol, named
free-function helpers grouped by concept, and these docstrings. The *arithmetic* is
held byte-for-byte faithful to that oracle (same expressions, groupings, epsilons,
evaluation order, and the ``_jsround`` half-up rounding), because the acceptance
gate is bit-identity: re-associating a single float op or merging a branch would
drift the last bit and fail the differential fuzz.

ONE DELIBERATE DIVERGENCE (2026-07-02): :func:`add_non_polar_base_vertices` fixes a
real bug that ``igeo7/hex7.py`` (and the JS port it descends from) carries — the
eC's i1/i2 interruption-frame reconciliation (``RI7H.ec:1991-2004``) was dropped,
so ~0.3% of cells (those whose boundary spans a rhombus interruption) got wrong /
degenerate "far-side" vertices, which cascade into wrong k-ring neighbours. py4dggs
ports the eC faithfully instead, so for those cells its vertices/neighbours match
DGGAL (pydggal) — the ground truth — NOT the buggy igeo7-py. Everything else
(quantize, centroid, hierarchy, text-id) stays bit-identical to igeo7-py.

Zero runtime dependencies. Historically this arithmetic was verified against a
frozen Python port (``igeo7-py``, a separate single-grid sibling repo) as a dev-
only test dependency; that sibling-repo dependency was removed (2026-07-07) --
verification now runs entirely against **pydggal** (the DGGAL engine itself),
like every other grid in this library (see ``tests/test_isea7h_fuzz.py``,
``tests/test_conformance.py``). No test in this repo imports ``igeo7`` anymore.
"""
from __future__ import annotations

import math
from typing import Any, NamedTuple

from py4dggs.types import PlanarPoint

# --------------------------------------------------------------------------- #
# I7H <-> Z7 mapping tables (RI7H_Z7.ec)
# --------------------------------------------------------------------------- #
C_MAP = [0, 3, 1, 5, 4, 6, 2]          # I7H child position -> Z7 digit
INV_C_MAP = [0, 2, 6, 1, 4, 3, 5]      # Z7 digit -> I7H child position
ROOT_MAP = [1, 6, 2, 7, 3, 8, 4, 9, 5, 10, 0, 11]      # I7H root rhombus -> Z7 base cell
INV_ROOT_MAP = [10, 0, 2, 4, 6, 8, 1, 3, 5, 7, 9, 11]  # Z7 base cell -> I7H root rhombus

# Powers of 7 (cells along a rhombus edge at l49r); the table covers the resolutions
# the grid reaches, with a rounded fallback past it.
POW7 = [1, 7, 49, 343, 2401, 16807, 117649, 823543, 5764801, 40353607, 282475249]

# The first Z7 level with no geometry. The 64-bit packing holds 20 direction
# digits, so level 20 is representable, but `to7H` (RI7H_Z7.ec:348-353) "does
# not support level 20 zones" and returns nullZone — hence DGGAL's highest real
# level is 19 (`Z7Indexing.max_resolution`). Kept here rather than imported from
# the indexing: it is a property of the Z7->7H *geometry* conversion, which is
# this module's job. See `HexAperture7Topology.is_null_geometry`.
_NULL_GEOMETRY_LEVEL = 20


# --------------------------------------------------------------------------- #
# Rounding & power helpers
# --------------------------------------------------------------------------- #
def _jsround(x: float) -> int:
    """Round half toward +infinity (``floor(x + 0.5)``), matching JS ``Math.round``
    and the eC. Python's built-in ``round`` is banker's rounding, which would pick
    a different cell on exact half-boundaries — so the whole pipeline uses this."""
    return math.floor(x + 0.5)


def pow7(n: int) -> int:
    """7**n — table lookup within range, rounded ``7**n`` beyond it."""
    return POW7[n] if n < len(POW7) else round(7 ** n)


# --------------------------------------------------------------------------- #
# The I7H zone + its derived predicates (RI7H.ec I7HZone)
# --------------------------------------------------------------------------- #
class _Z(NamedTuple):
    """An I7H zone. ``l49r`` is the aperture-49 resolution (two Z7 levels per step);
    ``root`` is the rhombus (0..9) or a polar sentinel (0xA north, 0xB south);
    ``row``/``col`` index within the rhombus; ``sub_hex`` is the odd-level child
    selector (0 = even level). ``level = 2*l49r + (sub_hex > 0)``."""
    l49r: int
    root: int
    row: int
    col: int
    sub_hex: int


def zone_level(z: _Z) -> int:
    """Z7 level of the zone (two levels per aperture-49 step, +1 if odd sub-hex)."""
    return 2 * z.l49r + (1 if z.sub_hex > 0 else 0)


def zone_npoints(z: _Z) -> int:
    """Number of cell vertices: 5 for the rhombus-corner pentagon, else 6."""
    if z.sub_hex > 1:
        return 6
    return 5 if (z.row == 0 and z.col == 0) else 6


def zone_eq(a, b) -> bool:
    """Zone equality that treats ``None`` as a distinct empty value."""
    if a is None or b is None:
        return a is b
    return a == b


def zone_is_edge_hex(z: _Z) -> bool:
    """True for an even-level hexagon lying on a rhombus seam (north edge of a
    north rhombus, south edge of a south rhombus) — the ``isEdgeHex`` property
    (RI7H.ec). Edge hexagons need the rotation fix-ups in child enumeration."""
    if zone_npoints(z) != 6:
        return False
    if z.sub_hex != 0:
        return False
    if z.root & 1:
        return z.col == 0   # South
    return z.row == 0       # North


# --------------------------------------------------------------------------- #
# 5x6 coordinate wrapping / interruption
# (eC 5x6 helpers live in ``projections/ri5x6.ec``, NOT ``dggrs/RI7H.ec``; spec §11)
# --------------------------------------------------------------------------- #
def canonicalize5x6(cx, cy):
    """Fold a 5x6 planar point into the canonical fundamental domain, handling the
    polar diagonals and the rhombus interruptions. Epsilons disambiguate points
    sitting exactly on a seam.

    NB: this ports the *inline* domain-fold the eC does inside ``fromCentroid``
    (``RI7H.ec:1439-1458``), NOT the eC function that happens to share this name
    (``canonicalize5x6``, ``ri5x6.ec:1562``) — that one is a different, more
    elaborate routine (pole-snap to {1,0}/{4,6}, ``cross5x6Interruption``, the
    ``x==5 -> {0, y-5}`` case) used on the vertices path. Same name, different job."""
    x, y = cx, cy
    if abs(x - y - 1) < 1e-10:
        return [x, y]  # north pole diagonal
    if abs(y - x - 2) < 1e-10:
        return [x, y]  # south pole diagonal
    if y < -1e-11 and x > -1e-11:
        x -= y; y = 0
    elif math.floor(x + 1e-11) > math.floor(y + 1e-11):
        iy = min(5, math.floor(y + 1e-11))
        x += (iy + 1 - y); y = iy + 1
    elif math.floor(y + 1e-11) - math.floor(x + 1e-11) > 1:
        ix = min(4, math.floor(x + 1e-11))
        y += (ix + 1 - x); x = ix + 1
    elif x < -1e-11 or y < -1e-11:
        x += 5; y += 5
    return [x, y]


def point_line_side(px, py, ax, ay, bx, by):
    """Signed side of point ``(px,py)`` relative to the directed line A->B
    (positive = left). Used to break ties when ``fromCentroid`` lands a point near
    a hex edge."""
    dx = bx - ax; dy = by - ay
    return dy * (px - ax) - dx * (py - ay)


def move5x6_vertex(cx, cy, dx, dy):
    """Offset a vertex from centroid ``(cx,cy)`` by ``(dx,dy)``, re-expressing the
    offset across a rhombus interruption when it crosses one. Ports
    ``move5x6Vertex`` (``ri5x6.ec:1388``) — NOT ``move5x6Vertex2`` (``ri5x6.ec:1470``),
    a genuinely different routine (with ``crossEarly``/dent logic) that the eC
    itself warns "does not generate correct geometry at level 2" (RI7H.ec:2954).
    The four ``ivx/ivy`` branches mirror the four crossing directions in the eC
    and intentionally share bodies — left as-is for faithfulness."""
    vx = cx + dx
    vy = cy + dy
    icx = math.floor(cx + 1e-11)
    icy = math.floor(cy + 1e-11)
    sgn_dx = 1 if dx > 0 else (-1 if dx < 0 else 1)
    sgn_dy = 1 if dy > 0 else (-1 if dy < 0 else 1)
    ivx = math.floor(cx + dx - sgn_dx * 1e-11)
    ivy = math.floor(cy + dy - sgn_dy * 1e-11)

    if (((ivx != icx and abs(vy - ivy) > 1e-11) or
         (ivy != icy and abs(vx - ivx) > 1e-11)) and
            (ivy - ivx > 1 or ivy < ivx)):
        fcx = cx - icx
        fcy = cy - icy
        if ivx < icx:
            vx = icx - fcy + dx - dy; vy = icy + dx
        elif ivx > icx:
            vx = icx - fcy + dx - dy; vy = icy + dx
        elif ivy < icy:
            vx = icx + dy; vy = icy - fcx - dx + dy
        elif ivy > icy:
            vx = icx + dy; vy = icy - fcx - dx + dy
    return [vx, vy]


def rotate5x6_offset(dx, dy, clockwise):
    """Rotate a 5x6 offset vector by 60 degrees (spec §11). The oblique basis makes
    a 60 deg turn a pure integer-combination of the components."""
    if clockwise:
        return [dx - dy, dx]      # 60 deg CW
    return [dy, dy - dx]          # 60 deg CCW


def crosses5x6_interruption_v2(c_in_x, c_in_y, dx, dy):
    """Walk an offset ``(dx,dy)`` from a 5x6 point and, if it crosses a rhombus
    interruption (the icosahedron's cut edges), return the crossing in/out points
    and hemisphere. Returns ``None`` if it stays within one rhombus. Ports
    ``crosses5x6InterruptionV2`` (``ri5x6.ec:1611``); drives vertex tracing in
    :func:`add_non_polar_base_vertices`."""
    cx0, cy0 = c_in_x, c_in_y
    if cx0 < 0 and cy0 < 1 + 1e-11:
        cx0 += 5; cy0 += 5

    cdx, cdy = cx0, cy0
    north = cdx - cdy - 1e-11 > 0
    if north:
        cdx -= 1e-11; cdy += 1e-11
    else:
        cdx += 1e-11; cdy -= 1e-11

    if cdx < 0 and cdy < 1 + 1e-11:
        cdx += 5; cdy += 5; cx0 += 5; cy0 += 5
    if cdx > 5 and cdy > 5 - 1e-11:
        cdx -= 5; cdy -= 5; cx0 -= 5; cy0 -= 5

    icx = math.floor(cdx)
    icy = math.floor(cdy)

    px = max(icx - cx0, dx) if dx < 0 else min(icx + 1 - cx0, dx)
    py = max(icy - cy0, dy) if dy < 0 else min(icy + 1 - cy0, dy)

    if dx and dy:
        pkx = px / dx; pky = py / dy
        if pkx < pky:
            py = pkx * dy
        elif pky < pkx:
            px = pky * dx

    cx0 += px
    cy0 += py

    if abs(dx - px) < 1e-11 and abs(dy - py) < 1e-11:
        return None

    sgn_dx = 1 if dx > 0 else (-1 if dx < 0 else 1)
    sgn_dy = 1 if dy > 0 else (-1 if dy < 0 else 1)
    nx = math.floor(cx0 + 1e-11 * sgn_dx)
    ny = math.floor(cy0 + 1e-11 * sgn_dy)

    if (nx != icx or ny != icy) and (nx > icx or abs(dx - px) > 1e-11 or abs(dy - py) > 1e-11):
        root = icx + icy
        i_src = i_dst = None
        in_north = None

        if not (root & 1):  # North rhombus
            if ny == icy and nx == icx + 1:
                iy = math.floor(cx0 - 1 + 1e-11)
                i_src = [cx0, cy0]
                i_dst = [iy + 2 - (cy0 - iy), cx0]
                in_north = True
            elif nx == icx and ny == icy - 1:
                ix = math.floor(cy0 + 1e-11)
                i_src = [cx0, cy0]
                i_dst = [cy0, ix - (cx0 - ix)]
                in_north = True
        else:  # South rhombus
            if nx == icx and ny == icy + 1:
                ix = math.floor(cy0 - 2 + 1e-11)
                i_src = [cx0, cy0]
                i_dst = [cy0 - 1, ix + 3 - (cx0 - ix)]
                in_north = False
            elif ny == icy and nx == icx - 1:
                iy = math.floor(cx0 + 1 + 1e-11)
                i_src = [cx0, cy0]
                i_dst = [iy - 1 - (cy0 - iy), cx0 + 1]
                in_north = False

        if i_src:
            if i_dst[0] < 0 and i_dst[1] < 1 + 1e-11:
                i_dst[0] += 5; i_dst[1] += 5
            return {"i_src": i_src, "i_dst": i_dst, "in_north": in_north}
    return None


# --------------------------------------------------------------------------- #
# Hex vertex offset tables (RI7H.ec getVertices) — in units of 1/(7*p)
# --------------------------------------------------------------------------- #
EVEN_HEX_A = 7 / 3
EVEN_HEX_B = 14 / 3
ODD_HEX_A = 4 / 3
ODD_HEX_B = 5 / 3
ODD_HEX_C = 1 / 3

EVEN_HEX_VERTS = [
    [-EVEN_HEX_A, -EVEN_HEX_B], [-EVEN_HEX_B, -EVEN_HEX_A], [-EVEN_HEX_A, EVEN_HEX_A],
    [EVEN_HEX_A, EVEN_HEX_B], [EVEN_HEX_B, EVEN_HEX_A], [EVEN_HEX_A, -EVEN_HEX_A],
]
ODD_HEX_VERTS = [
    [-ODD_HEX_A, -ODD_HEX_B], [-ODD_HEX_B, -ODD_HEX_C], [-ODD_HEX_C, ODD_HEX_A],
    [ODD_HEX_A, ODD_HEX_B], [ODD_HEX_B, ODD_HEX_C], [ODD_HEX_C, -ODD_HEX_A],
]


# --------------------------------------------------------------------------- #
# Zone centroid & containment (RI7H.ec getCentroid / containsPoint)
# --------------------------------------------------------------------------- #
def zone_centroid(z: _Z):
    """The zone's centroid in 5x6 space (``I7HZone::getCentroid``). Even levels sit
    on the rhombus lattice; odd levels are nudged by one of six sub-hex offsets,
    re-expressed across an interruption via :func:`move5x6_vertex`. Edge hexagons
    and the south-pole corner rotate the sub-hex index first (matching the children
    enumeration)."""
    l49r = z.l49r
    p = pow7(l49r)
    oop = 1.0 / p
    root = z.root

    if root == 0xA:
        vx = 1; vy = 0
        if z.sub_hex > 1:
            vx += z.sub_hex - 2 - 2 * oop / 7
            vy += z.sub_hex - 2 + 1 * oop / 7
    elif root == 0xB:
        vx = 4; vy = 6
        if z.sub_hex > 1:
            vx -= z.sub_hex - 2 - 2 * oop / 7
            vy -= z.sub_hex - 2 + 1 * oop / 7
    else:
        cx = root >> 1
        cy = cx + (root & 1)
        vx = cx + z.col * oop
        vy = cy + z.row * oop

    if z.sub_hex and root < 10:
        sh = z.sub_hex
        south = bool(root & 1)
        if z.row == 0 and z.col == 0 and south and sh >= 4:
            sh += 1
        elif sh >= 2:
            parent_z = _Z(z.l49r, z.root, z.row, z.col, 0)
            if zone_is_edge_hex(parent_z):
                sh = (sh + (-1 if south else 3)) % 6 + 2

        oop_div7 = oop / 7
        ddx = 0.0; ddy = 0.0
        if sh == 2:
            ddx = -1 * oop_div7; ddy = -3 * oop_div7
        elif sh == 3:
            ddx = -3 * oop_div7; ddy = -2 * oop_div7
        elif sh == 4:
            ddx = -2 * oop_div7; ddy = 1 * oop_div7
        elif sh == 5:
            ddx = 1 * oop_div7; ddy = 3 * oop_div7
        elif sh == 6:
            ddx = 3 * oop_div7; ddy = 2 * oop_div7
        elif sh == 7:
            ddx = 2 * oop_div7; ddy = -1 * oop_div7
        if sh != 1:
            vx, vy = move5x6_vertex(vx, vy, ddx, ddy)

    if abs(vy) < 1e-6:
        vy = 0
    elif abs(vx - 5) < 1e-6:
        vx = 5
    if vx > 5 - 1e-6 or vy > 6 + 1e-6:
        vx -= 5; vy -= 5
    if vx < -1e-6:
        vx += 5; vy += 5

    return [vx, vy]


def contains_point(z: _Z, vx, vy):
    """Inside/boundary test for the zone's hexagon/pentagon polygon
    (``containsPoint``, RI7H.ec). The polygon is built from the centroid + vertex
    offset table; the seam-wrap fix-ups bring a point sitting across an
    interruption into the polygon's frame before the half-plane checks."""
    c = zone_centroid(z)
    level = zone_level(z)
    p = pow7(math.floor(level / 2))
    oonp = 1.0 / (7 * p)
    is_odd = level & 1
    hex_verts = ODD_HEX_VERTS if is_odd else EVEN_HEX_VERTS
    n_pts = zone_npoints(z)

    poly = []
    for i in range(n_pts):
        # The else-arm is unreachable: zone_npoints returns 5 or 6 and both
        # vertex tables hold 6 entries. Left explicit rather than simplified to
        # hex_verts[i] because it mirrors the eC's own wrap-tolerant indexing.
        vi = i if i < len(hex_verts) else 0
        px = c[0] + hex_verts[vi][0] * oonp
        py = c[1] + hex_verts[vi][1] * oonp
        poly.append([px, py])

    tvx, tvy = vx, vy
    if abs(tvx - c[0]) > 3 or abs(tvy - c[1]) > 3:
        if tvx > c[0] + 3:
            tvx -= 5; tvy -= 5
        elif tvx < c[0] - 3:
            tvx += 5; tvy += 5

    for i in range(len(poly)):
        j = (i + 1) % len(poly)
        ax, ay = poly[i][0], poly[i][1]
        bx, by = poly[j][0], poly[j][1]
        if abs(ax - tvx) > 3:
            if ax > 3 and ay > 3:
                ax -= 5; ay -= 5
            else:
                ax += 5; ay += 5
        if abs(bx - tvx) > 3:
            if bx > 3 and by > 3:
                bx -= 5; by -= 5
            else:
                bx += 5; by += 5
        dx = bx - ax; dy = by - ay
        side = dy * (tvx - ax) - dx * (tvy - ay)
        if side < -1e-11:
            return False
    return True


def add_non_polar_base_vertices(cx, cy, n_points, v):
    """Trace the cell boundary for a non-polar hex/pentagon, rotating the offset
    around the centroid and re-expressing it across each rhombus interruption it
    crosses (RI7H.ec ``getVertices`` non-polar path; spec §11). ``v`` is the scaled
    vertex offset table; returns the boundary points in 5x6 space.

    The crossing branch ports the eC ``addNonPolarBaseVertices`` faithfully,
    INCLUDING the four i1/i2 wrap-reconciliation blocks (``RI7H.ec:1991-2004``) that
    the igeo7-py oracle dropped — see the inline comment. Without them, the far-side
    vertices of interruption-spanning cells land out of the 5x6 band (→ degenerate
    (0,0) / wrong lat-lon on inverse). This is the one place py4dggs deliberately
    diverges from igeo7-py to match DGGAL/pydggal (the module docstring explains)."""
    start = 0
    for i in range(6):
        tx = cx + v[i][0]; ty = cy + v[i][1]
        itx = math.floor(tx + 1e-11)
        if not (ty - itx > 2 or ty < itx):
            start = i
            break

    point = [cx + v[start][0], cy + v[start][1]]
    prev = (start + 5) % 6
    direction = [point[0] - (cx + v[prev][0]), point[1] - (cy + v[prev][1])]

    vertices = []
    for i in range(start, start + n_points):
        direction = rotate5x6_offset(direction[0], direction[1], False)  # 60 CCW
        n = [point[0] + direction[0], point[1] + direction[1]]

        p = [point[0], point[1]]
        if p[0] > 5 and p[1] > 5:
            p[0] -= 5; p[1] -= 5
        if p[0] < 0 or p[1] < 0:
            p[0] += 5; p[1] += 5

        cross = crosses5x6_interruption_v2(p[0], p[1], direction[0], direction[1])
        if cross:
            # Reconcile the crossing in/out points (i1 = i_src, i2 = i_dst) into a
            # consistent frame before using them: first shift both into the
            # unwrapped ``point`` frame, then bring i2 into i1's frame. Dropping
            # these four wrap corrections (eC RI7H.ec:1991-2004) is the igeo7-py
            # port bug that gave interruption-spanning cells wrong / degenerate
            # (0,0) far-side vertices — the remainder direction and the new point
            # ``i2 + d`` were computed in the wrong (off-by-5) frame.
            i1 = list(cross["i_src"])
            i2 = list(cross["i_dst"])
            north = cross["in_north"]
            if point[0] - p[0] > 4:
                i1[0] += 5; i1[1] += 5; i2[0] += 5; i2[1] += 5
            if p[0] - point[0] > 4:
                i1[0] -= 5; i1[1] -= 5; i2[0] -= 5; i2[1] -= 5
            if i2[1] - i1[1] > 4:
                i2[0] -= 5; i2[1] -= 5
            if i1[1] - i2[1] > 4:
                i2[0] += 5; i2[1] += 5

            crossing_left = (i2[0] < i1[0]) if north else (i2[0] > i1[0])
            rdx, rdy = rotate5x6_offset(direction[0] - (i1[0] - point[0]),
                                        direction[1] - (i1[1] - point[1]), not crossing_left)
            n = [i2[0] + rdx, i2[1] + rdy]
            direction = rotate5x6_offset(direction[0], direction[1], not crossing_left)
            vertices.append(point)
            point = n
        else:
            # Inert by construction (the loop header already bounds i), and
            # FAITHFUL: the eC carries the identical dead guard --
            # `else if(i < start + nPoints)` at RI7H.ec:2014 against the loop
            # `for(i = start + 0; i < start + nPoints; i++)` at :1973. Kept per
            # the byte-faithfulness convention; do not "simplify" it away.
            if i < start + n_points:
                vertices.append(point)
            point = n
    return vertices


# --------------------------------------------------------------------------- #
# Hierarchy: candidate parents, children, parent walk (RI7H.ec; spec §9)
# --------------------------------------------------------------------------- #
def calc_candidate_parent(l49r, root, row, col, add_col, add_row):
    """Map a perturbed (row,col) back to its containing even-level zone, wrapping
    across rhombus seams and poles (``calcCandidateParent``, RI7H.ec). Returns
    ``None`` if the perturbed cell falls outside any rhombus."""
    p = pow7(l49r)
    r = root
    c = col + add_col
    rw = row + add_row
    south = bool(r & 1)

    if c == p and rw < p and not south:
        c = p - rw; rw = 0; r += 2
    elif rw == p and c < p and south:
        rw = p - c; c = 0; r += 2
    else:
        if rw < 0 and c < 0:
            rw += p; c += p; r -= 2
        elif rw < 0:
            rw += p; r -= 1
        elif c < 0:
            c += p; r -= 1
        elif c >= p and rw >= p:
            rw -= p; c -= p; r += 2
        elif rw >= p:
            rw -= p; r += 1
        elif c >= p:
            c -= p; r += 1
    if r < 0:
        r += 10
    elif r > 9:
        r -= 10

    south = bool(r & 1)
    if not south and rw == 0 and c == p:
        return _Z(l49r, 0xA, 0, 0, 0)
    if south and rw == p and c == 0:
        return _Z(l49r, 0xB, 0, 0, 0)

    if rw < 0 or rw >= p or c < 0 or c >= p:
        return None
    return _Z(l49r, r, rw, c, 0)


def get_odd_level_centroid_child_root_row_col(z: _Z):
    """Root/row/col of the even-level centroid child of an odd-level zone — the
    aperture-49 child sitting on the parent's centroid (``getOddLevelCentroid
    ChildRootRowCol``, RI7H.ec). Handles the poles, the south-corner and edge-hex
    sub-hex rotations, then the cp-scaled seam wrapping."""
    p = pow7(z.l49r)
    cp = p * 7
    root = z.root

    if root == 0xA:
        if z.sub_hex == 1:
            return {"root": 0xA, "row": 0, "col": 0, "cp": cp}
        return {"root": 2 * (z.sub_hex - 2), "row": 1, "col": cp - 2, "cp": cp}
    if root == 0xB:
        if z.sub_hex == 1:
            return {"root": 0xB, "row": 0, "col": 0, "cp": cp}
        return {"root": 9 - 2 * (z.sub_hex - 2), "row": cp - 1, "col": 2, "cp": cp}

    sh = z.sub_hex
    c_rhombus = root
    south = bool(c_rhombus & 1)
    row = 7 * z.row
    col = 7 * z.col

    if z.row == 0 and z.col == 0 and south and sh >= 4:
        sh += 1
    elif sh >= 2:
        parent_z = _Z(z.l49r, z.root, z.row, z.col, 0)
        if zone_is_edge_hex(parent_z):
            sh = (sh + (-1 if south else 3)) % 6 + 2

    if sh == 2:
        row -= 3; col -= 1
    elif sh == 3:
        row -= 2; col -= 3
    elif sh == 4:
        row += 1; col -= 2
    elif sh == 5:
        row += 3; col += 1
    elif sh == 6:
        row += 2; col += 3
    elif sh == 7:
        row -= 1; col += 2

    if col == cp and row < cp and not south:
        col = cp - row; row = 0; c_rhombus += 2
    elif row == cp and col < cp and south:
        row = cp - col; col = 0; c_rhombus += 2
    elif col > 0 and col < cp and row < 0 and not south:
        ncol = cp + row; nrow = cp - col; col = ncol; row += nrow; c_rhombus -= 2
    elif row > 0 and row < cp and col < 0 and south:
        nrow = cp + col; ncol = cp - row; row = nrow; col += ncol; c_rhombus -= 2

    if row < 0 and col < 0:
        row += cp; col += cp; c_rhombus -= 2
    elif row < 0:
        row += cp; c_rhombus -= 1
    elif col < 0:
        col += cp; c_rhombus -= 1
    elif col >= cp and row >= cp:
        row -= cp; col -= cp; c_rhombus += 2
    elif row >= cp:
        row -= cp; c_rhombus += 1
    elif col >= cp:
        col -= cp; c_rhombus += 1

    if c_rhombus < 0:
        c_rhombus += 10
    elif c_rhombus > 9:
        c_rhombus -= 10

    if 0 <= row < cp and 0 <= col < cp:
        return {"root": c_rhombus, "row": row, "col": col, "cp": cp}
    return None


def get_primary_children(z: _Z):
    """The 6 or 7 primary children of a zone (``getPrimaryChildren``, spec §9).
    Even -> odd just enumerates sub-hexes 1..7 (6 for the centre pentagon). Odd ->
    even places the centroid child plus the ring of six offset children, wrapping
    each across rhombus seams (with the polar fan-out and the edge-hex rotation)."""
    if z is None:
        return []
    l49r, root, z_row, z_col, sub_hex = z.l49r, z.root, z.row, z.col, z.sub_hex
    if l49r > 9 or (l49r == 9 and sub_hex):
        return []

    if sub_hex == 0:
        children = [_Z(l49r, root, z_row, z_col, sh) for sh in range(1, 7)]
        if z_row != 0 or z_col != 0:
            children.append(_Z(l49r, root, z_row, z_col, 7))
        return children

    cc = get_odd_level_centroid_child_root_row_col(z)
    if not cc:
        return []
    c_root, crow, ccol, cp = cc["root"], cc["row"], cc["col"], cc["cp"]
    c_child = _Z(l49r + 1, c_root, crow, ccol, 0)
    children = [c_child]

    if c_root == 0xA:
        for i in range(5):
            children.append(_Z(l49r + 1, 2 * i, 0, cp - 1, 0))
    elif c_root == 0xB:
        for i in range(5):
            children.append(_Z(l49r + 1, 9 - 2 * i, cp - 1, 0, 0))
    else:
        n_points = zone_npoints(z)
        south = bool(c_root & 1)
        c_offsets = [[-1, 0], [-1, -1], [0, -1], [1, 0], [1, 1], [0, 1]]

        edge_hex_fix = False
        if n_points == 6 and sub_hex == 1:
            parent_z = _Z(l49r, root, z_row, z_col, 0)
            if zone_is_edge_hex(parent_z):
                edge_hex_fix = True

        for i in range(6):
            ii = (i + (1 if south else 5)) % 6 if edge_hex_fix else i

            if n_points == 5:
                if south and i == 2:
                    continue
                if not south and i == 0:
                    continue

            row = crow + c_offsets[ii][0]
            col = ccol + c_offsets[ii][1]
            c_rhombus = c_root

            if col == cp and row < cp and not south:
                col = cp - row; row = 0; c_rhombus += 2
            elif row == cp and col < cp and south:
                row = cp - col; col = 0; c_rhombus += 2
            elif col > 0 and col < cp and row < 0 and not south:
                ncol = cp + row; nrow = cp - col; col = ncol; row += nrow; c_rhombus -= 2
            elif row > 0 and row < cp and col < 0 and south:
                nrow = cp + col; ncol = cp - row; row = nrow; col += ncol; c_rhombus -= 2
            else:
                if row < 0 and col < 0:
                    row += cp; col += cp; c_rhombus -= 2
                elif row < 0:
                    row += cp; c_rhombus -= 1
                elif col < 0:
                    col += cp; c_rhombus -= 1
                elif col >= cp and row >= cp:
                    row -= cp; col -= cp; c_rhombus += 2
                elif row >= cp:
                    row -= cp; c_rhombus += 1
                elif col >= cp:
                    col -= cp; c_rhombus += 1

            if c_rhombus < 0:
                c_rhombus += 10
            elif c_rhombus > 9:
                c_rhombus -= 10

            if 0 <= row < cp and 0 <= col < cp:
                children.append(_Z(l49r + 1, c_rhombus, row, col, 0))
    return children


def get_centroid_child(z: _Z):
    """The single child sharing the parent's centroid (sub-hex 1 for an even zone;
    the even-level centroid child for an odd zone). Used by the rotation offset."""
    if z is None:
        return None
    if z.sub_hex == 0:
        return _Z(z.l49r, z.root, z.row, z.col, 1)
    result = get_odd_level_centroid_child_root_row_col(z)
    if not result:
        return None
    return _Z(z.l49r + 1, result["root"], result["row"], result["col"], 0)


def from_even_level_primary_child(child: _Z):
    """Find the odd-level parent of an even-level zone by inverting the centroid
    geometry: quantize the child's centroid back to a rhombus cell, enumerate
    candidate parents, and return the one whose grandchildren include ``child``
    (RI7H.ec even-level parent search). Returns ``None`` if no parent matches."""
    l49r = child.l49r - 1
    if child.sub_hex or l49r < 0:
        return None

    cx, cy = zone_centroid(child)
    x, y = cx, cy
    p = pow7(l49r)
    oop = 1.0 / p

    if abs(x - y - 1) < 1e-10:
        pass
    elif abs(y - x - 2) < 1e-10:
        pass
    elif y < -1e-11 and x > -1e-11:
        x -= y; y = 0
    elif math.floor(x + 1e-11) > math.floor(y + 1e-11):
        iy = min(5, math.floor(y + 1e-11))
        x += (iy + 1 - y); y = iy + 1
    elif math.floor(y + 1e-11) - math.floor(x + 1e-11) > 1:
        ix = min(4, math.floor(x + 1e-11))
        y += (ix + 1 - x); x = ix + 1
    elif x < -1e-11 or y < -1e-11:
        x += 5; y += 5

    if x > 5 - 1e-11 and y > 5 - 1e-11 and x + y > 10 - oop - 1e-11:
        x -= 5; y -= 5

    ix = min(4, math.floor(x + 1e-11))
    iy = min(5, math.floor(y + 1e-11))
    root = ix + iy
    fx = x - ix; fy = y - iy
    col = max(0, _jsround(fx * p))
    row = max(0, _jsround(fy * p))
    south = y - x - 1e-11 > 1
    north = x - y - 1e-11 > 0
    north_pole = north and abs(x - y - 1) < 1e-11
    south_pole = south and abs(y - x - 2) < 1e-11

    if north_pole:
        return _Z(l49r, 0xA, 0, 0, 1)
    if south_pole:
        return _Z(l49r, 0xB, 0, 0, 1)

    candidate_parents = [None] * 7
    if north and row == 0 and col == p:
        candidate_parents[0] = _Z(l49r, 0xA, 0, 0, 0)
    elif south and row == p and col == 0:
        candidate_parents[0] = _Z(l49r, 0xB, 0, 0, 0)
    else:
        candidate_parents[0] = calc_candidate_parent(l49r, root, row, col, 0, 0)
    candidate_parents[1] = calc_candidate_parent(l49r, root, row, col, 0, -1)
    candidate_parents[2] = calc_candidate_parent(l49r, root, row, col, 0, 1)
    candidate_parents[3] = calc_candidate_parent(l49r, root, row, col, 1, 0)
    candidate_parents[4] = calc_candidate_parent(l49r, root, row, col, -1, 0)
    candidate_parents[5] = calc_candidate_parent(l49r, root, row, col, -1, -1)
    candidate_parents[6] = calc_candidate_parent(l49r, root, row, col, 1, 1)

    for cand in candidate_parents:
        if cand is None:
            continue
        for ch in get_primary_children(cand):
            for gch in get_primary_children(ch):
                if zone_eq(gch, child):
                    return ch
    return None


def get_parent0(z: _Z):
    """The immediate parent of a zone: drop the sub-hex for an odd zone, else find
    the odd parent of an even zone via :func:`from_even_level_primary_child`.
    Returns ``None`` at the base level."""
    if z.sub_hex > 0:
        return _Z(z.l49r, z.root, z.row, z.col, 0)
    level = zone_level(z)
    if level == 0:
        return None
    return from_even_level_primary_child(z)


def compute_parents(zone: _Z):
    """The full parent chain from ``zone`` up to (but excluding) the base zone."""
    parents = []
    z = zone
    while zone_level(z) > 0:
        parent = get_parent0(z)
        if parent is None:
            break
        parents.append(parent)
        z = parent
    return parents


# --------------------------------------------------------------------------- #
# fromCentroid: quantize a 5x6 point to an I7H zone (RI7H.ec; spec §6)
# --------------------------------------------------------------------------- #
def from_centroid(level, cx, cy):
    """Quantize a 5x6 planar point to the I7H zone whose cell contains it, at
    ``level`` (``I7HZone::fromCentroid``, spec §6).

    Even levels (§6 even) snap to the rhombus lattice, then nudge the (row,col) by
    one cell using the sub-cell offset ``(dx,dy)`` and the four hex-edge half-plane
    tests. Odd levels (§6 odd, brute-force) enumerate the candidate even parents'
    children and pick the one whose polygon contains the point, falling back to the
    finer even level's parent and finally the rhombus-corner centre."""
    l49r = math.floor(level / 2)
    p = pow7(l49r)
    oop = 1.0 / p

    x, y = canonicalize5x6(cx, cy)

    if x > 5 - 1e-11 and y > 5 - 1e-11 and x + y > 10 - oop - 1e-11:
        x -= 5; y -= 5

    ix = min(4, math.floor(x + 1e-11))
    iy = min(5, math.floor(y + 1e-11))
    root = ix + iy
    fx = x - ix; fy = y - iy
    col = max(0, _jsround(fx * p))
    row = max(0, _jsround(fy * p))
    dx = fx * p + 0.5 - col
    dy = fy * p + 0.5 - row
    south_rhombus = bool(root & 1)
    south = y - x - 1e-11 > 1
    north = x - y - 1e-11 > 0
    north_pole = north and abs(x - y - 1) < 1e-11
    south_pole = south and abs(y - x - 2) < 1e-11

    if level & 1:
        if north_pole:
            return _Z(l49r, 0xA, 0, 0, 1)
        if south_pole:
            return _Z(l49r, 0xB, 0, 0, 1)

        candidate_parents = [None] * 7
        if north and row == 0 and col == p:
            candidate_parents[0] = _Z(l49r, 0xA, 0, 0, 0)
        elif south and row == p and col == 0:
            candidate_parents[0] = _Z(l49r, 0xB, 0, 0, 0)
        else:
            candidate_parents[0] = calc_candidate_parent(l49r, root, row, col, 0, 0)
        candidate_parents[1] = calc_candidate_parent(l49r, root, row, col, 0, -1)
        candidate_parents[2] = calc_candidate_parent(l49r, root, row, col, 0, 1)
        candidate_parents[3] = calc_candidate_parent(l49r, root, row, col, 1, 0)
        candidate_parents[4] = calc_candidate_parent(l49r, root, row, col, -1, 0)
        candidate_parents[5] = calc_candidate_parent(l49r, root, row, col, -1, -1)
        candidate_parents[6] = calc_candidate_parent(l49r, root, row, col, 1, 1)

        for cand in candidate_parents:
            if cand is None:
                continue
            for ch in get_primary_children(cand):
                if contains_point(ch, x, y):
                    return ch

        even_child = from_centroid(level + 1, cx, cy)
        if even_child is not None:
            odd_parent = from_even_level_primary_child(even_child)
            if odd_parent is not None:
                return odd_parent
        return _Z(l49r, root, 0, 0, 1)

    if north_pole:
        return _Z(l49r, 0xA, 0, 0, 0)
    if south_pole:
        return _Z(l49r, 0xB, 0, 0, 0)

    if dx > 1 - dy:
        if dx > dy:
            if point_line_side(dx, dy, 1.0, 0.5, 5 / 6, 1 / 6) < 0:
                col += 1
        else:
            if point_line_side(dx, dy, 1 / 6, 5 / 6, 0.5, 1.0) < 0:
                row += 1
    else:
        if dx > dy:
            if point_line_side(dx, dy, 5 / 6, 1 / 6, 0.5, 0.0) < 0:
                row -= 1
        else:
            if point_line_side(dx, dy, 0.0, 0.5, 1 / 6, 5 / 6) < 0:
                col -= 1

    if north and col == p and row == 0:
        return _Z(l49r, 0xA, 0, 0, 0)
    if south and col == 0 and row == p:
        return _Z(l49r, 0xB, 0, 0, 0)

    if col == p and row < p and not south_rhombus:
        col = p - row; row = 0; root += 2
    elif row == p and col < p and south_rhombus:
        row = p - col; col = 0; root += 2
    else:
        if row < 0 and col < 0:
            row += p; col += p; root -= 2
        elif row < 0:
            row += p; root -= 1
        elif col < 0:
            col += p; root -= 1
        elif col >= p and row >= p:
            row -= p; col -= p; root += 2
        elif row >= p:
            row -= p; root += 1
        elif col >= p:
            col -= p; root += 1
    if root < 0:
        root += 10
    elif root > 9:
        root -= 10

    if row < 0 or row >= p or col < 0 or col >= p:
        # Audit A6: eC fromCentroid returns nullZone here (RI7H.ec:1680); this
        # port returns a corner cell instead. Verified UNREACHABLE via the
        # forward/quantize path (0 hits over 200k random points across all
        # resolutions), so the divergence is inert -- kept as-is rather than
        # perturb a hot, otherwise byte-faithful path. (It is also why the
        # even-level fromCentroid never returns None, making the `is not None`
        # guard on its odd-level caller a defensive no-op.)
        return _Z(l49r, root, 0, 0, 0)

    return _Z(l49r, root, row, col, 0)


# --------------------------------------------------------------------------- #
# I7H <-> Z7 conversion (RI7H_Z7.ec)
# --------------------------------------------------------------------------- #
def get_child_position(parent: _Z, zone: _Z):
    """Index of ``zone`` among ``parent``'s primary children (0-based). Odd levels
    read it straight off the sub-hex; even levels scan the children list."""
    level = zone_level(zone)
    if level & 1:
        return zone.sub_hex - 1
    children = get_primary_children(parent)
    for i, ch in enumerate(children):
        if zone_eq(ch, zone):
            return i
    return 0


def adjust_z7_pentagon_child_position(i, level, p_root):
    """Re-index a child position around a pentagon base cell so the Z7 digit
    ordering is consistent (``adjustZ7PentagonChildPosition``, RI7H_Z7.ec)."""
    if not i:
        return 0
    south_p_rhombus = bool(p_root & 1)
    odd_level = bool(level & 1)
    if p_root == 10:
        i = ((i + 1) % 5) + 1
    elif p_root == 11:
        i = ((i + (3 if odd_level else 4)) % 5) + 1
    elif not odd_level and not south_p_rhombus:
        i = ((i + 5) % 5) + 1
    if south_p_rhombus and i >= 3:
        i += 1
    return i


def deadjust_z7_pentagon_child_position(i, level, p_root):
    """Inverse of :func:`adjust_z7_pentagon_child_position`
    (``deadjustZ7PentagonChildPosition``, RI7H_Z7.ec) — used when reconstructing
    the I7H child position from a Z7 digit."""
    if not i:
        return 0
    south_p_rhombus = bool(p_root & 1)
    odd_level = bool(level & 1)
    if south_p_rhombus and i >= 4:
        i -= 1
    if p_root == 10:
        i = ((i - 1 + 3) % 5) + 1
    elif p_root == 11:
        i = ((i - 1 + (1 if odd_level else 5)) % 5) + 1
    elif not odd_level and not south_p_rhombus:
        i = ((i - 1 + 4) % 5) + 1
    return i


def get_level_rotation_offset(l, prev_i, zone, parent, grand_parent):
    """Per-level rotation correction accumulated when walking parent->child, so
    the Z7 digit ring stays aligned through pentagons and edge hexagons
    (``getLevelRotationOffset``, RI7H_Z7.ec). The cascade of cases mirrors the eC
    exactly, including its child-position sentinel: the eC opens with
    ``if(i == -1) i = getChildPosition(parent, zone);`` (RI7H_Z7.ec:200-201), so
    ``prev_i`` is an INPUT -- callers that already know the position pass it
    (:377 passes ``prevI``, :330 passes ``prevCIX``) and only ``-1`` asks the
    function to look it up (:276). This port used to recompute unconditionally,
    which dropped the eC's branch and left the parameter inert. Behaviour is
    unchanged: the supplied value is a cache of the same lookup (measured over
    36,163 calls, 0 divergences), so this restores the structure, not the result."""
    if parent is None:
        return 0
    i = prev_i
    if i == -1:
        i = get_child_position(parent, zone)
    if not i:
        return 0

    offset = 0
    p_root = parent.root
    pn_points = zone_npoints(parent)
    odd_level = bool(l & 1)
    south_p_rhombus = bool(p_root & 1)
    is_edge_hex = (not odd_level) and zone_is_edge_hex(zone)
    p_edge_hex = odd_level and zone_is_edge_hex(parent)
    gp_edge_hex = (not odd_level) and grand_parent is not None and zone_is_edge_hex(grand_parent)

    if pn_points == 5:
        i = adjust_z7_pentagon_child_position(i, l, p_root)

    if p_root >= 10:
        if pn_points == 5:
            offset += i + (
                (0 if south_p_rhombus else 3) if odd_level
                else (5 if south_p_rhombus else 2)
            )
        elif is_edge_hex and (not south_p_rhombus or not zone_eq(zone, get_centroid_child(parent))):
            offset += 5

    if south_p_rhombus and is_edge_hex:
        offset += 1

    if p_edge_hex:
        if not south_p_rhombus and i >= 4:
            offset += 1
        elif south_p_rhombus and (i == 0 or (3 <= i <= 5)):
            offset += 5
    elif gp_edge_hex:
        pc = get_primary_children(grand_parent)
        c = get_primary_children(parent)
        if south_p_rhombus:
            if (len(pc) > 1 and zone_eq(pc[1], parent) and
                    len(c) > 2 and len(c) > 5 and c[2].root != c[5].root and
                    (i == 4 or i == 5)):
                offset += 5
        else:
            if len(pc) > 4 and zone_eq(pc[4], parent) and (i == 1 or i == 2):
                offset += 5

        if zone_eq(parent, get_centroid_child(grand_parent)):
            if south_p_rhombus:
                if i > 2:
                    offset += 5
            else:
                if i == 5 or i == 6:
                    offset += 1

        if south_p_rhombus and is_edge_hex:
            if i == 4 or i == 5:
                offset += 5
    return offset


def from_7h(zone: _Z):
    """Convert an I7H zone to its Z7 address: base pentagon + a 20-level packed
    ancestry of 3-bit digits (``Z7Zone::from7H``, RI7H_Z7.ec). Walks the parent
    chain top-down, applying the rotation offset and pentagon re-indexing at each
    level, then maps each I7H child position through ``C_MAP``. Returns ``None``
    above the representable level."""
    if zone is None:
        return None
    level = zone_level(zone)
    if level > 20:
        return None

    parents = compute_parents(zone)
    ancestry = 0
    offset = 0
    prev_i = 0

    for l in range(1, level + 1):
        p_index = level - l
        z = zone if l == level else parents[p_index - 1]
        parent = parents[p_index]
        grand_parent = parents[p_index + 1] if l > 1 else None
        great_grand_parent = parents[p_index + 2] if l > 2 else None

        i = get_child_position(parent, z)
        offset = (offset + get_level_rotation_offset(l - 1, prev_i, parent, grand_parent, great_grand_parent)) % 6
        prev_i = i

        if i:
            if zone_npoints(parent) == 5:
                i = adjust_z7_pentagon_child_position(i, l, parent.root)
            i = ((i - 1) + offset) % 6 + 1

        shift = (19 - (l - 1)) * 3
        ancestry |= C_MAP[i] << shift

    for l in range(level + 1, 21):
        shift = (19 - (l - 1)) * 3
        ancestry |= 7 << shift

    top_parent = zone if level == 0 else parents[-1]
    root_pentagon = ROOT_MAP[top_parent.root]

    return {"root_pentagon": root_pentagon, "ancestry": ancestry, "level": level}


def _zone_from_steps(base_cell, directions):
    """Reconstruct the I7H zone from a Z7 base cell + direction digits
    (``Z7Zone::to7H``, RI7H_Z7.ec) — the inverse of :func:`from_7h`. Walks the
    digits top-down, undoing the rotation offset and pentagon re-indexing, stepping
    even->odd (sub-hex) and odd->even (primary children) alternately."""
    z7_root = INV_ROOT_MAP[base_cell]
    zone = _Z(0, z7_root, 0, 0, 0)
    offset = 0
    prev_cix = 0

    for l in range(len(directions)):
        b = directions[l]
        if b == 7:
            break
        cix = INV_C_MAP[b]
        n_points = zone_npoints(zone)
        parent = get_parent0(zone) if l > 0 else None
        grand_parent = get_parent0(parent) if (l > 1 and parent) else None

        if cix or l < 19:
            offset = (offset + get_level_rotation_offset(l, prev_cix, zone, parent, grand_parent)) % 6
        if cix:
            cix = cix - 1 - offset
            if cix < 0:
                cix += 6
            cix += 1
            if n_points == 5:
                cix = deadjust_z7_pentagon_child_position(cix, l + 1, zone.root)
        prev_cix = cix

        if not (l & 1):
            zone = _Z(zone.l49r, zone.root, zone.row, zone.col, 1 + cix)
        else:
            children = get_primary_children(zone)
            if cix < len(children):
                zone = children[cix]
            else:
                break
    return zone


# --------------------------------------------------------------------------- #
# Public planar geometry: centroid & vertices in 5x6 space (RI7H.ec)
# --------------------------------------------------------------------------- #
def _centroid5x6(zone):
    """The zone centroid in 5x6 space, folded into the canonical domain
    (``getCentroid`` + the final wrap). Returns ``(cx, cy)``."""
    cx, cy = zone_centroid(zone)
    if cy > 6 + 1e-9 or cx > 5 + 1e-9:
        cx -= 5; cy -= 5
    elif cx < 0:
        cx += 5; cy += 5
    return cx, cy


def _verts5x6(zone, cx, cy, level):
    """The zone's 5/6 vertices in 5x6 space around centroid ``(cx,cy)``
    (``I7HZone::getVertices``). Polar base cells (roots 0xA/0xB) fan five vertices
    along the pole diagonal (or trace the odd-level hex); all other cells trace the
    boundary via :func:`add_non_polar_base_vertices`. A final seam wrap brings each
    vertex into the canonical domain."""
    p = pow7(math.floor(level / 2))
    oonp = 1.0 / (7 * p)
    is_odd = level & 1
    n_pts = zone_npoints(zone)
    verts5x6 = None

    if zone.root == 0xA or zone.root == 0xB:
        is_pole = zone.sub_hex <= 1
        if zone.root == 0xA:
            bx = by = None
            if (not is_odd) or is_pole:
                A = ODD_HEX_C if is_odd else EVEN_HEX_A
                bx = 1 - oonp * A
                by = 0 + oonp * (ODD_HEX_A if is_odd else EVEN_HEX_A)
            else:
                v = [[ODD_HEX_VERTS[i][0] * oonp, ODD_HEX_VERTS[i][1] * oonp] for i in range(6)]
                verts5x6 = add_non_polar_base_vertices(cx, cy, n_pts, v)
            if verts5x6 is None:
                verts5x6 = [[bx + i, by + i] for i in range(5)]
        else:
            bx = by = None
            if (not is_odd) or is_pole:
                A = ODD_HEX_C if is_odd else EVEN_HEX_A
                bx = 4 + oonp * A
                by = 6 - oonp * (ODD_HEX_A if is_odd else EVEN_HEX_A)
            else:
                v = [[ODD_HEX_VERTS[i][0] * oonp, ODD_HEX_VERTS[i][1] * oonp] for i in range(6)]
                verts5x6 = add_non_polar_base_vertices(cx, cy, n_pts, v)
            if verts5x6 is None:
                verts5x6 = [[bx - i, by - i] for i in range(5)]
    else:
        hex_verts = ODD_HEX_VERTS if is_odd else EVEN_HEX_VERTS
        v = [[hex_verts[i][0] * oonp, hex_verts[i][1] * oonp] for i in range(6)]
        verts5x6 = add_non_polar_base_vertices(cx, cy, n_pts, v)

    out = []
    for vx, vy in verts5x6:
        if vx > 5 and vy > 5:
            vx -= 5; vy -= 5
        elif vx < 0 and vy < 1:
            vx += 5; vy += 5
        out.append([vx, vy])
    return out


# --------------------------------------------------------------------------- #
# The Topology facade
# --------------------------------------------------------------------------- #
class HexAperture7Topology:
    """Aperture-7 hexagonal topology on the 5x6 ISEA grid — a ``Topology`` impl.

    Maps a projection's planar point to a Z7 ``(base cell, direction digits)``
    address and back to planar cell geometry. ``geom`` is accepted for protocol
    parity (the 5x6 layout is fixed across ISEA-family configs, so this topology
    does not read it); the *Grid* lifts the planar centroid/vertices to lat/lon via
    the Projection. Output is bit-identical to the ``igeo7`` oracle by construction.
    """

    aperture = 7

    def quantize(self, geom: Any, p: PlanarPoint, res: int) -> tuple[int, list[int]]:
        """Quantize planar point ``p`` to a Z7 cell at resolution ``res``, returning
        ``(base_cell, directions)``. Like the oracle, ``p.face`` is unused — the
        cell is determined from ``(p.x, p.y)`` alone via ``fromCentroid`` then the
        I7H->Z7 conversion; the packed ancestry is unpacked into the first ``res``
        direction digits (stopping at the level-7 terminator)."""
        i7h = from_centroid(res, p.x, p.y)
        z7 = from_7h(i7h)
        if not z7:
            # from_7h returns None only above the representable level (level > 20,
            # :1077-1078), and Grid.zone_from_geo rejects res > max_resolution = 19
            # first, so nothing in this library reaches here. Unlike hex_a3 -- whose
            # quantize returns (NULL_ZONE, []) -- the Z7 (base, digits) protocol has
            # no null sentinel to return: base 0 with no digits IS a real zone (base
            # cell 0 at resolution 0), so this arm would silently claim a point
            # quantized to the null island if it ever became reachable. Adding a Z7
            # nullZone sentinel is the proper fix and a protocol change; until then
            # the guard is Grid.zone_from_geo's resolution bound, one layer up.
            return 0, []

        base_cell = z7["root_pentagon"]
        ancestry = z7["ancestry"]
        directions: list[int] = []
        for l in range(res):
            shift = (19 - l) * 3
            d = (ancestry >> shift) & 7
            if d == 7:
                break
            directions.append(d)
        return base_cell, directions

    def is_null_geometry(self, geom: Any, base: int, digits: list[int]) -> bool:
        """True when this address has NO geometry — DGGAL's ``nullZone`` outcome.

        Z7 level 20 is representable (the 64-bit packing has 20 direction-digit
        slots) but is not a valid DGGAL zone: the Z7 -> 7H conversion
        (``to7H``, RI7H_Z7.ec:348-353, "does not support level 20 zones")
        returns nullZone, so ``getZoneWGS84Centroid``/``Vertices`` degenerate to
        the zero point. Verified against pydggal for all three aperture-7 grids:
        a level-20 zone yields centroid (0,0) and six (0,0) vertices, for both
        hexagon- and pentagon-paths, while level 19 is normal (and a
        pentagon-path level-19 zone correctly yields 5 vertices).

        ``Grid`` consults this before calling ``planar_centroid``/
        ``planar_vertices``; without it this topology happily computes a
        plausible-looking coordinate for a zone DGGAL says does not exist."""
        return len(digits) >= _NULL_GEOMETRY_LEVEL

    def planar_centroid(self, geom: Any, base: int, digits: list[int]) -> PlanarPoint:
        """The cell's centroid as a planar ``PlanarPoint`` in 5x6 space. ``face`` is
        the ``-1`` sentinel ("derive from x,y"); the Grid runs the projection's
        inverse to get lat/lon."""
        cx, cy = _centroid5x6(_zone_from_steps(base, digits))
        return PlanarPoint(face=-1, x=cx, y=cy)

    def planar_vertices(self, geom: Any, base: int, digits: list[int]) -> list[PlanarPoint]:
        """The cell's boundary vertices (6 for a hexagon, 5 for a pentagon) as
        planar ``PlanarPoint``s in 5x6 space, each with the ``-1`` face sentinel.
        The resolution (Z7 level) is ``len(digits)``."""
        zone = _zone_from_steps(base, digits)
        cx, cy = _centroid5x6(zone)
        level = len(digits)
        verts = _verts5x6(zone, cx, cy, level)
        return [PlanarPoint(face=-1, x=vx, y=vy) for vx, vy in verts]
