"""icovertex — the shared icosahedral great-circle equal-area projection kernel.

This is the projection half of the ISEA7H / IVEA7H family: it maps a geographic
lat/lon onto a point ``(face, x, y)`` in DGGAL's oblique 5x6 grid (the planar
layout the hexagon topology then quantizes), and back. It is equal-area by
construction — Snyder's 1992 great-circle vector method, refined in DGGAL's
``icoVertexGreatCircle.ec``.

DGGAL's ``SliceAndDiceGreatCircleIcosahedralProjection`` is one kernel with a
``radialVertex`` enum (``isea`` / ``ivea`` / ``rtea``) selecting the per-variant
spherical-triangle edge angles and the A-vs-BC sub-triangle assignment. This
module mirrors that: :class:`_IcoVertexProjection` holds ALL the shared math and
is parameterized by :attr:`~_IcoVertexProjection.radial_vertex`; the thin variant
classes (:mod:`py4dggs.projections.isea`, :mod:`py4dggs.projections.ivea`) just bind it
to a variant. The variant-specific state is two things (per the eC property setter
``icoVertexGreatCircle.ec:118-158``): the triangle-edge constants
(``AB``/``AC``/``BC`` and their trig), computed by :func:`variant_consts`, and the
``b_is_a`` flag, which XORs ``(radial_vertex == "ivea")`` into the sub-triangle test.

Provenance (the *why* of each block is cited inline against these sources):
  - ``projections/ri5x6.ec``           — orientation, the 12 icosahedron vertices,
                                          the 20-face / 5x6 layout, ``find_face``.
  - ``projections/icoVertexGreatCircle.ec`` — the equal-area forward/inverse vector
                                          projection (``sqrtOneMinusDotOver2``, the
                                          midpoint-cross ``sphericalTriArea``, the
                                          barycentric sub-triangle decomposition,
                                          the ``radialVertex`` variant selection).
  - ``projections/authalic.ec``        — geodetic<->authalic latitude (Karney 2022
                                          Clenshaw series over the WGS84 ellipsoid).
  - ``projections/barycentric5x6.ec``  — a DIFFERENT projection class (Goldberg);
                                          the barycentric / spherical-area / slerp
                                          helpers this module ports
                                          (``cartesianToBary``, ``baryToCartesian``,
                                          ``sphericalTriArea``, ``slerpAngle``) actually
                                          live in ``ri5x6.ec`` — barycentric5x6.ec only
                                          *calls* them, so cite ``ri5x6.ec`` for them.

It is re-expressed from ``igeo7/isea.py`` (itself a bit-identical port of the eC).
Readability lives in the structure — an :class:`_IcoVertexProjection` class
implementing the ``Projection`` protocol, a frozen :class:`_Geometry` of precomputed
per-config state, named free-function helpers, and these docstrings. The *arithmetic*
is kept as close to the oracle as practical (same expressions, groupings, epsilons and
evaluation order) so the differential fuzz passes. A FEW spots deliberately substitute
an analytically-equal form (e.g. the authalic-triangle constant for a runtime
scalar-triple-product, or ``sin`` via a double-angle expansion) that differs only at
~1e-16; the acceptance gate is therefore a tight *tolerance*, not literal last-bit
identity — so don't "fix" a harmless reassociation on the strength of a bit-identity
that isn't actually being claimed.
"""
from __future__ import annotations

import math
from collections import namedtuple
from dataclasses import dataclass
from typing import Any

from py4dggs.types import GeoPoint, GridConfig, PlanarPoint

DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi
TWO_PI = 2.0 * math.pi
PHI = (1.0 + math.sqrt(5.0)) / 2.0  # golden ratio

# --- Authalic latitude conversion (Karney 2022) over WGS84 (authalic.ec) ---
WGS84_A = 6378137.0
WGS84_B = 6356752.314245179

Cxiphi = [  # geodetic -> authalic (21 coefficients, flat)
    -4/3,    -4/45,    88/315,       538/4725,     20824/467775,      -44732/2837835,
             34/45,     8/105,     -2482/14175,    -37192/467775,   -12467764/212837625,
                     -1532/2835,     -898/14175,     54968/467775,   100320856/1915538625,
                                     6007/14175,     24496/467775,    -5884124/70945875,
                                                    -23356/66825,     -839792/19348875,
                                                                    570284222/1915538625,
]
Cphixi = [  # authalic -> geodetic (21 coefficients, flat)
    4/3,  4/45,   -16/35,  -2582/14175,  60136/467775,    28112932/212837625,
         46/45,  152/945, -11966/14175, -21016/51975,   251310128/638512875,
                3044/2835,   3802/14175, -94388/66825,    -8797648/10945935,
                             6059/4725,  41072/93555, -1472637812/638512875,
                                        768272/467775,  -455935736/638512875,
                                                       4210684958/1915538625,
]

# --- Face vertex indices (DGGAL icoIndices[20][3], ri5x6.ec) ---
FACE_VERTICES = [
    [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1],
    [6, 2, 1], [7, 3, 2], [8, 4, 3], [9, 5, 4], [10, 1, 5],
    [2, 6, 7], [3, 7, 8], [4, 8, 9], [5, 9, 10], [1, 10, 6],
    [11, 7, 6], [11, 8, 7], [11, 9, 8], [11, 10, 9], [11, 6, 10],
]

# --- 5x6 planar layout coordinates for each face (ri5x6.ec) ---
FACE_5X6 = [
    [[1, 0], [0, 0], [1, 1]], [[2, 1], [1, 1], [2, 2]], [[3, 2], [2, 2], [3, 3]], [[4, 3], [3, 3], [4, 4]], [[5, 4], [4, 4], [5, 5]],
    [[0, 1], [1, 1], [0, 0]], [[1, 2], [2, 2], [1, 1]], [[2, 3], [3, 3], [2, 2]], [[3, 4], [4, 4], [3, 3]], [[4, 5], [5, 5], [4, 4]],
    [[1, 1], [0, 1], [1, 2]], [[2, 2], [1, 2], [2, 3]], [[3, 3], [2, 3], [3, 4]], [[4, 4], [3, 4], [4, 5]], [[5, 5], [4, 5], [5, 6]],
    [[0, 2], [1, 2], [0, 1]], [[1, 3], [2, 3], [1, 2]], [[2, 4], [3, 4], [2, 3]], [[3, 5], [4, 5], [3, 4]], [[4, 6], [5, 6], [4, 5]],
]

# Variant-independent constants (same for isea/ivea/rtea).
parallelepipedV = math.sqrt((5 - 2 * math.sqrt(5)) / 15)

# SDT area = 6 degrees = Pi/30 radians (one icosahedral sub-triangle's area).
SDT_AREA = 6 * DEG2RAD


# --- Per-variant spherical-triangle edges (icoVertexGreatCircle.ec:118-158) ---
# The ``radialVertex`` property setter assigns the three great-circle edge angles
# (and their cached trig) of the base sub-triangle. Only these differ between
# variants; the surrounding vector math is shared.
_VariantConsts = namedtuple("_VariantConsts", "AB AC BC cosAB cosAC sinAC cosBC")


# Per-variant vertex ordering (va, vb, vc) — indices into the SDT vertex array
# [mid(0), corner(1), center(2)]. The eC assigns a=vb, b=bIsA?va:vc, c=bIsA?vc:va.
# isea makes the CENTER the radial vertex A (vb=2); ivea makes the CORNER A (vb=1).
# (eC icoVertexGreatCircle.ec radialVertex property setter.)
_VERTEX_ORDER = {"isea": (0, 2, 1), "ivea": (0, 1, 2), "rtea": (1, 0, 2)}


def variant_consts(radial_vertex: str) -> "_VariantConsts":
    """Select the spherical-triangle edge angles for ``radial_vertex`` (the eC
    ``radialVertex`` cases). ``isea`` reproduces the historical module constants
    exactly (byte-identical); ``ivea`` swaps AB<->AC; ``rtea`` permutes all three
    (AB<-isea's AC, AC<-isea's BC, BC<-isea's AB) (``icoVertexGreatCircle.ec``)."""
    if radial_vertex == "isea":
        AB = math.acos(math.sqrt((PHI + 1) / 3))
        AC = math.atan(1 / PHI)
        BC = math.atan(2 / (PHI * PHI))
    elif radial_vertex == "ivea":
        AB = math.atan(1 / PHI)
        AC = math.acos(math.sqrt((PHI + 1) / 3))
        BC = math.atan(2 / (PHI * PHI))
    elif radial_vertex == "rtea":
        AB = math.atan(1 / PHI)
        AC = math.atan(2 / (PHI * PHI))
        BC = math.acos(math.sqrt((PHI + 1) / 3))
    else:
        raise ValueError(f"unknown radial_vertex {radial_vertex!r}")
    return _VariantConsts(
        AB, AC, BC,
        math.cos(AB), math.cos(AC), math.sin(AC), math.cos(BC),
    )


@dataclass(frozen=True)
class _Geometry:
    """Precomputed per-config geometry, built once by :meth:`_IcoVertexProjection.build_geometry`.

    Mirrors the oracle's ``_GridGeometry`` (rotated icosahedron vertices, face
    centroid tables, the authalic Clenshaw coefficient pair) and additionally
    carries ``azimuth_deg`` and ``authalic`` from the :class:`GridConfig`. The
    oracle reads those two off the ``grid`` at call time; here ``forward`` /
    ``inverse`` only receive the ``geom``, so the geometry must thread them. The
    ``radial_vertex`` + ``consts`` fields likewise thread the projection variant
    (eC ``radialVertex``) so the free-function kernel needs no module globals.
    """
    authalic_cp: tuple        # (geodetic->authalic coeffs, authalic->geodetic coeffs)
    ico_vertices: list        # 12 vertices [x, y, z] (DGGAL Y-up)
    face_centroids: list      # 20 face centroids [x, y, z]
    ico3rd_centroids: list
    ico6th_centroids: list
    ico3rd_mids: list
    ico56_center: list
    ico56_mids: list
    azimuth_deg: float        # vertex-2 azimuth offset (GridConfig)
    authalic: bool            # apply geodetic<->authalic conversion? (GridConfig)
    radial_vertex: str        # projection variant ("isea"|"ivea"|"rtea")
    consts: "_VariantConsts"  # per-variant spherical-triangle edges (variant_consts)


# --------------------------------------------------------------------------- #
# Authalic latitude (Karney 2022 Clenshaw series, authalic.ec)
# --------------------------------------------------------------------------- #
def precompute_coefficients(a, b, C):
    """Fold the 21 flat series coefficients into 6 Clenshaw coefficients for one
    ellipsoid (third flattening ``n``). Done once per geometry, per direction."""
    n = (a - b) / (a + b)
    d = n
    cp = [0.0] * 6
    cp[0] = (((((C[5]*n + C[4])*n + C[3])*n + C[2])*n + C[1])*n + C[0]) * d; d *= n
    cp[1] = ((((C[10]*n + C[9])*n + C[8])*n + C[7])*n + C[6]) * d; d *= n
    cp[2] = (((C[14]*n + C[13])*n + C[12])*n + C[11]) * d; d *= n
    cp[3] = ((C[17]*n + C[16])*n + C[15]) * d; d *= n
    cp[4] = (C[19]*n + C[18]) * d; d *= n
    cp[5] = C[20] * d
    return cp


def apply_coefficients(cp, phi):
    """Evaluate the Clenshaw series at latitude ``phi`` (the converted latitude)."""
    szeta = math.sin(phi)
    czeta = math.cos(phi)
    X = 2 * (czeta - szeta) * (czeta + szeta)  # 2*cos(2*phi) — keep factored form
    u0 = X * cp[5] + cp[4]
    u1 = X * u0 + cp[3]
    u0 = X * u1 - u0 + cp[2]
    u1 = X * u0 - u1 + cp[1]
    u0 = X * u1 - u0 + cp[0]
    return phi + 2 * szeta * czeta * u0


# --------------------------------------------------------------------------- #
# Vector math (3-component, DGGAL Y-up convention)
# --------------------------------------------------------------------------- #
def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def normalize(v):
    length = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if length < 1e-15:
        return [0.0, 0.0, 0.0]
    return [v[0] / length, v[1] / length, v[2] / length]


def vec_length(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


# --------------------------------------------------------------------------- #
# Quaternion math — used to rotate the base icosahedron into the configured
# orientation (yaw = orientation lon, pitch = orientation lat).
# --------------------------------------------------------------------------- #
def quaternion_from_yaw_pitch(yaw, pitch):
    cy = math.cos(yaw / 2)
    sy = math.sin(yaw / 2)
    cp = math.cos(pitch / 2)
    sp = math.sin(pitch / 2)
    return [cy * sp, sy * cp, -sy * sp, cy * cp]  # [x, y, z, w]


def quaternion_rotate(q, v):
    qx, qy, qz, qw = q
    vx, vy, vz = v
    tx = 2 * (qy * vz - qz * vy)
    ty = 2 * (qz * vx - qx * vz)
    tz = 2 * (qx * vy - qy * vx)
    return [
        vx + qw * tx + (qy * tz - qz * ty),
        vy + qw * ty + (qz * tx - qx * tz),
        vz + qw * tz + (qx * ty - qy * tx),
    ]


# --------------------------------------------------------------------------- #
# Coordinate conversions (DGGAL Y-up cartesian <-> geographic radians)
# --------------------------------------------------------------------------- #
def dggal_to_geo(x, y, z):
    p = math.sqrt(x * x + z * z)
    return [math.atan2(-y, p), math.atan2(x, -z)]


def geo_to_dggal_cart(lat, lon):
    clat = math.cos(lat)
    return [math.sin(lon) * clat, -math.sin(lat), -math.cos(lon) * clat]


def compute_vertices(orientation_lat_rad, orientation_lon_rad):
    """Build the 12 icosahedron vertices (two poles + two staggered pentagons),
    then rotate them into the configured orientation (ri5x6.ec)."""
    t = math.atan(0.5)
    sin_t = math.sin(t)
    cos_t = math.cos(t)
    step = TWO_PI / 5
    verts = [None] * 12

    verts[0] = [0.0, -1.0, 0.0]
    verts[11] = [0.0, 1.0, 0.0]

    for i in range(5):
        ta = (-270 + 72 * i) * DEG2RAD
        ba = ta + step / 2
        verts[1 + i] = [math.cos(ta) * cos_t, -sin_t, math.sin(ta) * cos_t]
        verts[6 + i] = [math.cos(ba) * cos_t, sin_t, math.sin(ba) * cos_t]

    q = quaternion_from_yaw_pitch(orientation_lon_rad, orientation_lat_rad)
    for i in range(12):
        verts[i] = quaternion_rotate(q, verts[i])

    return verts


def sqrt_one_minus_dot_over2(a, b):
    """sin(angle/2) between unit vectors ``a`` and ``b``, computed stably via the
    midpoint-cross trick (icoVertexGreatCircle.ec ``sqrtOneMinusDotOver2``)."""
    mid = normalize([(a[0]+b[0])/2, (a[1]+b[1])/2, (a[2]+b[2])/2])
    c = cross(a, mid)
    D = vec_length(c)
    if D < 1e-8:
        D = vec_length([a[0]-b[0], a[1]-b[1], a[2]-b[2]]) / 2
    return D


def spherical_tri_area(A, B, C):
    """Signed spherical excess of triangle ABC via the midpoint-cross form
    (``sphericalTriArea``, ``ri5x6.ec:224``) — numerically stable for the tiny
    SDT triangles."""
    mid_ab = normalize([(A[0]+B[0])/2, (A[1]+B[1])/2, (A[2]+B[2])/2])
    mid_bc = normalize([(B[0]+C[0])/2, (B[1]+C[1])/2, (B[2]+C[2])/2])
    mid_ca = normalize([(C[0]+A[0])/2, (C[1]+A[1])/2, (C[2]+A[2])/2])
    cr = cross(mid_bc, mid_ca)
    return math.asin(max(-1, min(1, dot(mid_ab, cr)))) * 2


# --------------------------------------------------------------------------- #
# Planar barycentric helpers (``cartesianToBary``/``baryToCartesian``,
# ``ri5x6.ec:211``/``:237`` — NOT ``barycentric5x6.ec``, a separate Goldberg
# projection class that only calls them, with a different ``knownOneOverDet``)
# --------------------------------------------------------------------------- #
def cartesian_to_bary(p, p1, p2, p3, known_one_over_det):
    """Barycentric coords of planar ``p`` in triangle (p1,p2,p3). ``known_one_over_det``
    short-circuits the determinant when the caller knows it (e.g. -6 or -1)."""
    d31x = p1[0] - p3[0]; d31y = p1[1] - p3[1]
    d23x = p3[0] - p2[0]; d23y = p3[1] - p2[1]
    d3px = p[0] - p3[0]; d3py = p[1] - p3[1]
    o_det = known_one_over_det if known_one_over_det else (1 / (d23x * d31y - d23y * d31x))
    b0 = (d23x * d3py - d23y * d3px) * o_det
    b1 = (d31x * d3py - d31y * d3px) * o_det
    return [b0, b1, 1 - b0 - b1]


def bary_to_cartesian(b, p1, p2, p3):
    return [
        b[0] * p1[0] + b[1] * p2[0] + b[2] * p3[0],
        b[0] * p1[1] + b[1] * p2[1] + b[2] * p3[1],
    ]


def slerp_angle(p0, p1, distance, movement):
    """Spherical interpolation by ``movement`` radians along the ``distance``-long
    great-circle arc from ``p0`` toward ``p1`` (``slerpAngle``, ``ri5x6.ec:176``).
    The ``abs(s_distance) < 1e-15`` short-circuit is a safety net this port adds;
    the eC divides unconditionally (its distances are fixed SDT edges, never ~0)."""
    s_distance = math.sin(distance)
    if abs(s_distance) < 1e-15:
        return list(p0)
    o_o_sin = 1 / s_distance
    l0 = math.sin(distance - movement)
    l1 = math.sin(movement)
    return [
        (l0 * p0[0] + l1 * p1[0]) * o_o_sin,
        (l0 * p0[1] + l1 * p1[1]) * o_o_sin,
        (l0 * p0[2] + l1 * p1[2]) * o_o_sin,
    ]


# --------------------------------------------------------------------------- #
# Face & sub-triangle finding.
# NOTE: this nearest-centroid + early-accept-threshold method is NOT how vendored
# DGGAL v0.06 finds the face/sub-triangle -- that iterates faces calling the
# plane-sign tests ``vertexWithinSphericalTriPlanes``/``vertexWithinSphericalTri``
# (``ri5x6.ec:362``/``:333``). This centroid-nearest form (and the three
# ``early_accept`` constants below, which are absent from the v0.06 tree) is ported
# from the newer DGGAL / ``igeo7`` lineage this library descends from. Behaviourally
# equivalent (both identify the containing region; verified by the fuzz), but a
# reader checking against v0.06 will find a different algorithm -- so cite the
# lineage, not v0.06, and treat the constants as empirically-tuned thresholds.
# --------------------------------------------------------------------------- #
# Hoisted to module scope 2026-08-21: these three are constants, and find_face /
# get_face are the hottest functions in the library -- forward() calls find_face
# once per geographic point and Grid.neighbors' k-ring calls forward() once per
# cell edge, so rebuilding the literals per call allocated thousands of throwaway
# lists/dicts per sweep. Read-only, so hoisting is behaviour-identical.
_N_MAP = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
_S_MAP = (19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5)
_RHOMBUS_FACES = {
    0: (0, 5), 2: (1, 6), 4: (2, 7), 6: (3, 8), 8: (4, 9),
    1: (10, 15), 3: (11, 16), 5: (12, 17), 7: (13, 18), 9: (14, 19),
}


def find_face(geom, P):
    """Pick the icosahedron face whose centroid is closest (max dot) to point ``P``.
    The pole-sign split chooses the 15-face hemisphere to scan; the early-accept
    threshold short-circuits once a face is unambiguously the nearest. See the
    section note above on provenance (newer lineage, not v0.06's plane-test)."""
    best = 0
    best_dot = -math.inf

    pole_dot = geom.ico_vertices[0][0] * P[0] + geom.ico_vertices[0][1] * P[1] + geom.ico_vertices[0][2] * P[2]
    m = _N_MAP if pole_dot > 0 else _S_MAP

    early_accept = 0.934172358962715696451118623548

    for i in range(15):
        face = m[i]
        c = geom.face_centroids[face]
        d = c[0] * P[0] + c[1] * P[1] + c[2] * P[2]
        if d > best_dot:
            best_dot = d
            best = face
            if d > early_accept:
                break
    return best


def find_sub_tri(geom, face, v):
    """Locate the SDT sub-triangle (0..5) within ``face`` containing direction ``v``:
    first the nearest third, then the nearest half within it. Same provenance
    caveat as :func:`find_face` (newer lineage; the two ``early_accept`` constants
    are empirically-tuned and absent from vendored v0.06)."""
    best3rd = 0
    best6th = 0
    best_dot = -math.inf

    early_accept_3rd = 0.9802668134226932631948092150332116
    early_accept_6th = 0.9829160426524585629980328985973081873244

    for i in range(3):
        c = geom.ico3rd_centroids[face][i]
        d = c[0] * v[0] + c[1] * v[1] + c[2] * v[2]
        if d > best_dot:
            best_dot = d
            best3rd = i
            if d > early_accept_3rd:
                break

    best_dot = -math.inf
    for i in range(2):
        c = geom.ico6th_centroids[face][best3rd][i]
        d = c[0] * v[0] + c[1] * v[1] + c[2] * v[2]
        if d > best_dot:
            best_dot = d
            best6th = i
            if d > early_accept_6th:
                break
    return 2 * best3rd + best6th


# --------------------------------------------------------------------------- #
# The equal-area vector projection itself (icoVertexGreatCircle.ec)
# --------------------------------------------------------------------------- #
def forward_vector(v, A, B, C, pA, pB, pC):
    """Forward Snyder map of unit vector ``v`` in spherical triangle (A,B,C) to
    planar barycentrics in triangle (pA,pB,pC). Equal-area: the planar area
    fraction equals the spherical area fraction (``area_abp / SDT_AREA``)."""
    c1 = cross(A, v)
    c2 = cross(B, C)
    p = normalize(cross(c1, c2))

    if dot(p, B) < 0:
        p = [-p[0], -p[1], -p[2]]

    area_abp = max(0.0, spherical_tri_area(A, B, p))
    h = sqrt_one_minus_dot_over2(A, v) / sqrt_one_minus_dot_over2(A, p)

    b0 = 1 - h
    b2 = min(h, h * area_abp / SDT_AREA)
    b1 = h - b2

    return bary_to_cartesian([b0, b1, b2], pA, pB, pC)


def inverse_vector(pi, pA, pB, pC, A, B, C, b_is_a, consts):
    """Inverse Snyder map: planar point ``pi`` (in triangle pA,pB,pC) back to a unit
    vector in spherical triangle (A,B,C). The closed-form solves for the foot ``p``
    on edge BC then slerps from A; a degenerate fallback handles the near-edge /
    near-singular cases (icoVertexGreatCircle.ec). ``consts`` supplies the per-variant
    spherical-triangle edges (see :func:`variant_consts`)."""
    b = cartesian_to_bary(pi, pA, pB, pC, -6)

    if b[0] > 1 - 1e-15:
        return list(A)
    if b[1] > 1 - 1e-15:
        return list(B)
    if b[2] > 1 - 1e-15:
        return list(C)

    h = 1 - b[0]
    b2oh = b[2] / h
    b2oh_abc = b2oh * SDT_AREA
    half_c = math.sin(b2oh_abc / 2)
    half_c2 = half_c * half_c
    CC = 2 * half_c2
    S = 2 * half_c * math.sqrt(1 - half_c2)
    c01 = consts.cosAB if b_is_a else consts.cosBC
    c12 = consts.cosAC
    c20 = consts.cosBC if b_is_a else consts.cosAB
    s12 = consts.sinAC
    f = S * parallelepipedV + CC * (c01 * c12 - c20)
    g = CC * s12 * (1 + c01)
    f2 = f * f; g2 = g * g; gf = g * f
    numerator = s12 * (f2 - g2) - 2 * gf * c12
    divisor = s12 * (f2 + g2)

    if abs(numerator) > 1e-9 and abs(divisor) > 1e-9:
        o_o_divisor = 1.0 / divisor
        ap = max(0.0, numerator * o_o_divisor)
        bp = min(1.0, 2 * gf * o_o_divisor)
        px = ap * B[0] + bp * C[0]
        py = ap * B[1] + bp * C[1]
        pz = ap * B[2] + bp * C[2]

        av = A[0] * px + A[1] * py + A[2] * pz
        bv = 1 + h * h * (av - 1)
        bvp = h * math.sqrt((1 + bv) / (1 + av))
        avp = bv - av * bvp

        return [avp * A[0] + bvp * px, avp * A[1] + bvp * py, avp * A[2] + bvp * pz]
    else:
        # Degenerate fallback: the eC recomputes the AB/AC edge angles here. These
        # MUST come from the variant ``consts`` (isea and ivea swap AB<->AC) — for
        # isea consts.AB/AC equal these formulas exactly (byte-identical), for ivea
        # they are swapped. ``alpha`` (vertex angle at A) is pi/2 for isea and ivea.
        beta = consts.AB
        gamma = consts.AC
        alpha = math.pi / 2
        b1pb2 = b[1] + b[2]
        up_over_up_pvp = 0 if b1pb2 < 1e-11 else (b[1] if b_is_a else b[2]) / b1pb2
        rho_plus_delta = beta + gamma - up_over_up_pvp * SDT_AREA
        area_abd = rho_plus_delta + alpha - math.pi

        if abs(area_abd) < 1e-11:
            D = list(B) if b_is_a else list(C)
            BD = consts.AB
        elif abs(area_abd - SDT_AREA) < 1e-13:
            D = list(C) if b_is_a else list(B)
            BD = consts.BC
        else:
            AD = 2 * math.atan2(g, f)
            D = slerp_angle(B, C, consts.AC, AD)
            BD = math.acos(max(-1, min(1, dot(A, D))))

        x = 2 * math.asin(min(1, (1 - b[0]) * math.sin(BD / 2)))
        return slerp_angle(A, D, BD, x)


def _sub_tri_frame(geom, face, sub_tri):
    """Resolve SDT sub-triangle ``sub_tri`` on ``face`` to its projection frame:
    the spherical-vertex triple and the 5x6-vertex triple, both reordered into the
    eC ``(a, b, c)`` roles for this variant, plus the ``b_is_a`` flag.

    forward_ico_face and inverse_ico_face differ ONLY in how they find ``sub_tri``
    (spherical nearest-centroid vs planar barycentrics); everything downstream is
    this shared frame selection. Pure list construction + integer index selection
    (no float arithmetic), so ISEA output stays bit-identical. eC: ``a = vb;
    b = bIsA?va:vc; c = bIsA?vc:va`` (va/vb/vc per variant, isea 0,2,1 / ivea 0,1,2)
    with ``bIsA = (radial_vertex=='ivea') ^ (sub_tri in {0,3,4})``."""
    ii0, ii1, ii2 = FACE_VERTICES[face]
    v1 = geom.ico_vertices[ii0]; v2 = geom.ico_vertices[ii1]; v3 = geom.ico_vertices[ii2]
    p1 = FACE_5X6[face][0]; p2 = FACE_5X6[face][1]; p3 = FACE_5X6[face][2]
    tri3rd = sub_tri >> 1

    p5x6 = [
        geom.ico56_mids[face][tri3rd],
        p3 if sub_tri in (0, 2) else (p2 if sub_tri in (1, 4) else p1),
        geom.ico56_center[face],
    ]
    v3d = [
        geom.ico3rd_mids[face][tri3rd],
        v3 if sub_tri in (0, 2) else (v2 if sub_tri in (1, 4) else v1),
        geom.face_centroids[face],
    ]

    b_is_a = (geom.radial_vertex == "ivea") ^ (sub_tri in (0, 3, 4))
    va, vb, vc = _VERTEX_ORDER[geom.radial_vertex]
    a = vb
    b_idx = va if b_is_a else vc
    c_idx = vc if b_is_a else va
    return [v3d[a], v3d[b_idx], v3d[c_idx]], [p5x6[a], p5x6[b_idx], p5x6[c_idx]], b_is_a


def forward_ico_face(geom, face, v):
    """Project direction ``v`` onto ``face``: find its SDT, then map within the
    centre-anchored sub-triangle to 5x6 planar coords (icoVertexGreatCircle.ec)."""
    sub_tri = find_sub_tri(geom, face, v)
    v3d, p5x6, _ = _sub_tri_frame(geom, face, sub_tri)
    return forward_vector(v, v3d[0], v3d[1], v3d[2], p5x6[0], p5x6[1], p5x6[2])


def inverse_ico_face(geom, face, v56):
    """Inverse of :func:`forward_ico_face`: from a 5x6 planar point on ``face``,
    pick the sub-triangle by planar barycentrics, then run the inverse vector map."""
    p1 = FACE_5X6[face][0]; p2 = FACE_5X6[face][1]; p3 = FACE_5X6[face][2]

    b = cartesian_to_bary(v56, p1, p2, p3, -1)
    if b[0] <= b[1] and b[0] <= b[2]:
        sub_tri = 0 if b[1] < b[2] else 1
    elif b[1] <= b[0] and b[1] <= b[2]:
        sub_tri = 2 if b[0] < b[2] else 3
    else:
        sub_tri = 4 if b[0] < b[1] else 5

    v3d, p5x6, b_is_a = _sub_tri_frame(geom, face, sub_tri)
    return inverse_vector(v56, p5x6[0], p5x6[1], p5x6[2], v3d[0], v3d[1], v3d[2], b_is_a, geom.consts)


class _IcoVertexProjection:
    """The shared icosahedral great-circle equal-area projection — a ``Projection`` impl.

    One kernel, parameterized by :attr:`radial_vertex` (eC ``radialVertex``): the thin
    variant subclasses (:class:`~py4dggs.projections.isea.ISEAProjection`,
    :class:`~py4dggs.projections.ivea.IVEAProjection`) set that class attribute and inherit
    everything here. ``build_geometry`` precomputes the rotated icosahedron + authalic
    coefficients + the variant edge constants once per :class:`GridConfig`; ``forward`` /
    ``inverse`` are pure given that geometry. Output is bit-identical to the ``igeo7``
    oracle by construction.
    """

    radial_vertex: str = "isea"  # variant selector; overridden by subclasses

    def build_geometry(self, config: GridConfig) -> Any:
        """Compute the icosahedron + authalic data for ``config`` (= oracle's
        ``build_geometry(grid)``). The authalic coefficients are WGS84-fixed; the
        ``authalic`` flag only toggles whether the conversion is applied, in
        forward/inverse. ``azimuth_deg`` and ``authalic`` are carried on the geom
        so the pure forward/inverse can read them without the config."""
        authalic_cp = (
            precompute_coefficients(WGS84_A, WGS84_B, Cxiphi),
            precompute_coefficients(WGS84_A, WGS84_B, Cphixi),
        )

        ico_vertices = compute_vertices(
            config.orientation_lat_deg * DEG2RAD, config.orientation_lon_deg * DEG2RAD
        )

        face_centroids = [None] * 20
        for f in range(20):
            i0, i1, i2 = FACE_VERTICES[f]
            v0 = ico_vertices[i0]; v1 = ico_vertices[i1]; v2 = ico_vertices[i2]
            face_centroids[f] = normalize([
                v0[0] + v1[0] + v2[0],
                v0[1] + v1[1] + v2[1],
                v0[2] + v1[2] + v2[2],
            ])

        ico3rd_centroids = [None] * 20
        ico6th_centroids = [None] * 20
        ico3rd_mids = [None] * 20
        ico56_center = [None] * 20
        ico56_mids = [None] * 20

        for f in range(20):
            ii0, ii1, ii2 = FACE_VERTICES[f]
            v1 = ico_vertices[ii0]; v2 = ico_vertices[ii1]; v3 = ico_vertices[ii2]
            c = face_centroids[f]
            p1 = FACE_5X6[f][0]; p2 = FACE_5X6[f][1]; p3 = FACE_5X6[f][2]

            ico56_center[f] = [(p1[0]+p2[0]+p3[0])/3, (p1[1]+p2[1]+p3[1])/3]

            ico3rd_centroids[f] = [
                normalize([c[0]+v2[0]+v3[0], c[1]+v2[1]+v3[1], c[2]+v2[2]+v3[2]]),
                normalize([c[0]+v3[0]+v1[0], c[1]+v3[1]+v1[1], c[2]+v3[2]+v1[2]]),
                normalize([c[0]+v1[0]+v2[0], c[1]+v1[1]+v2[1], c[2]+v1[2]+v2[2]]),
            ]
            ico3rd_mids[f] = [
                normalize([(v2[0]+v3[0])/2, (v2[1]+v3[1])/2, (v2[2]+v3[2])/2]),
                normalize([(v3[0]+v1[0])/2, (v3[1]+v1[1])/2, (v3[2]+v1[2])/2]),
                normalize([(v1[0]+v2[0])/2, (v1[1]+v2[1])/2, (v1[2]+v2[2])/2]),
            ]
            ico56_mids[f] = [
                [(p2[0]+p3[0])/2, (p2[1]+p3[1])/2],
                [(p3[0]+p1[0])/2, (p3[1]+p1[1])/2],
                [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2],
            ]

            ico6th_centroids[f] = [None] * 3
            for tt in range(3):
                m = ico3rd_mids[f][tt]
                if tt == 0:
                    va, vb = v3, v2
                elif tt == 1:
                    va, vb = v3, v1
                else:
                    va, vb = v2, v1
                ico6th_centroids[f][tt] = [
                    normalize([c[0]+m[0]+va[0], c[1]+m[1]+va[1], c[2]+m[2]+va[2]]),
                    normalize([c[0]+m[0]+vb[0], c[1]+m[1]+vb[1], c[2]+m[2]+vb[2]]),
                ]

        return _Geometry(
            authalic_cp=authalic_cp,
            ico_vertices=ico_vertices,
            face_centroids=face_centroids,
            ico3rd_centroids=ico3rd_centroids,
            ico6th_centroids=ico6th_centroids,
            ico3rd_mids=ico3rd_mids,
            ico56_center=ico56_center,
            ico56_mids=ico56_mids,
            azimuth_deg=config.azimuth_deg,
            authalic=config.authalic,
            radial_vertex=self.radial_vertex,
            consts=variant_consts(self.radial_vertex),
        )

    def geodetic_to_authalic(self, geom: Any, lat_rad: float) -> float:
        """Geodetic latitude -> authalic latitude (authalic.ec; oracle:
        ``lat_geodetic_to_authalic(geom, phi)``). Karney coeffs live on ``geom``."""
        return apply_coefficients(geom.authalic_cp[0], lat_rad)

    def authalic_to_geodetic(self, geom: Any, lat_rad: float) -> float:
        """Authalic latitude -> geodetic latitude (authalic.ec; oracle:
        ``lat_authalic_to_geodetic(geom, phi)``). Karney coeffs live on ``geom``."""
        return apply_coefficients(geom.authalic_cp[1], lat_rad)

    def get_face(self, geom: Any, x: float, y: float) -> int:
        """Derive the face id from a 5x6 planar point alone (oracle: ``get_face``).
        Nudges edge/corner points inward by ``epsilon`` so the rhombus lookup is
        unambiguous; returns -1 if the point lies outside the 5x6 layout."""
        epsilon = 1e-11
        if x < 0 or (y > x and x < 5 - epsilon):
            x += epsilon
        elif x > 5 or (y < x and x > 0 + epsilon):
            x -= epsilon
        if y < 0 or (x > y and y < 6 - epsilon):
            y += epsilon
        elif y > 6 or (x < y and y > 0 + epsilon):
            y -= epsilon

        if 0 <= x <= 5 and 0 <= y <= 6:
            ix = max(0, min(4, math.floor(x)))
            iy = max(0, min(5, math.floor(y)))
            if iy == ix or iy == ix + 1:
                rhombus = ix + iy
                top = x - ix > y - iy
                if rhombus in _RHOMBUS_FACES:
                    t, bt = _RHOMBUS_FACES[rhombus]
                    return t if top else bt
        return -1

    def forward(self, geom: Any, lat: float, lon: float) -> PlanarPoint:
        """Geographic (degrees) -> planar ``PlanarPoint(face, x, y)`` in 5x6 space
        (oracle: ``project_point``). Applies the vertex-2 azimuth offset, the
        geodetic->authalic conversion (if enabled), finds the face, then the
        equal-area face projection."""
        lon_rad = (lon - geom.azimuth_deg) * DEG2RAD          # vertex2Azimuth (ri5x6.ec:397)
        lat_rad = lat * DEG2RAD
        lat_conv = self.geodetic_to_authalic(geom, lat_rad) if geom.authalic else lat_rad
        v3d = geo_to_dggal_cart(lat_conv, lon_rad)
        face = find_face(geom, v3d)
        x, y = forward_ico_face(geom, face, v3d)
        return PlanarPoint(face, x, y)

    def inverse(self, geom: Any, p: PlanarPoint) -> GeoPoint:
        """Planar ``PlanarPoint(face, x, y)`` -> geographic ``GeoPoint`` (degrees)
        (oracle: ``unproject_point``). A negative ``face`` means "derive from x,y"
        (matches the oracle's ``face=None`` path). Edge-wrap fixups, the inverse
        face projection, then the authalic->geodetic conversion + azimuth offset.

        Deliberate omission: the eC ``inverse`` calls ``fixPoles`` (``ri5x6.ec:463-509``)
        to snap the four near-pole 5x6 points to an exact +-90 deg with a chosen
        longitude; this port does not, because the topology never feeds an exact
        pole point to the projection. This is the mechanism behind the documented
        IVEA/RTEA pole / vertex-0-meridian conformance xfails (a point landing
        EXACTLY on that singularity would get the raw ``atan2`` longitude here
        rather than the snapped one). Likewise the authalic gate and the post-wrap
        ``azimuth_deg`` add assume the six registered grids' config (authalic on,
        azimuth 0); both are latent if the kernel is reused with a different config."""
        vx = p.x
        vy = p.y
        if vx < 0 and vy < 1:
            vx += 5; vy += 5
        elif vx > 5 and vy > 5:
            vx -= 5; vy -= 5
        elif vx > 5 and vy > 5 - 1e-10:
            vx -= 5; vy = 0
        elif vy > 6 and vx > 5 - 1e-10:
            vy -= 5; vx = 0
        elif vy < 0 and vx < 1e-10:
            vx = 5; vy += 5

        face = p.face
        if face < 0:
            face = self.get_face(geom, vx, vy)
        if face < 0:
            return GeoPoint(0.0, 0.0)

        p3d = inverse_ico_face(geom, face, [vx, vy])
        lat, lon = dggal_to_geo(p3d[0], p3d[1], p3d[2])

        lat_geo = self.authalic_to_geodetic(geom, lat) if geom.authalic else lat

        lon_n = lon
        if lon_n > math.pi:
            lon_n -= TWO_PI
        if lon_n < -math.pi:
            lon_n += TWO_PI

        return GeoPoint(lat_geo * RAD2DEG, lon_n * RAD2DEG + geom.azimuth_deg)  # + vertex2Azimuth (ri5x6.ec:547)
