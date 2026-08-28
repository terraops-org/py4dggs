"""I3HIndexing — the minimal indexing for the ISEA3H geometry slice.

The I3H cell (levelI9R, rootRhombus, rhombusIX, subHex) is bit-packed into the
Zone value using DGGAL's exact I3HZone uint64 layout (RI3H.ec:855-861), so the
value EQUALS pydggal's DGGRSZone int. Geometry only: hierarchy (parent/children)
and the canonical rhombic text-id ("C2-23-C") are deferred to sub-project A. The
digit-path protocol is honoured by carrying the whole cell as ``base`` with empty
``digits``.
"""
from __future__ import annotations

import re

from py4dggs.types import NULL_TEXT, NULL_ZONE, InvalidZoneError

_RIX_BITS = 51
_RIX_MASK = (1 << _RIX_BITS) - 1
# RI3H.ec:855-861 declares I3HZone as C bitfields -- `uint levelI9R:5:57`,
# `uint rootRhombus:4:53`, `uint64 rhombusIX:51:2`, `uint subHex:2:0` -- and a C
# bitfield TRUNCATES on assignment: writing 16 to a 4-bit field stores 0. It does
# NOT carry into the neighbouring field. Masking only `rix` (as this module did
# until 2026-08-21) let an over-range level/root/sub_hex corrupt an ADJACENT
# field instead, fabricating a plausible-looking zone -- e.g. pack_i3h(0, 16, 0, 0)
# unpacked as levelI9R=1. Mask every field so the packing matches the eC.
_LEVEL_MASK = (1 << 5) - 1
_ROOT_MASK = (1 << 4) - 1
_SUB_HEX_MASK = (1 << 2) - 1

_TEXTID_RE = re.compile(r"([A-Za-z])([0-9A-Fa-f]+)-([0-9A-Fa-f]+)-([A-Za-z])")


def pack_i3h(level_i9r: int, root: int, rix: int, sub_hex: int) -> int:
    """(levelI9R, rootRhombus, rhombusIX, subHex) -> DGGAL-layout uint64.

    Every field is truncated to its declared bit width, exactly as the eC
    bitfield does (see the mask constants above) -- an out-of-range field
    wraps within itself and never disturbs its neighbour.
    """
    return (
        ((level_i9r & _LEVEL_MASK) << 57)
        | ((root & _ROOT_MASK) << 53)
        | ((rix & _RIX_MASK) << 2)
        | (sub_hex & _SUB_HEX_MASK)
    )


def unpack_i3h(value: int) -> tuple[int, int, int, int]:
    """Inverse of :func:`pack_i3h`."""
    sub_hex = value & 0x3
    rix = (value >> 2) & _RIX_MASK
    root = (value >> 53) & 0xF
    level_i9r = (value >> 57) & 0x1F
    return level_i9r, root, rix, sub_hex


def _ilrc_from_lrti(level, root, ix):
    """iLRCFromLRtI (RI9R.ec:710-729): (levelI9R, root, ix) -> global (row, col).
    Returns (level, row, col) or None for the eC's -1 sentinel."""
    if 0 <= level <= 16 and 0 <= root <= 9:
        p = 3 ** level
        if 0 <= ix < p * p:
            row_op = (root + 1) >> 1
            col_op = root >> 1
            ix_op = ix // p
            row = row_op * p + ix_op
            col = (col_op - ix_op) * p + ix
            return level, row, col
    return None


def _from_i9r(level, row, col, sub_hex, pole):
    """fromI9R (RI3H.ec:951-963): global (row, col) -> canonical (levelI9R, root,
    rhombusIX, subHex). ``pole`` is the pole root (10/11) or 0. None on nullZone."""
    p = 3 ** level
    row_op = row // p
    col_op = col // p
    root = pole if pole else row_op + col_op
    y = row - row_op * p
    x = col - col_op * p
    ix = 0 if pole else y * p + x
    if not (0 <= sub_hex <= 3) or root > 11 or (
        root < 10 and (row_op < col_op or row_op - col_op > 1 or y >= p or x >= p)
    ):
        return None
    return (level, root, ix, sub_hex)


def _validate_lrc(level_i9r, root_rhombus, row, col, sub_hex):
    """validate (RI3H.ec:965-975): the guard fromZoneID applies before accepting."""
    p = 3 ** level_i9r
    row_op = row // p
    col_op = col // p
    y = row - row_op * p
    x = col - col_op * p
    root = row_op + col_op
    if not (0 <= sub_hex <= 3) or root_rhombus > 11 or (
        root_rhombus <= 9 and (root != root_rhombus or row_op < col_op or row_op - col_op > 1 or y >= p or x >= p)
    ):
        return False
    return True


def is_pentagon_cell(rix: int, sub_hex: int) -> bool:
    """Shared pentagon predicate: rix==0 (the root-rhombus's polar corner
    cell) with sub_hex in {0,1} (the A/B centroid children of that corner)
    is a pentagon. Used by I3HIndexing.is_pentagon and by hex_a3.py's
    sub-zone dispatch/count, which both need the same rule."""
    return rix == 0 and sub_hex <= 1


def absolute_level(level_i9r: int, sub_hex: int) -> int:
    """The absolute I3H level (resolution) a cell sits at: sub_hex==0 (A)
    is on the even level `level_i9r*2`; sub_hex in {1,2,3} (B/C/D) are one
    level deeper. Used by I3HIndexing.resolution and by hex_a3.py's
    sub-zone code, which both need the same "parent's own absolute level"
    quantity."""
    return level_i9r * 2 + (1 if sub_hex > 0 else 0)


class I3HIndexing:
    """Minimal ``Indexing`` for ISEA3H geometry (see module docstring)."""

    max_resolution = 33  # DGGAL max: levelI9R<=16 -> level<=33 (rix p*p<2^51); res 34 -> nullZone

    def encode(self, base: int, digits: list[int]) -> int:
        if digits:
            raise InvalidZoneError("I3H carries the whole cell in base; digits must be empty")
        return base

    def decode(self, value: int) -> tuple[int, list[int]]:
        return value, []

    def resolution(self, value: int) -> int:
        level_i9r, _, _, sub_hex = unpack_i3h(value)
        return absolute_level(level_i9r, sub_hex)

    def base_cell(self, value: int) -> int:
        return unpack_i3h(value)[1]

    def is_pentagon(self, value: int) -> bool:
        _, _, rix, sub_hex = unpack_i3h(value)
        return is_pentagon_cell(rix, sub_hex)

    def to_text(self, value: int) -> str:
        """Canonical rhombic id (``getZoneID``, RI3H.ec:935-949):
        ``{levelChar}{root:X}-{rhombusIX:X}-{subHexChar}`` — levelChar = 'A'+levelI9R
        (the I9R level, not the DGGS resolution), root/ix uppercase hex (poles
        print root 10/11 as 'A'/'B'), subHexChar = 'A'+subHex.

        The ``nullZone`` sentinel prints as ``"(null)"``, exactly as the eC does
        (``getZoneTextID`` special-cases it before decoding the bit fields).
        Without this branch the sentinel's all-ones fields stringify into a
        plausible-looking id (```F-7FFFFFFFFFFFF-D``), which silently hides the
        one signal DGGAL gives that a zone does not exist — and this sentinel is
        genuinely reachable from ``zone_from_geo`` very near the poles at the
        coarsest resolutions (see ``hex_a3.py``'s note at ``NULL_ZONE``)."""
        if value == NULL_ZONE:
            return NULL_TEXT
        l9r, root, rix, sub_hex = unpack_i3h(value)
        return "%c%X-%X-%c" % (ord("A") + l9r, root, rix, ord("A") + sub_hex)

    # NOTE: the I3H hierarchy is NON-congruent and geometric, so it is NOT a
    # digit-path on Indexing — it lives on the topology (`hex_a3.parents`/
    # `children`/`centroid_parent`/`is_centroid_child`) and is reached via
    # Grid-level dispatch (A2). Hence no parent/num_children/child_digits here.

    def from_text(self, text: str) -> int:
        """Parse a canonical rhombic id (``fromZoneID``, RI3H.ec:908-931): parse the
        4 tokens, resolve (row,col) via iLRCFromLRtI (non-polar) or the level-letter
        (polar), validate, re-derive canonical (root,ix) via fromI9R, then require the
        result to re-serialize to the *exact* input (round-trip canonicalization, the
        eC strcmp). Raises InvalidZoneError on malformed / out-of-range / non-canonical.

        ``"(null)"`` round-trips back to the ``nullZone`` sentinel, so
        ``from_text(to_text(v)) == v`` holds for every IN-RANGE value (levelI9R
        <= 16, i.e. resolution <= max_resolution).

        It does NOT hold past max_resolution. ``children()`` can emit such zones
        because DGGAL itself does -- pydggal's ``getZoneChildren`` on a res-33
        ISEA3H zone returns seven res-34 zones, matching ours exactly (verified
        2026-08-21) -- and their text ids carry an out-of-range level letter that
        pydggal's own ``getZoneFromTextID`` also rejects. Mirroring the engine is
        the policy; the round-trip claim is scoped rather than the behaviour changed.
        Note the eC is laxer here: DGGAL's ``getZoneFromTextID`` returns nullZone
        for *any* unparseable string (``"null"``, ``""``, garbage). This port
        deliberately keeps raising ``InvalidZoneError`` for those instead — the
        same established deviation as everywhere else in this class, where DGGAL
        signals failure with a sentinel return and we signal it with an
        exception. Only the exact canonical spelling maps to the sentinel."""
        if text == NULL_TEXT:
            return NULL_ZONE
        m = _TEXTID_RE.fullmatch(text)
        if not m:
            raise InvalidZoneError(f"not a valid I3H text-id: {text!r}")
        level_char = m.group(1).upper()
        root = int(m.group(2), 16)
        ix = int(m.group(3), 16)
        sub_hex = ord(m.group(4).upper()) - ord("A")
        if root < 10:
            r = _ilrc_from_lrti(ord(level_char) - ord("A"), root, ix)
            if r is None:
                raise InvalidZoneError(f"invalid I3H text-id: {text!r}")
            l9r, row, col = r
        elif root <= 12 and ix == 0 and "A" <= level_char <= "Q":
            l9r, row, col = ord(level_char) - ord("A"), 0, 0
        else:
            raise InvalidZoneError(f"invalid I3H text-id: {text!r}")
        if not _validate_lrc(l9r, root, row, col, sub_hex):
            raise InvalidZoneError(f"invalid I3H text-id: {text!r}")
        res = _from_i9r(l9r, row, col, sub_hex, root if root > 9 else 0)
        if res is None:
            raise InvalidZoneError(f"invalid I3H text-id: {text!r}")
        value = pack_i3h(*res)
        if self.to_text(value) != text:  # round-trip canonicalization (eC strcmp)
            raise InvalidZoneError(f"non-canonical I3H text-id: {text!r}")
        return value

