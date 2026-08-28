# src/py4dggs/indexings/z7.py
"""Z7 indexing: a 64-bit cell id = 4-bit base cell (0-11) + up to twenty 3-bit
direction digits (0-6; 7 = unused). Resolution is implicit (first 7 digit).

Bit layout (verified against repos/dggal-v0.06/src/dggrs/RI7H_Z7.ec):
  - bits 63-60 : 4 bits   base cell (rootPentagon field, :4:60 in eC struct)
  - bits 59-57 : 3 bits   direction digit at resolution 1
  - bits 56-54 : 3 bits   direction digit at resolution 2
  - ...
  - bits 2-0   : 3 bits   direction digit at resolution 20

Each direction digit is 0-6. Digits beyond the cell's resolution are filled
with 7 (sentinel/unused). Resolution is found by scanning for the first 7.

Re-expressed from:
  - igeo7/_bits.py         — encode/decode/resolution/hierarchy
  - igeo7/zone.py          — text id and pentagon canonical check
"""
from __future__ import annotations

from py4dggs.types import InvalidZoneError

# Default id: base cell 0, all 20 direction digits = 7.
# 4 bits of 0 followed by 60 bits of 1 = 0x0FFFFFFFFFFFFFFF.
_DEFAULT = 0x0FFFFFFFFFFFFFFF
_MAX_RES = 20


def _bit_offset(res: int) -> int:
    """Bit offset (in the full 64-bit integer) of direction digit at *res*.

    The eC iterates ``shift = 19 * 3`` down to ``0`` inside the 60-bit
    ``ancestry`` field; since that field sits at bit 0 of the 64-bit int,
    the same formula applies directly to the integer.
    """
    return (20 - res) * 3


def _pentagon_omit(base: int) -> int:
    """The direction digit that is **invalid** as a child of a pentagon rooted
    at *base*.

    Pentagons live only on the all-zero centre path; the hemisphere of the
    root cell fixes the one direction digit that was removed when the grid was
    constructed: 2 for north-hemisphere bases (0-5), 5 for south (6-11).
    Verified against the eC's ``fromTextID`` rejection check:
      ``c == (south ? '5' : '2')``  where  ``south = root >= 6``.
    """
    return 2 if base <= 5 else 5


class Z7Indexing:
    """Pure-Python implementation of the Indexing protocol for the Z7 scheme.

    This is a re-expression of the bit-packing in DGGAL's ``RI7H_Z7.ec`` /
    ``igeo7/_bits.py`` behind the canonical ``Indexing`` protocol interface.
    All operations are pure stdlib with no runtime dependencies.
    """

    # DGGAL supports Z7 levels 0..19 (pydggal getMaxDGGRSZoneLevel() == 19).
    # The 64-bit packing has 20 direction-digit slots (_MAX_RES), so level 20 is
    # *representable* but is NOT a valid DGGAL zone: its geometric conversion
    # (to7H, RI7H_Z7.ec:348-353 "does not support level 20 zones") yields
    # nullZone, and getZoneFromWGS84Centroid(20, ...) returns nullZone. So the
    # highest resolution a zone can actually take is 19, distinct from the
    # packing width _MAX_RES used by resolution()/encode()/from_text().
    max_resolution: int = 19

    # ------------------------------------------------------------------
    # Bit-level primitives (same arithmetic as igeo7/_bits.py)
    # ------------------------------------------------------------------

    def _get_dir(self, bits: int, res: int) -> int:
        """Extract the 3-bit direction digit at resolution *res* (1-based)."""
        return (bits >> _bit_offset(res)) & 0x7

    def _set_dir(self, bits: int, res: int, d: int) -> int:
        """Return *bits* with the direction digit at *res* replaced by *d*."""
        off = _bit_offset(res)
        return (bits & ~(0x7 << off)) | (d << off)

    # ------------------------------------------------------------------
    # Indexing protocol — read operations
    # ------------------------------------------------------------------

    def base_cell(self, value: int) -> int:
        """Extract the 4-bit base cell (0-11) from bits 63-60."""
        return (value >> 60) & 0xF

    def resolution(self, value: int) -> int:
        """Return the resolution (0-20) by finding the first digit-7 slot.

        Mirrors the eC ``level`` property and ``igeo7/_bits.get_resolution``.
        """
        for r in range(1, _MAX_RES + 1):
            if self._get_dir(value, r) == 7:
                return r - 1
        return _MAX_RES

    def is_pentagon(self, value: int) -> bool:
        """True iff every direction digit in the path is 0 (centre child only).

        Resolution-0 cells (the 12 base cells) are always pentagons.
        """
        res = self.resolution(value)
        return res == 0 or all(self._get_dir(value, r) == 0 for r in range(1, res + 1))

    def num_children(self, value: int) -> int:
        """Return 6 for pentagons (one direction deleted), 7 for hexagons."""
        return 6 if self.is_pentagon(value) else 7

    def child_digits(self, value: int) -> list[int]:
        """Valid child direction digits: all 7 for a hexagon; the pentagon's
        deleted child (2 north / 5 south) omitted for a pentagon."""
        if self.is_pentagon(value):
            omit = _pentagon_omit(self.base_cell(value))
            return [d for d in range(7) if d != omit]
        return list(range(7))

    def parent(self, value: int) -> int | None:
        """Return the *congruent* parent cell (one resolution lower), or None at res 0.

        Obtained by dropping the last meaningful direction digit (set it back to
        7, the unused sentinel). This is a library-defined digit-path relation
        (exactly one parent), matching ``igeo7/_bits.get_parent``. It has NO
        DGGAL eC counterpart on purpose: DGGAL's Z7 hierarchy is the *geometric*,
        non-congruent one (2 parents / 13 children per interior cell). The eC
        ``from7H`` is the unrelated I7H<->Z7 address conversion, not a parent op.
        """
        res = self.resolution(value)
        return None if res == 0 else self._set_dir(value, res, 7)

    # ------------------------------------------------------------------
    # Indexing protocol — encode / decode
    # ------------------------------------------------------------------

    def encode(self, base: int, digits: list[int]) -> int:
        """Pack *(base, digits)* into a Z7 64-bit integer.

        Starts from ``_DEFAULT`` (all direction slots = 7), replaces the
        base-cell nibble, then writes each direction digit in order.
        """
        # Swap out the base-cell nibble (bits 63-60)
        bits = (_DEFAULT & ~(0xF << 60)) | (base << 60)
        # Write each direction digit into its slot
        for i, d in enumerate(digits):
            bits = self._set_dir(bits, i + 1, d)
        return bits

    def decode(self, value: int) -> tuple[int, list[int]]:
        """Unpack a Z7 integer into *(base_cell, [digit …])*.

        Inverse of ``encode``; resolution is inferred from the first 7 digit.
        """
        res = self.resolution(value)
        return (
            self.base_cell(value),
            [self._get_dir(value, r) for r in range(1, res + 1)],
        )

    # ------------------------------------------------------------------
    # Indexing protocol — text id (mirrors igeo7/zone.py)
    # ------------------------------------------------------------------

    def to_text(self, value: int) -> str:
        """Serialise *value* to its canonical Z7 text id.

        Format: two-digit zero-padded base cell, followed by one ASCII digit
        per direction digit. Matches the eC ``getTextID`` and
        ``OracleZone.text_id``.
        """
        base, digits = self.decode(value)
        return f"{base:02d}" + "".join(str(d) for d in digits)

    def from_text(self, text: str) -> int:
        """Parse a Z7 text id and return the 64-bit integer.

        Validates format, digit range, and the pentagon-canonical rule
        (a pentagon's deleted direction digit must never appear as a child
        digit along the all-zero prefix). Raises ``InvalidZoneError`` for
        any malformed or non-canonical input.
        """
        if len(text) < 2 or not (text.isascii() and text.isdigit()):
            raise InvalidZoneError(f"bad Z7 text id {text!r}")
        base = int(text[:2])
        digits = [int(c) for c in text[2:]]
        if len(digits) > _MAX_RES:
            raise InvalidZoneError(f"bad Z7 text id {text!r}: too many digits")
        if base > 11 or any(d > 6 for d in digits):
            raise InvalidZoneError(f"bad Z7 text id {text!r}")
        self._check_canonical(base, digits)
        return self.encode(base, digits)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _check_canonical(self, base: int, digits: list[int]) -> None:
        """Raise ``InvalidZoneError`` if *digits* contains a non-canonical
        pentagon-child digit.

        A cell's path is "pentagon-mode" while every direction digit so far
        is 0 (the centre-child path).  The first non-zero digit must NOT be
        the hemisphere-deleted digit (2 for north bases, 5 for south).
        Once a non-zero digit appears the cell is a hexagon and no further
        check is needed.

        Matches the eC ``fromTextID`` guard and ``igeo7/zone.py _check_canonical``.
        """
        omit = _pentagon_omit(base)
        for d in digits:
            if d == 0:
                continue   # still on the all-zero (pentagon) prefix
            if d == omit:
                raise InvalidZoneError(
                    f"non-canonical pentagon child: digit {d} invalid under base {base}"
                )
            break          # first non-zero non-omit digit → hexagon, done
