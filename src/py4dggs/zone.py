# src/py4dggs/zone.py
"""Immutable, hashable, grid-bound cell. Geometry delegates to its Grid, so a
Zone carries enough to answer any query about itself."""
from __future__ import annotations
from py4dggs.types import GeoPoint

class Zone:
    __slots__ = ("_grid", "_value")

    def __init__(self, grid, value: int):
        object.__setattr__(self, "_grid", grid)
        object.__setattr__(self, "_value", int(value))

    def __setattr__(self, *_):
        raise AttributeError("Zone is immutable")

    @property
    def value(self) -> int: return self._value
    @property
    def text_id(self) -> str: return self._grid.indexing.to_text(self._value)
    @property
    def resolution(self) -> int: return self._grid.indexing.resolution(self._value)
    @property
    def is_pentagon(self) -> bool: return self._grid.indexing.is_pentagon(self._value)

    @property
    def centroid(self) -> GeoPoint: return self._grid.centroid(self._value)
    @property
    def vertices(self): return self._grid.vertices(self._value)
    @property
    def neighbors(self):
        return tuple(Zone(self._grid, v) for v in self._grid.neighbors(self._value))

    @property
    def parent(self):
        """The primary parent (parents[0]), or None at the root. For congruent
        grids this is the single parent; for I3H it is DGGAL's ``parent0``."""
        ps = self._grid.parents(self._value)
        return Zone(self._grid, ps[0]) if ps else None

    @property
    def parents(self):
        """All parents (1 for congruent grids; 1 or 3 for the non-congruent I3H)."""
        return tuple(Zone(self._grid, v) for v in self._grid.parents(self._value))

    @property
    def children(self):
        return tuple(Zone(self._grid, v) for v in self._grid.children(self._value))

    @property
    def centroid_parent(self):
        """The parent whose refinement centre this cell is (I3H), or the single
        parent (congruent). None at the root."""
        cp = self._grid.centroid_parent(self._value)
        return None if cp is None else Zone(self._grid, cp)

    @property
    def is_centroid_child(self) -> bool:
        """Whether this cell has a single parent (always True for congruent grids)."""
        return self._grid.is_centroid_child(self._value)

    def is_immediate_child_of(self, other) -> bool:
        return any(p == other for p in self.parents)

    def sub_zones(self, relative_depth: int):
        return tuple(Zone(self._grid, v) for v in self._grid.sub_zones(self._value, relative_depth))

    def count_sub_zones(self, relative_depth: int) -> int:
        return self._grid.count_sub_zones(self._value, relative_depth)

    def first_sub_zone(self, relative_depth: int):
        return Zone(self._grid, self._grid.first_sub_zone(self._value, relative_depth))

    def sub_zone_index(self, sub_zone) -> int:
        """Index of ``sub_zone`` within this zone's sub-zone order, or -1.

        The grid check is not redundant: ISEA3H/IVEA3H/RTEA3H share an identical
        I3H packing, so comparing raw ints alone cannot tell a genuine sub-zone
        from the same value on a different grid, and the foreign zone would get
        a plausible index. ``Zone.__eq__`` already checks grid identity; this
        keeps the two consistent.
        """
        if sub_zone._grid is not self._grid:
            return -1
        return self._grid.sub_zone_index(self._value, sub_zone.value)

    def sub_zone_at_index(self, relative_depth: int, index: int):
        return Zone(self._grid, self._grid.sub_zone_at_index(self._value, relative_depth, index))

    def is_ancestor_of(self, other) -> bool:
        """Whether ``self`` is an ancestor of ``other`` (walks the parent DAG up)."""
        seen, frontier = set(), {other}
        while frontier:
            nxt = set()
            for z in frontier:
                for p in z.parents:
                    if p == self:
                        return True
                    if p not in seen:
                        seen.add(p); nxt.add(p)
            frontier = nxt
        return False

    def is_sibling_of(self, other) -> bool:
        """Whether ``self`` and ``other`` share a parent (and a resolution)."""
        if self.resolution != other.resolution:
            return False
        return bool({p.value for p in self.parents} & {p.value for p in other.parents})

    def disk(self, k: int = 1):
        seen = {self}; frontier = {self}
        for _ in range(k):
            nxt = set()
            for z in frontier: nxt |= set(z.neighbors)
            nxt -= seen; seen |= nxt; frontier = nxt
        return frozenset(seen)

    def __eq__(self, o): return isinstance(o, Zone) and o._value == self._value and o._grid is self._grid
    def __hash__(self): return hash((id(self._grid), self._value))
    def __repr__(self): return f"Zone({self.text_id!r})"
