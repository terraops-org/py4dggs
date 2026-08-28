# src/py4dggs/grid.py
"""A Grid is one DGGRS = Projection + Topology + Indexing. Operations needing
all three (notably the grid-agnostic geographic edge-crossing k-ring) live
here."""
from __future__ import annotations
import math
from dataclasses import dataclass
from py4dggs.interfaces import Projection, Topology, Indexing
from py4dggs.types import GeoPoint, GridConfig, InvalidZoneError

def _ll_to_xyz(lat, lon):
    la, lo = math.radians(lat), math.radians(lon); cla = math.cos(la)
    return [cla * math.cos(lo), cla * math.sin(lo), math.sin(la)]

def _xyz_to_ll(v):
    n = math.sqrt(sum(k * k for k in v)); x, y, z = (k / n for k in v)
    return math.degrees(math.asin(max(-1.0, min(1.0, z)))), math.degrees(math.atan2(y, x))


# DGGAL's geometry for a zone that has none: getZoneWGS84Centroid/Vertices on a
# nullZone give the zero point and six zero vertices (six regardless of whether
# the address is a pentagon path — verified against pydggal on all three
# aperture-7 grids). See `Topology.is_null_geometry` in interfaces.py.
_NULL_GEOPOINT = GeoPoint(0.0, 0.0)
_NULL_VERTICES = tuple(GeoPoint(0.0, 0.0) for _ in range(6))


@dataclass(frozen=True)
class Grid:
    projection: Projection
    topology: Topology
    indexing: Indexing
    config: GridConfig
    name: str

    def __post_init__(self):
        object.__setattr__(self, "_geom", self.projection.build_geometry(self.config))

    # --- indexing entry points ---
    def zone_from_geo(self, lat: float, lon: float, res: int):
        if not 0 <= res <= self.indexing.max_resolution:
            raise InvalidZoneError(f"resolution {res} out of range 0..{self.indexing.max_resolution}")
        p = self.projection.forward(self._geom, lat, lon)
        base, digits = self.topology.quantize(self._geom, p, res)
        from py4dggs.zone import Zone
        return Zone(self, self.indexing.encode(base, digits))

    def zone_from_text(self, text: str):
        """Construct a Zone from its canonical Z7 text id (validated via the indexing)."""
        from py4dggs.zone import Zone
        return Zone(self, self.indexing.from_text(text))

    # --- geometry: planar (topology) -> geographic (projection) ---
    def _is_null_geometry(self, base: int, digits: list[int]) -> bool:
        """Whether the topology declares this address geometry-less (DGGAL's
        `nullZone`). Optional capability — absent means "always has geometry",
        which is right for every topology that has no degenerate level."""
        fn = getattr(self.topology, "is_null_geometry", None)
        return fn is not None and fn(self._geom, base, digits)

    def centroid(self, value: int) -> GeoPoint:
        base, digits = self.indexing.decode(value)
        if self._is_null_geometry(base, digits):
            return _NULL_GEOPOINT
        p = self.topology.planar_centroid(self._geom, base, digits)
        return self.projection.inverse(self._geom, p)

    def vertices(self, value: int) -> tuple[GeoPoint, ...]:
        base, digits = self.indexing.decode(value)
        if self._is_null_geometry(base, digits):
            return _NULL_VERTICES
        ps = self.topology.planar_vertices(self._geom, base, digits)
        return tuple(self.projection.inverse(self._geom, p) for p in ps)

    # --- neighbours: exact topology override if provided, else grid-agnostic k-ring ---
    def neighbors(self, value: int) -> tuple[int, ...]:
        base, digits = self.indexing.decode(value)
        # A zone with no geometry has no neighbours: pydggal's getZoneNeighbors
        # returns [] for a Z7 level-20 zone. Checked before the k-ring below,
        # which would otherwise reflect the degenerate (0,0) vertices back
        # through quantize and raise a bare IndexError (it did so even before
        # is_null_geometry existed, when the vertices were merely meaningless
        # rather than zero).
        if self._is_null_geometry(base, digits):
            return ()
        # A Topology MAY expose an exact `neighbors` (aperture-3 does — it ports
        # DGGAL's exact adjacency); prefer it. Aperture-7 has none, so it falls
        # through to the grid-agnostic edge k-ring below (proven exact for 7H).
        nb_fn = getattr(self.topology, "neighbors", None)
        if nb_fn is not None:
            return tuple(self.indexing.encode(b, d) for (b, d) in nb_fn(self._geom, base, digits))
        # --- grid-agnostic geographic edge-crossing k-ring ---
        res = self.indexing.resolution(value)
        c = self.centroid(value); verts = self.vertices(value)
        c3 = _ll_to_xyz(c.lat, c.lon); n = len(verts)
        origin, out, seen = value, [], set()
        for i in range(n):
            v0 = _ll_to_xyz(*verts[i]); v1 = _ll_to_xyz(*verts[(i + 1) % n])
            t = [v0[k] + v1[k] - c3[k] for k in range(3)]
            lat, lon = _xyz_to_ll(t)
            p = self.projection.forward(self._geom, lat, lon)
            base, digits = self.topology.quantize(self._geom, p, res)
            nv = self.indexing.encode(base, digits)
            if nv == origin or nv in seen:
                continue
            seen.add(nv); out.append(nv)
        return tuple(out)

    # --- hierarchy: exact topology override if provided, else congruent digit-path ---
    def parents(self, value: int) -> tuple[int, ...]:
        fn = getattr(self.topology, "parents", None)
        if fn is not None:
            base, digits = self.indexing.decode(value)
            return tuple(self.indexing.encode(b, d) for (b, d) in fn(self._geom, base, digits))
        p = self.indexing.parent(value)
        return () if p is None else (p,)

    def children(self, value: int) -> tuple[int, ...]:
        fn = getattr(self.topology, "children", None)
        base, digits = self.indexing.decode(value)
        if fn is not None:
            return tuple(self.indexing.encode(b, d) for (b, d) in fn(self._geom, base, digits))
        # congruent digit-path default (Z7): a max-resolution zone has no
        # children -- you cannot refine past max_resolution (a further digit
        # does not fit the packing). Mirrors parent() -> () at resolution 0.
        if self.indexing.resolution(value) >= self.indexing.max_resolution:
            return ()
        return tuple(self.indexing.encode(base, digits + [d]) for d in self.indexing.child_digits(value))

    def centroid_parent(self, value: int) -> int | None:
        fn = getattr(self.topology, "centroid_parent", None)
        if fn is not None:
            base, digits = self.indexing.decode(value)
            r = fn(self._geom, base, digits)
            return None if r is None else self.indexing.encode(*r)
        return self.indexing.parent(value)  # congruent: the one parent IS the centroid parent

    def is_centroid_child(self, value: int) -> bool:
        fn = getattr(self.topology, "is_centroid_child", None)
        if fn is not None:
            base, digits = self.indexing.decode(value)
            return fn(self._geom, base, digits)
        return True  # congruent: every cell has exactly one parent

    # --- sub-zones: required topology override (no congruent-digit default exists) ---
    #: Ceiling on how many sub-zones :meth:`sub_zones` will materialise into a
    #: tuple. `count_sub_zones` is a cheap closed form, but `sub_zones` builds
    #: the whole sequence, and `relative_depth` was validated only against
    #: `max_resolution` -- never against result *cardinality*. On an I3H grid at
    #: resolution 0, depth 33 counts 4,632,550,579,746,406 sub-zones and passed
    #: every guard, so a caller-supplied depth (the README's DGGS-as-tile-store
    #: pattern) was a one-line denial of service. 4 million keeps every
    #: realistic tiling request working -- an I3H cell has 6 sub-zones at depth
    #: 1 and ~4.6M at depth 16 -- while refusing the unservable ones outright.
    #: Callers who want the count without the cost still call `count_sub_zones`.
    MAX_MATERIALISED_SUB_ZONES = 4_000_000

    def _sub_zone_fn(self, name: str, value: int, relative_depth: int):
        """Shared guard for count_sub_zones/first_sub_zone/sub_zones: resolve
        the topology override (or raise), then reject the two invalid-input
        cases that would otherwise reach the depth-dependent geometric
        generators, which aren't designed for them and misbehave silently
        (a negative depth turns `3**((depth-1)//2)`-style exponents into
        floats; an over-deep depth produces a zone past `max_resolution`
        with no error, since only `zone_from_geo` -- not the sub-zone path --
        validates against it)."""
        fn = getattr(self.topology, name, None)
        if fn is None:
            raise NotImplementedError("this grid has no sub-zone order")
        if relative_depth < 0:
            raise ValueError(f"relative_depth must be >= 0, got {relative_depth}")
        target_res = self.indexing.resolution(value) + relative_depth
        if target_res > self.indexing.max_resolution:
            raise InvalidZoneError(
                f"relative_depth {relative_depth} would reach resolution {target_res}, "
                f"past max_resolution {self.indexing.max_resolution}"
            )
        return fn

    def count_sub_zones(self, value: int, relative_depth: int) -> int:
        fn = self._sub_zone_fn("count_sub_zones", value, relative_depth)
        if relative_depth == 0:
            # I3H's own override treats depth 0 as "the zone itself" (eC
            # getI3HSubZoneCentroids, I3HSubZones.ec:1787-1790: "if(rDepth > 0)
            # ... else centroids[0] = zone.centroid" -- the base DGGRS class's
            # getSubZoneCRSCentroids is a null stub, not a generic dispatcher).
            # We hoist that semantic to the Grid level (deliberate
            # generalization, not a claim about where DGGAL implements it --
            # "depth 0 = same level = the zone itself" is universal to any
            # DGGS) because it's true regardless of grid/cell-class (verified
            # against pydggal's getSubZones for interior hex, non-polar
            # pentagon and polar pentagon -- all three return exactly [value]
            # at depth 0, even though pydggal's own raw per-grid
            # getFirstSubZone(zone, 0) is NOT self-consistent for the
            # pentagon/pole cases, an oracle quirk this generic short-circuit
            # deliberately does not replicate). Handled here, above the topology
            # call, so it never reaches the depth-dependent geometric generators
            # (which are not designed for depth 0 and would misbehave).
            return 1
        base, digits = self.indexing.decode(value)
        return fn(self._geom, base, digits, relative_depth)

    def first_sub_zone(self, value: int, relative_depth: int) -> int:
        fn = self._sub_zone_fn("first_sub_zone", value, relative_depth)
        if relative_depth == 0:
            return value  # see count_sub_zones' depth-0 note
        base, digits = self.indexing.decode(value)
        b, d = fn(self._geom, base, digits, relative_depth)
        return self.indexing.encode(b, d)

    def sub_zones(self, value: int, relative_depth: int) -> tuple:
        fn = self._sub_zone_fn("sub_zones", value, relative_depth)
        if relative_depth == 0:
            return (value,)  # see count_sub_zones' depth-0 note
        n = self.count_sub_zones(value, relative_depth)
        if n > self.MAX_MATERIALISED_SUB_ZONES:
            raise InvalidZoneError(
                f"relative_depth {relative_depth} yields {n} sub-zones, past the "
                f"{self.MAX_MATERIALISED_SUB_ZONES} materialisation limit; use "
                f"count_sub_zones() for the count, or refine in smaller steps"
            )
        base, digits = self.indexing.decode(value)
        return tuple(self.indexing.encode(b, d) for (b, d) in fn(self._geom, base, digits, relative_depth))

    def sub_zone_index(self, value: int, sub_zone_value: int) -> int:
        """Generic (dggrs.ec:115-131): build the full ordered list, find the
        position. Not the generators' internal index>=0 fast-forward path
        (explicitly out of scope -- see the A3 design spec)."""
        # caller must know the sub-zone's relative depth; derive it from level
        sub_level = self.indexing.resolution(sub_zone_value)
        parent_level = self.indexing.resolution(value)
        relative_depth = sub_level - parent_level
        if relative_depth == 0:
            # depth 0 = "the zone itself": sub_zones(v, 0) == (v,),
            # count_sub_zones(v, 0) == 1 and sub_zone_at_index(v, 0, 0) == v, so
            # a zone IS its own sub-zone at index 0. Folding this into the
            # `< 0` rejection reported a real sub-zone as not-a-descendant.
            return 0 if sub_zone_value == value else -1
        if relative_depth < 0:
            return -1
        subs = self.sub_zones(value, relative_depth)
        try:
            return subs.index(sub_zone_value)
        except ValueError:
            return -1

    def sub_zone_at_index(self, value: int, relative_depth: int, index: int) -> int:
        """Generic (dggrs.ec:133-149)."""
        if index < 0 or index >= self.count_sub_zones(value, relative_depth):
            raise IndexError(f"index {index} out of range for relative_depth {relative_depth}")
        if index == 0:
            return self.first_sub_zone(value, relative_depth)
        return self.sub_zones(value, relative_depth)[index]
