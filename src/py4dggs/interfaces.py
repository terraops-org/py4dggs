# src/py4dggs/interfaces.py
"""The three pluggable parts a Grid composes. A concrete DGGRS = one
Projection + one Topology + one Indexing. New grids add new parts; they never
edit existing ones.

`geom` is the opaque, precomputed geometry a Projection builds from a
GridConfig (icosahedron vertices, authalic coefficients, etc.); it is passed
back into the projection/topology calls so they stay pure and config-driven.
"""
from __future__ import annotations
from typing import Any, Protocol, runtime_checkable
from py4dggs.types import GeoPoint, PlanarPoint, GridConfig


@runtime_checkable
class Projection(Protocol):
    """Geographic sphere/ellipsoid <-> planar face coordinates."""
    def build_geometry(self, config: GridConfig) -> Any: ...
    def forward(self, geom: Any, lat: float, lon: float) -> PlanarPoint: ...
    def inverse(self, geom: Any, p: PlanarPoint) -> GeoPoint: ...


@runtime_checkable
class Topology(Protocol):
    """Planar coords <-> (base cell, direction digits), plus PLANAR cell
    geometry. Carries the aperture and cell shape. All geometry here is in the
    projection's planar space; the Grid lifts it to geographic via the
    Projection.

    Note: ``@runtime_checkable`` isinstance() does NOT verify ``aperture``
    (only methods are checked) — every concrete Topology MUST set ``aperture``
    as a real class attribute."""
    aperture: int
    def quantize(self, geom: Any, p: PlanarPoint, res: int) -> tuple[int, list[int]]: ...
    def planar_centroid(self, geom: Any, base: int, digits: list[int]) -> PlanarPoint: ...
    def planar_vertices(self, geom: Any, base: int, digits: list[int]) -> list[PlanarPoint]: ...
    # OPTIONAL: a Topology MAY expose exact topological neighbours. When present,
    # Grid.neighbors prefers it over the grid-agnostic edge k-ring. hex_a3 defines
    # it (exact aperture-3 adjacency); hex_a7 omits it (the edge k-ring is exact
    # for aperture-7). Signature when present:
    #   def neighbors(self, geom: Any, base: int, digits: list[int]) -> list[tuple[int, list[int]]]: ...
    #
    # OPTIONAL (A2): a Topology MAY expose an exact NON-congruent geometric
    # hierarchy. When present, Grid.parents/children/centroid_parent/is_centroid_child
    # prefer these over the congruent digit-path default (which uses Indexing.parent/
    # child_digits). hex_a3 defines them (exact I3H hierarchy); hex_a7 omits them
    # (the congruent Z7 digit hierarchy is used). Signatures when present:
    #   def parents(self, geom, base, digits) -> list[tuple[int, list[int]]]: ...        # 1 or 3
    #   def children(self, geom, base, digits) -> list[tuple[int, list[int]]]: ...        # 6 or 7
    #   def centroid_parent(self, geom, base, digits) -> tuple[int, list[int]] | None: ...
    #   def is_centroid_child(self, geom, base, digits) -> bool: ...
    #
    # OPTIONAL (A3): a Topology MAY expose an ordered sub-zone enumeration (a
    # refinement to `relative_depth` levels below the cell, e.g. depth=2 skips
    # a resolution). When present, Grid.count_sub_zones/first_sub_zone/
    # sub_zones use these directly (there is no congruent-digit default -- a
    # grid either has a sub-zone order or Grid raises NotImplementedError).
    # hex_a3 defines them (I3H's exact sub-zone generators); hex_a7 omits them.
    # Signatures when present:
    #   def count_sub_zones(self, geom, base, digits, relative_depth) -> int: ...
    #   def first_sub_zone(self, geom, base, digits, relative_depth) -> tuple[int, list[int]]: ...
    #   def sub_zones(self, geom, base, digits, relative_depth) -> list[tuple[int, list[int]]]: ...
    #
    # OPTIONAL: a Topology MAY declare that a cell's geometry is DEGENERATE —
    # that the address is representable in the packing but has no geometry,
    # DGGAL's `nullZone` outcome. When present, Grid.centroid/vertices return
    # DGGAL's null geometry (lat/lon 0,0; six zero vertices) instead of calling
    # planar_centroid/planar_vertices, mirroring the eC. hex_a7 defines it (Z7
    # level 20: representable in the 20-slot packing, but `to7H`,
    # RI7H_Z7.ec:348-353, "does not support level 20 zones" and yields nullZone).
    # hex_a3 omits it: the I3H nullZone sentinel HAS meaningless-but-nonzero
    # geometry in DGGAL, which this port reproduces, so overriding it would be
    # the unfaithful choice. Signature when present:
    #   def is_null_geometry(self, geom, base, digits) -> bool: ...
    #
    # Grid pre-validates before calling any of the three above (a Topology's
    # sub-zone functions never see these cases): relative_depth == 0 is
    # intercepted at the Grid level ("the zone itself", universal to any DGGS,
    # not delegated since the geometric generators aren't designed for depth
    # 0); relative_depth < 0 raises ValueError; and a relative_depth that would
    # push the result past Indexing.max_resolution raises InvalidZoneError
    # (mirroring the same bound Grid.zone_from_geo enforces on the input side).
    # See Grid._sub_zone_fn.


@runtime_checkable
class Indexing(Protocol):
    """(base, digits) <-> integer id <-> text id.

    Note: ``@runtime_checkable`` isinstance() does NOT verify ``max_resolution``
    (only methods are checked) — every concrete Indexing MUST set
    ``max_resolution`` as a real class attribute."""
    max_resolution: int
    def encode(self, base: int, digits: list[int]) -> int: ...
    def decode(self, value: int) -> tuple[int, list[int]]: ...
    def resolution(self, value: int) -> int: ...
    def base_cell(self, value: int) -> int: ...
    def is_pentagon(self, value: int) -> bool: ...
    def to_text(self, value: int) -> str: ...
    def from_text(self, text: str) -> int: ...
    # OPTIONAL: an Indexing whose grid has a CONGRUENT digit-path hierarchy MAY
    # expose it here; Grid falls back to these when the Topology supplies no
    # geometric hierarchy override. Z7Indexing defines them; I3HIndexing does
    # NOT — the I3H hierarchy is non-congruent and geometric, so it lives on
    # hex_a3 (see the Topology hierarchy note above) and is never reached
    # through the Indexing. These are deliberately kept OUT of the Protocol
    # body: as required members they would make `isinstance(I3HIndexing(),
    # Indexing)` False, i.e. the declared contract would be one that a shipped,
    # correct implementation fails. Signatures when present:
    #   def parent(self, value: int) -> int | None: ...
    #   def num_children(self, value: int) -> int: ...
    #   def child_digits(self, value: int) -> list[int]: ...   # omits the deleted pentagon child
