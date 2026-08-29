# py4dggs

*Born 7th July 2026, Aveiro, Portugal @ 18:50 WET*

Canonical, readable, pure-Python multi-grid DGGS reference library (package name `py4dggs`, therefore Python for DGGS).

This library provides a clean-room Python implementation of Discrete Global Grid Systems (DGGS),
designed for clarity, correctness, and zero runtime dependencies.

Its reference point is **[DGGAL](https://dggal.org)**, the Discrete Global Grid Abstraction Library
from Ecere ([source](https://github.com/ecere/dggal)): the canonical engine, written in the eC
language, and the ground truth for every grid implemented here. Nothing is trusted in isolation, so
every grid, IGEO7 included, is verified against **[pydggal](https://github.com/ecere/pydggal)**,
DGGAL's official Python binding to the `libdggal` C library, installed from PyPI as the
[`dggal`](https://pypi.org/project/dggal/) package. That is a dev dependency only, so the runtime
stays dependency-free (and there are no sibling-repo dependencies either). See "Verification
against pydggal" below for what is checked and how.

## Quickstart

Python 3.12 or newer, and zero runtime dependencies, so nothing else gets
pulled in:

```bash
uv add py4dggs           # or: pip install py4dggs
```

To work on the library itself, or to run its test suite, install from a clone
instead:

```bash
git clone https://github.com/terraops-org/py4dggs.git
cd py4dggs
uv sync                  # creates .venv/ and installs the dev extras
```

`uv sync` also brings in the dev dependencies: `pytest`, and `dggal` (pydggal,
the DGGAL engine's own Python binding), which the suite uses as its correctness
oracle. If you would rather not use `uv`:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .         # add: pip install pytest dggal==0.0.6   to run the tests
```

Confirm the install, and that this machine reproduces the verified behaviour:

```bash
uv run python -c "import py4dggs; print(py4dggs.__all__)"
uv run pytest -q                          # the full suite, ~4300 tests in ~11s
uv run python examples/01_first_zone.py   # a first runnable example
```

From there, `uv run python` drops you into a REPL with the package importable,
or activate `.venv` and use plain `python`.

Runnable scripts live in **`examples/`** (`uv run python examples/01_first_zone.py`); see
`examples/README.md` for the list. They and every code block in this file are executed by the test
suite, so they cannot silently go stale.

See **`TUTORIAL.md`** for a longer, task-oriented walkthrough with worked examples (choosing a
grid, neighbours, hierarchy, sub-zones/tiling, error handling), and `ARCHITECTURE.md` for how the
pieces below fit together and how to read the source.

```python
from py4dggs import IGEO7, ISEA3H, get_grid

# Encode a geographic point into a zone at a given resolution
lisbon = IGEO7.zone_from_geo(lat=38.7223, lon=-9.1393, res=7)
print(lisbon.text_id)      # "006415654" -- IGEO7's canonical Z7 digit id
print(lisbon.value)        # 940643742434459647 -- packed int (== pydggal's DGGRSZone int)
print(lisbon.resolution)   # 7
print(lisbon.is_pentagon)  # False

# Geometry (always WGS84 GeoPoint(lat, lon))
print(lisbon.centroid)     # GeoPoint(lat=38.74412563021016, lon=-9.138855236620287)
print(len(lisbon.vertices))  # 6 (or 5 for a pentagon)

# Round-trip through a zone's canonical text id
same = IGEO7.zone_from_text(lisbon.text_id)
assert same == lisbon

# Topology: neighbours, hierarchy
print(len(lisbon.neighbors))       # 6 (or 5 for a pentagon)
parent = lisbon.parent             # one level coarser (None at resolution 0)
children = lisbon.children         # one level finer, 7 of them (6 for a pentagon)

# A different grid family: ISEA3H (aperture-3, rhombic addressing)
rome = ISEA3H.zone_from_geo(lat=41.9028, lon=12.4964, res=6)
print(rome.text_id)                # "D4-83-A" -- rhombic id, NOT a digit path
print(rome.parents)                # (Zone('C4-10-C'),) -- 1 or 3 Zones (non-congruent hierarchy)

# Sub-zones -- the DGGS-as-storage primitive (I3H grids only, see "Sub-zones" below)
depth = 2
tile_size = rome.count_sub_zones(depth)          # 13 -- fixed array length for this "tile"
slots = rome.sub_zones(depth)                    # ordered tuple[Zone, ...], len == tile_size
assert rome.sub_zone_index(slots[0]) == 0
assert rome.sub_zone_at_index(depth, 0) == slots[0]

# Grids can also be looked up by name (useful when the grid is chosen at runtime)
assert get_grid("ISEA3H") is ISEA3H
```

(All values above are actual output from running this exact snippet - not illustrative placeholders.)

All six registered grids: `IGEO7`, `IVEA7H`, `RTEA7H` (aperture-7, Z7 digit addressing) and
`ISEA3H`, `IVEA3H`, `RTEA3H` (aperture-3, I3H rhombic addressing) - see "Grids" below for what
distinguishes each. Every grid exposes the same `Zone`/`Grid` API shown above; what differs is
which pieces are exact vs. approximated (the "Verification against pydggal" section) and which
optional capabilities a grid's `Topology` provides (exact neighbours, geometric hierarchy,
sub-zones - all six now have all three, per "Verification against pydggal" and "Sub-zones" below).

## Grids

`py4dggs` is multi-grid. Each grid composes a `Projection` + `Topology` + `Indexing`:

- **`IGEO7`** - ISEA projection + aperture-7 hex topology + Z7 indexing (the IGEO7/ISEA7H grid).
- **`IVEA7H`** - the **IVEA** projection variant with the *same* aperture-7 hex topology + Z7
  indexing. Adding it changed **zero** of the topology/indexing/grid layers - only the projection
  (see below). Get either via `from py4dggs import IGEO7, IVEA7H` or `get_grid("IVEA7H")`.
- **`RTEA7H`** - the **RTEA** projection variant, the third `VGCRadialVertex` case, again on the
  *same* aperture-7 hex topology + Z7 indexing. Like `IVEA7H`, adding it changed only the
  projection (`variant_consts("rtea")` in `icovertex.py` + a thin `rtea.py` + one registry line) -
  zero changes to topology/indexing/grid. Get it via `from py4dggs import RTEA7H` or
  `get_grid("RTEA7H")`.
- **`ISEA3H`** - ISEA projection + **aperture-3** rhombic topology + packed **I3H** rhombic indexing.
  This is a second-topology grid, DGGAL-verified (against `pydggal` ISEA3H) for quantization,
  centroids, vertices, **exact topological neighbours** (sub-project A slice **A0**), the
  **canonical rhombic text-id** `"C2-23-C"` (slice **A1**), the **non-congruent geometric
  hierarchy** - `parents` (1 or 3), geometric `children`, `centroid_parent`/`is_centroid_child`
  (slice **A2**) - and **sub-zones** (OGC descendants-at-depth; slice **A3**, see "Sub-zones"
  below). Get it via `from py4dggs import ISEA3H` or `get_grid("ISEA3H")`.
- **`IVEA3H`** - the **IVEA** projection variant on the *same* aperture-3 rhombic topology + I3H
  indexing as `ISEA3H`. A pure combination - one registry line, zero new topology/indexing code
  (`HexAperture3Topology`/`I3HIndexing` are already projection-agnostic). DGGAL-verified
  (against `pydggal` IVEA3H) for quantization, centroids, vertices, exact neighbours (A0), the
  canonical rhombic text-id (A1), the non-congruent geometric hierarchy (A2), and sub-zones (A3).
  Get it via `from py4dggs import IVEA3H` or `get_grid("IVEA3H")`.
- **`RTEA3H`** - the **RTEA** projection variant on the *same* aperture-3 rhombic topology + I3H
  indexing as `ISEA3H`/`IVEA3H`. Again a pure combination - one registry line, zero new
  topology/indexing code. DGGAL-verified (against `pydggal` RTEA3H) for quantization,
  centroids, vertices, exact neighbours (A0), the canonical rhombic text-id (A1), the
  non-congruent geometric hierarchy (A2), and sub-zones (A3).
  Get it via `from py4dggs import RTEA3H` or `get_grid("RTEA3H")`. `RTEA3H` is the final cell of the
  projection × aperture matrix - **ISEA/IVEA/RTEA × 7H/3H**, all six grids now registered, and
  all six now carry the full sub-project A surface (exact neighbours, text-id, geometric
  hierarchy, sub-zones).

The projection layer is a shared kernel (`projections/icovertex.py`) parameterized by a
`radial_vertex` variant, with thin per-variant classes (`isea.py`, `ivea.py`, `rtea.py`) -
mirroring DGGAL's `SliceAndDiceGreatCircleIcosahedralProjection` + `VGCRadialVertex {isea, ivea,
rtea}`. All three variants (IGEO7/IVEA7H/RTEA7H) now drop into this same kernel.

## Verification against pydggal (DGGAL ground truth)

Every grid in this librar, IGEO7 included, is verified directly against **pydggal**, the DGGAL
Python binding, i.e. the canonical engine itself. Install it as a dev dependency:
```
uv add --dev dggal    # dggal==0.0.6
```

**Golden tables are vendored, so nothing external is required.** The golden-table tests below
(`test_conformance.py` and its per-grid siblings) read frozen fixtures from `tests/golden/`, which
ships in this repo along with the `generate.py` that produced them. So `git clone && uv sync && uv
run pytest` runs the complete suitem, no sibling repo, no credential, no configuration. The tables
originate from [igeo7-spec](https://github.com/terraops-org/igeo7-spec) (private) and were vendored
2026-07-29; `tests/golden/PROVENANCE.md` records the source commit and the regeneration commands.
The path stays overridable per-grid with `DGGS_GOLDEN_TABLES[_<GRID>]=...` for an out-of-tree table
set.

**`IGEO7`** is verified two ways (the tables ship in-repo; the live-oracle half skips cleanly if `dggal` is absent):

- **Live differential fuzz** (`tests/test_isea7h_fuzz.py`) - forward / centroid / vertices /
  neighbours vs pydggal; hierarchy (parent/children) is checked structurally instead, since Z7's
  congruent digit-path hierarchy has no pydggal equivalent (pydggal's own `getZoneParent`/
  `getZoneChildren` return DGGAL's *geometric* Z7 hierarchy, a different, non-congruent relation).
- **Golden-table conformance** (`tests/test_conformance.py` → **581 passed, 7 xfailed**) -
  against `tests/golden/`. Override the path with `export DGGS_GOLDEN_TABLES=...`. The
  7 xfails are benign pole/vertex-0 (lon 11.2°) boundary tie-breaks.

**`IVEA7H`** is verified two ways (the tables ship in-repo; the live-oracle half skips cleanly if `dggal` is absent):

- **Live differential fuzz** (`tests/test_ivea7h_fuzz.py`) - forward / centroid / vertices /
  neighbours vs pydggal. Discrete outputs match exactly; float geometry agrees to ~1e-9°
  (pydggal is a C engine, so agreement is numerical, not bit-identical).
- **Golden-table conformance** (`tests/test_ivea7h_conformance.py` -> **553 passed, 13 xfailed**) against `tests/golden/ivea7h/` (regenerate with `uv run python tests/golden/generate.py IVEA7H_Z7`).
  Override the path with `export DGGS_GOLDEN_TABLES_IVEA7H=...`. The 13 xfails are exact-singular
  boundary tie-breaks at the pole and the vertex-0 meridian (lon 11.2°).

**`RTEA7H`** is verified the same two ways (the tables ship in-repo; the live-oracle half skips cleanly if `dggal` is absent):

- **Live differential fuzz** (`tests/test_rtea7h_fuzz.py`) - forward / centroid / vertices /
  neighbours vs pydggal. Discrete outputs match exactly; float geometry agrees to ~1e-9°.
- **Golden-table conformance** (`tests/test_rtea7h_conformance.py` -> **552 passed, 14 xfailed**) -
  against `tests/golden/rtea7h/` (regenerate with `uv run python tests/golden/generate.py RTEA7H_Z7`).
  Override the path with `export DGGS_GOLDEN_TABLES_RTEA7H=...`. 13 of the 14 xfails are the same
  exact-singularity class as `IVEA7H` (pole + vertex-0 meridian, lon 11.2°); the 14th is an
  ordinary point (London, res 3) whose RTEA-projected cell happens to straddle a 5x6-layout
  rhombus interruption seam, confirmed adjacent to pydggal's cell (both engines agree exactly on
  both cells' vertices/centroids) and consistent with the low, adjacent-only mismatch rate seen
  under broader random fuzzing. See `tests/test_rtea7h_conformance.py` for the full derivation.

**`ISEA3H`** (point-keyed geometry) is verified two ways (the tables ship in-repo; the live-oracle half skips cleanly if `dggal` is absent):

- **Live differential fuzz** (`tests/test_isea3h_fuzz.py`) - quantize / centroid / vertices /
  neighbours vs pydggal. Discrete outputs match exactly; float geometry agrees to ~1e-9°.
- **Golden-table conformance** (`tests/test_isea3h_conformance.py` -> **770 passed**) - against
  `tests/golden/isea3h/` (regenerate with `uv run python tests/golden/generate.py ISEA3H`). Override
  the path with `export DGGS_GOLDEN_TABLES_ISEA3H=...`. The verification lever is DGGAL-exact
  `uint64` int packing: our Zone value equals pydggal's `DGGRSZone` int. Includes value-keyed
  text-id/hierarchy/sub-zones tables (154 cases each) alongside the point-keyed geometry tables.

**`IVEA3H`** (point-keyed geometry) is verified the same two ways (both skip cleanly if `dggal` is
absent), using the same DGGAL-exact `uint64` int-packing lever as `ISEA3H`:

- **Live differential fuzz** (`tests/test_ivea3h_fuzz.py`) - quantize / centroid / vertices /
  neighbours vs pydggal. Discrete outputs match exactly; float geometry agrees to ~1e-9°.
- **Golden-table conformance** (`tests/test_ivea3h_conformance.py` -> **758 passed, 12 xfailed**) -
  against `tests/golden/ivea3h/` (regenerate with `uv run python tests/golden/generate.py IVEA3H`).
  Override the path with `export DGGS_GOLDEN_TABLES_IVEA3H=...`. The 12 xfails are exact-boundary
  tie-breaks at the vertex-0 meridian (lon 11.2°, even resolutions only) - the same singularity
  class `IVEA7H` documents, confirmed benign (adjacent cell; both engines agree exactly on that
  cell's own geometry). `ISEA3H` has zero xfails here since ISEA lacks IVEA's extra
  pole-longitude sensitivity. Like text-id/hierarchy, the sub-zones table (A3) is value-keyed, so
  it is immune to this boundary tie-break too.

**`RTEA3H`** (point-keyed geometry) is verified the same two ways (both skip cleanly if `dggal` is
absent), using the same DGGAL-exact `uint64` int-packing lever as `ISEA3H`/`IVEA3H`:

- **Live differential fuzz** (`tests/test_rtea3h_fuzz.py`) - quantize / centroid / vertices /
  neighbours vs pydggal. Quantize/centroid/vertices match exactly (0 tolerance); the neighbour (k-ring) check is deliberately tolerant, see the note below.
- **Golden-table conformance** (`tests/test_rtea3h_conformance.py` -> **748 passed, 22 xfailed**)
  against `tests/golden/rtea3h/` (regenerate with `uv run python tests/golden/generate.py RTEA3H`).
  Override the path with `export DGGS_GOLDEN_TABLES_RTEA3H=...`. The 22 xfails are exact-boundary
  tie-breaks at the vertex-0 meridian (lon 11.2°) at **every** resolution 0-10 (broader than
  `IVEA3H`'s even-only set, because RTEA's own vertex-assignment permutation warps the cells
  straddling this exact meridian differently), confirmed benign the same way (adjacent cell; both
  engines agree exactly on that cell's own geometry). The sub-zones table (A3), like text-id/hierarchy, is value-keyed and unaffected by this tie-break.

**A note on neighbours (`Grid.neighbors`):** by default the k-ring is a grid-agnostic *geometric* edge-crossing construction (reflect the centroid through each edge midpoint, re-quantize the
reflected point). This is **exact for every aperture-7 grid** (IGEO7/IVEA7H/RTEA7H) - verified 0
mismatch vs pydggal's exact topological neighbours. It is **approximate for the aperture-3 grids** (~0.01% of cells, at root-rhombus boundaries / interruption seams / pentagons / polar rows, where
the edge-reflection overshoots), so sub-project A slice **A0** gave the `hex_a3` topology an **exact `neighbors` override**m a faithful port of DGGAL's `I3HZone::getNeighbor`/`getNeighbors`.
`Grid.neighbors` prefers a topology's override when present, else falls back to the geometric
k-ring. **So all six grids now have exact neighbours vs pydggal** (int-set equality, any seed);
`test_{isea3h,ivea3h,rtea3h}_fuzz.py` assert strict equality. The geometric k-ring remains the
grid-agnostic fallback for any future topology that doesn't supply its own.

**A note on the hierarchy (`Zone.parents`/`.children`/`.centroid_parent`):** each grid returns its
*native* hierarchy. The **I3H** grids have only a geometric one - non-congruent, `parents` = 1 or 3,
`children` = 6 (pentagon) or 7 (hexagon), ported exactly from DGGAL's `getZoneParents`/
`getZoneChildren` (`hex_a3` topology override, verified vs pydggal, slice A2). The **Z7** grids
return their congruent *digit* hierarchy (1 parent / 7 children - append/drop a Z7 digit), matching
`igeo7-py`. Note this deliberately differs from DGGAL's *geometric* hierarchy for `ISEA7H_Z7`, which
is itself non-congruent (2 parents / 13 children, incl. cells shared from neighbouring base rhombi);
reproducing that geometric view for the Z7 grids would be a separate slice. `Grid.parents`/
`children`/`centroid_parent`/`is_centroid_child` prefer a topology override (I3H) and otherwise use
the congruent digit-path default (Z7).

## Sub-zones

Sub-project A slice **A3** adds OGC-style *descendants-at-depth*, the set of all cells `relative_depth`
levels below a zone (not just its immediate children), as an ordered, indexable sequence, for the
**ISEA3H**/**IVEA3H**/**RTEA3H** (I3H) grids. This is DGGAL's `getSubZones`/`countSubZones`/
`getFirstSubZone`/`getSubZoneIndex`/`getSubZoneAtIndex` family, ported faithfully from
`I3HSubZones.ec` and verified 0-mismatch against `pydggal` across all four I3H cell classes (interior hexagon, edge hexagon straddling a rhombus interruption, non-polar pentagon, polar pentagon) plus a reference-vector replay of the eC source's own worked examples.

**Grid coverage:** all five methods below work on the three I3H grids
(`ISEA3H`/`IVEA3H`/`RTEA3H`) and raise `NotImplementedError("this grid has no
sub-zone order")` on the three Z7 grids (`IGEO7`/`IVEA7H`/`RTEA7H`). Z7 sub-zones
are a **planned roadmap item**, not a permanent non-goal - see `TUTORIAL.md`.

**API** (mirrored on `Grid` and `Zone`):

- `Grid.count_sub_zones(value, relative_depth) -> int` / `Zone.count_sub_zones(relative_depth)` -
  how many descendants a zone has `relative_depth` levels down (closed-form: 7/13/37/... for a
  hexagon at depth 1/2/3, fewer for a pentagon).
- `Grid.sub_zones(value, relative_depth) -> tuple[int, ...]` / `Zone.sub_zones(relative_depth) ->
  tuple[Zone, ...]` - the full ordered list of descendant values/`Zone`s.
- `Grid.first_sub_zone(value, relative_depth) -> int` / `Zone.first_sub_zone(relative_depth) ->
  Zone` - the *centroid* descendant (index 0), computed directly without building the whole list.
- `Grid.sub_zone_index(value, sub_zone_value) -> int` / `Zone.sub_zone_index(sub_zone) -> int` -
  the position of a known descendant within its parent's ordered sub-zone list (`-1` if it is not
  a descendant at a valid depth). A zone is its own sub-zone at index **0** (`relative_depth == 0`),
  matching `sub_zones(v, 0) == (v,)` and pydggal's own `getSubZoneIndex(v, v)`; an *ancestor* gives
  `-1`. On `Zone`, a zone belonging to a **different grid** also gives `-1` - the I3H grids share an
  identical packing, so the grid identity, not just the int, decides.
- `Grid.sub_zone_at_index(value, relative_depth, index) -> int` / `Zone.sub_zone_at_index(
  relative_depth, index) -> Zone` - the descendant at a given position (`IndexError` if out of
  range).

Both `sub_zone_index` and `sub_zone_at_index` are the **generic** `dggrs.ec`-style implementation
(build the ordered list, then index into or search it); DGGAL's internal `index >= 0`
*fast-forward* short-circuit - computing an arbitrary index's cell directly, without materializing
the whole list is deliberately **out of scope** for this port (see the A3 design spec). Callers needing that performance characteristic for very deep sub-zone sets should be aware both methods are currently `O(count_sub_zones)`.

**Bound on materialisation:** `count_sub_zones` is a cheap closed form, but `sub_zones` builds the
whole tuple, so `relative_depth` is additionally checked against `Grid.MAX_MATERIALISED_SUB_ZONES`
(4,000,000) and raises `InvalidZoneError` past it. Without that bound a caller-supplied depth was
unbounded - depth 33 on a resolution-0 I3H cell counts 4.6x10^15 sub-zones and would try to build
them all. Ask `count_sub_zones` first when the depth comes from a request.

**Why this matters - DGGS-as-storage:** sub-zones are the mechanism that turns a DGGS into a fixed size **tile store**. Pick a coarse zone as a "tile" (e.g. an ISEA3H cell at resolution 6) and a `relative_depth` (e.g. 4): `count_sub_zones` gives you the tile's fixed array length up front, `sub_zones`/`sub_zone_index`/`sub_zone_at_index` give you a stable, deterministic mapping between
"array slot" and "fine-resolution cell" - exactly the raster-band/array-index model a format like

**Golden-table conformance:** like text-id/hierarchy (A1/A2), sub-zones are checked against a frozen,
value-keyed golden table (`tests/golden/{isea3h,ivea3h,rtea3h}/subzones.json` -
`{name, res, value, subZonesByDepth: {"1": [...], "2": [...], "3": [...]}}`, regenerate with
`python generate.py {ISEA3H,IVEA3H,RTEA3H}`) via `test_subzones` in each grid's
`tests/test_*3h_conformance.py`, pinning correctness even without a live `pydggal` install. Being
value-keyed, it is immune to the vertex-0-meridian forward-quantization tie-break the point-keyed
geometry tables xfail on (see above) `ISEA3H`/`IVEA3H`/`RTEA3H` all have **zero** `test_subzones`
xfails.

