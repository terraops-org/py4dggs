# Tutorial — using `py4dggs` to work with real cells

This is a task-oriented walkthrough, distinct from the other two docs:

- **`README.md`** — what's verified and the terse API reference.
- **`ARCHITECTURE.md`** — how the code is structured, for studying/extending it.
- **This file** — how to *do things* with it, with worked, runnable examples.

Every snippet below was actually run against this exact library version, the
printed output is real, not illustrative. Run any of them yourself with

`uv run python -c "..."` from the repo root (or a fresh `python` after
`pip install py4dggs`/`uv add py4dggs`).

## Install

```bash
uv add py4dggs           # or: pip install py4dggs
```

Both the PyPI distribution name and the import name are `py4dggs` (there's already an unrelated
`dggs` name on PyPI, and the established DGGAL Python wrapper is `dggs4py` — this project is the
pure-Python inverse of that, no C library required, hence the name).

**Not yet published to PyPI.** Until it is, install from a local clone instead:

```bash
git clone https://github.com/terraops-org/py4dggs.git
cd py4dggs && uv sync    # or: pip install -e .
```

Zero runtime dependencies — nothing else gets pulled in.

## Your first zone

The one thing every DGGS operation starts from: turn a geographic point into a
**zone** (a discrete cell) at a chosen **resolution** (finer resolution =
smaller cell).

```python
from py4dggs import IGEO7

tokyo = IGEO7.zone_from_geo(lat=35.6762, lon=139.6503, res=6)

print(tokyo.text_id)      # "05460005" -- IGEO7's canonical Z7 digit id
print(tokyo.value)        # 6449181054673616895 -- packed int (== pydggal's DGGRSZone int)
print(tokyo.resolution)   # 6
print(tokyo.is_pentagon)  # False
print(tokyo.centroid)     # GeoPoint(lat=35.71138423764992, lon=139.70703969781812)
print(len(tokyo.vertices))  # 6 (5 for a pentagon)
```

`tokyo.centroid` is the cell's *own* center (not your input lat/lon) — a
`GeoPoint(lat, lon)` in WGS84 degrees. `tokyo.vertices` gives you the cell's
boundary as a list of `GeoPoint`s (6 for a hexagon, 5 for one of the 12
pentagons every DGGS grid has).

## Converting between representations

A zone has three interchangeable forms — a `Zone` object, its canonical
**text id**, and its packed **integer value**. All three round-trip exactly:

```python
from py4dggs import IGEO7, Zone

tokyo = IGEO7.zone_from_geo(lat=35.6762, lon=139.6503, res=6)

# text id -> Zone
same = IGEO7.zone_from_text(tokyo.text_id)
print(same == tokyo)              # True

# int value -> Zone (Zone() takes the grid + the packed int directly)
same2 = Zone(IGEO7, tokyo.value)
print(same2 == tokyo)              # True
```

`Zone` is a plain immutable `(grid, value)` pair — cheap to construct, safe to
put in a `set`/use as a `dict` key (it's hashable), and comparing two `Zone`s
only makes sense if they came from the *same* grid.

## Choosing a grid

`py4dggs` registers **six** grids — every combination of 3 projections ×
2 topology/indexing families:

| | **Z7** (aperture-7, digit-addressed) | **I3H** (aperture-3, rhombic-addressed) |
|---|---|---|
| **ISEA** projection | `IGEO7` | `ISEA3H` |
| **IVEA** projection | `IVEA7H` | `IVEA3H` |
| **RTEA** projection | `RTEA7H` | `RTEA3H` |

**Which one should you use?**

- **Start with `IGEO7`** unless you have a specific reason not to — it's the
  most tested grid (the original one this library was built to reproduce
  bit-for-bit) and its aperture-7 hexagons are the shape most people picture
  when they think "hex grid".
- **Pick a projection variant (IVEA/RTEA) only if you need it** — they exist
  because DGGAL supports all three `VGCRadialVertex` orientations; if nothing
  in your project already commits you to IVEA or RTEA, `ISEA`-based grids
  (`IGEO7`/`ISEA3H`) are the default choice.
- **Pick aperture-3 (`*3H`) if you need one of these things `*7H` doesn't
  have**: sub-zones/tiling (see below — I3H-only), or you specifically want
  the rhombic addressing scheme (text ids like `"D8-1B8-A"` rather than
  Z7's digit strings like `"05460005"`).
- Every grid exposes the **exact same `Zone`/`Grid` API** shown throughout
  this tutorial — switching grids never changes your code shape, only which
  grid object you called it on:

```python
from py4dggs import get_grid

def describe(grid_name, lat, lon, res):
    grid = get_grid(grid_name)
    z = grid.zone_from_geo(lat, lon, res)
    print(f"{grid_name} @ res {res}: {z.text_id}  centroid={z.centroid}  "
          f"pentagon={z.is_pentagon}  neighbors={len(z.neighbors)}  children={len(z.children)}")

for name in ("IGEO7", "IVEA7H", "RTEA7H", "ISEA3H", "IVEA3H", "RTEA3H"):
    describe(name, 35.6762, 139.6503, 5)
```
```
IGEO7 @ res 5: 0546000  centroid=GeoPoint(lat=35.925883664238306, lon=139.78487902277107)  pentagon=False  neighbors=6  children=7
IVEA7H @ res 5: 0546005  centroid=GeoPoint(lat=35.408436551508785, lon=139.62086704524344)  pentagon=False  neighbors=6  children=7
RTEA7H @ res 5: 0546000  centroid=GeoPoint(lat=35.93541373871177, lon=139.61625379293372)  pentagon=False  neighbors=6  children=7
ISEA3H @ res 5: C8-2F-C  centroid=GeoPoint(lat=35.83525263481736, lon=141.2182032548381)  pentagon=False  neighbors=6  children=7
IVEA3H @ res 5: C8-2F-C  centroid=GeoPoint(lat=35.85201309402165, lon=140.97250301545193)  pentagon=False  neighbors=6  children=7
RTEA3H @ res 5: C8-2F-C  centroid=GeoPoint(lat=35.844051255213806, lon=141.09014750256864)  pentagon=False  neighbors=6  children=7
```

`get_grid("IGEO7")` is the same object as `IGEO7` — use `get_grid(name)` when
the grid is chosen at runtime (e.g. from a config value or CLI flag), and the
plain import (`from py4dggs import IGEO7`) when it's a compile-time choice.

## Neighbours (the k-ring)

Every zone knows its immediate neighbours — 6 for an ordinary hexagon, 5 for
one of the 12 pentagons:

```python
from py4dggs import IGEO7

tokyo = IGEO7.zone_from_geo(lat=35.6762, lon=139.6503, res=6)
print(sorted(n.text_id for n in tokyo.neighbors))
```
```
['05460000', '05460001', '05460004', '05460016', '05460052', '05460053']
```

```python
pentagon = IGEO7.zone_from_text("060")   # a known pentagon
print(pentagon.is_pentagon, len(pentagon.neighbors))
```
```
True 5
```

## Hierarchy: parent, children, and ancestors

Every zone at resolution `r > 0` has one or more **parents** at resolution
`r - 1`, and every zone has **children** at resolution `r + 1`.

```python
from py4dggs import IGEO7

tokyo = IGEO7.zone_from_geo(lat=35.6762, lon=139.6503, res=6)
print(tokyo.parent.text_id)                       # one level coarser
print(sorted(c.text_id for c in tokyo.children))  # one level finer

# walk up to a coarser ancestor by following .parent repeatedly
z = tokyo
for _ in range(3):
    z = z.parent
print(z.text_id, z.resolution)
```
```
0546000
['054600050', '054600051', '054600052', '054600053', '054600054', '054600055', '054600056']
05460 3
```

**Congruent vs. non-congruent hierarchy — a real gotcha to know about.** The
two topology families answer "how many parents does a cell have?" differently:

- **Z7 grids** (`IGEO7`/`IVEA7H`/`RTEA7H`) are *congruent*: every cell has
  **exactly one** parent (`.parent` is never ambiguous) and 7 children (6 for
  a pentagon) — a parent/child relationship is just appending/dropping one
  digit from the text id.
- **I3H grids** (`ISEA3H`/`IVEA3H`/`RTEA3H`) are *not* congruent: a cell can
  have **1 or 3** parents. `.parent` always gives you the *primary* one
  (`parents[0]`), but `.parents` gives you all of them, and
  `.is_centroid_child`/`.centroid_parent` tell you which of the (possibly 3)
  parents is the "centroid" one:

```python
from py4dggs import ISEA3H

z = ISEA3H.zone_from_text("B4-6-A")
print(z.parents)                                    # 3 parents, not 1
print([p.text_id for p in z.parents])
print(z.is_centroid_child)
print(z.centroid_parent.text_id)
print(sorted(c.text_id for c in z.children))
```
```
(Zone('A4-0-D'), Zone('A3-0-C'), Zone('A5-0-B'))
['A4-0-D', 'A3-0-C', 'A5-0-B']
False
A5-0-B
['B3-5-C', 'B3-5-D', 'B3-8-C', 'B4-3-D', 'B4-6-B', 'B4-6-C', 'B4-6-D']
```

If you're writing grid-agnostic code (working with whichever of the 6 grids
the caller passed in), use `.parents`/`.is_centroid_child` rather than
assuming `.parent` is the *only* parent — it's always safe to call, but on an
I3H grid it's only ever showing you one of up to three true parents.

## Sub-zones — generating many cells at once (I3H grids only)

This is the closest thing in this library to "generate a grid": instead of
computing one zone at a time, **sub-zones** give you *every* descendant of a
zone `relative_depth` levels down, as one ordered, indexable batch — the OGC
"descendants-at-depth" operation. It is currently implemented on the I3H grids
(`ISEA3H`/`IVEA3H`/`RTEA3H`); on the Z7 grids (`IGEO7`/`IVEA7H`/`RTEA7H`) these
five methods raise `NotImplementedError` for now.

Z7 sub-zones are **planned, not a non-goal** — a roadmap item, not a scope
boundary. `Grid` already dispatches sub-zones through its optional-capability
`getattr` pattern, so adding them means writing the three methods on `hex_a7`
following the `hex_a3` precedent; no architectural change is needed.

```python
from py4dggs import ISEA3H

tile = ISEA3H.zone_from_geo(lat=35.6762, lon=139.6503, res=6)
depth = 3

n = tile.count_sub_zones(depth)          # how many descendants, without generating them
print(n)                                  # 37

subs = tile.sub_zones(depth)              # the actual ordered tuple[Zone, ...]
print(len(subs), subs[0].text_id)

# the batch is a stable, indexable array: slot <-> cell in both directions
print(tile.first_sub_zone(depth) == subs[0])
print(tile.sub_zone_index(subs[5]) == 5)
print(tile.sub_zone_at_index(depth, 5) == subs[5])
```
```
37
37 E8-EF5-B
True
True
True
```

**Why this is the "DGGS-as-storage" pattern:** pick a coarse zone as a fixed
"tile" and a `relative_depth`, and `count_sub_zones` tells you the tile's
exact array length *before* you generate anything — exactly the contract a
raster/array storage format needs (Cloud-Optimized GeoTIFF, Zarr, ...), except
the "pixels" are equal-area DGGS cells instead of a row/column grid. A coarse
zone id plus a slot index stands in for a fine cell id without ever computing
or storing that fine id directly.

`sub_zone_index`/`sub_zone_at_index` currently build the full `sub_zones()`
list internally (`O(count_sub_zones)`, not `O(1)`) — fine for the tile sizes
this pattern is meant for, but worth knowing if you're reaching for a single
slot out of a very deep tile in a hot loop.

## Handling bad input

Every entry point raises `py4dggs.InvalidZoneError` (a `ValueError` subclass) on
malformed input, rather than silently returning nonsense:

```python
from py4dggs import IGEO7, InvalidZoneError

try:
    IGEO7.zone_from_geo(38.7, -9.1, res=999)
except InvalidZoneError as e:
    print(e)

try:
    IGEO7.zone_from_text("not-a-zone")
except InvalidZoneError as e:
    print(e)
```
```
resolution 999 out of range 0..19
bad Z7 text id 'not-a-zone'
```

## Runnable versions of everything here

Every example above also ships as a standalone script under `examples/`, so you
can run them instead of copying them:

```bash
python examples/01_first_zone.py     # this tutorial's first sections
python examples/05_tile_store.py     # the sub-zones/tile-store walkthrough
```

See `examples/README.md` for the full list. Both those scripts and every code
block in this file are executed by `tests/test_docs_examples.py`, so if something
here stops working, the test suite says so.

## Where to go next

- **`README.md`** — the "Verification against pydggal" section explains
  exactly what's proven correct (and how) for each of the 6 grids, and the
  "Grids"/"Sub-zones" sections are the terse reference version of what this
  tutorial walked through with examples.
- **`ARCHITECTURE.md`** — start here once you want to *read* the source: how
  a `Grid` is composed from a `Projection`+`Topology`+`Indexing`, and how to
  find the eC/DGGAL source line a given Python function ports.
