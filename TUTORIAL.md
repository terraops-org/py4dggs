# Tutorial - using `py4dggs` to work with real cells

This is a task-oriented walkthrough, distinct from the other two docs:

- **`README.md`**  what's verified and the terse API reference.
- **`ARCHITECTURE.md`**  how the code is structured, for studying/extending it.
- **This file**  how to *do things* with it, with worked, runnable examples.

Every snippet below was actually run against this exact library version, the
printed output is real, not illustrative. Run them yourself with `uv run python`
from the repo root, pasting the blocks in order: this is a walkthrough, so some
blocks reuse names that earlier ones defined.

`uv run` is what puts the package on the import path. A bare `python` will not
find it, because `uv sync` installs into `.venv/` rather than your system
interpreter; for a plain REPL, activate that environment first with
`source .venv/bin/activate`.

## Install

```bash
uv add py4dggs           # or: pip install py4dggs
```

Zero runtime dependencies - nothing else gets pulled in.

That is the way in for *using* the library. It does not work from inside a clone
of this repo, though: `uv add py4dggs` there fails with a self-dependency error,
because the project you are standing in is itself named `py4dggs`. To follow this
tutorial against a clone, or to run the test suite, install it as a project:

```bash
git clone https://github.com/terraops-org/py4dggs.git
cd py4dggs && uv sync    # or: pip install -e .
```

Both the PyPI distribution name and the import name are `py4dggs`, read as "Python for
DGGS". The other package you will see alongside it is [`dggal`](https://pypi.org/project/dggal/)
(pydggal), DGGAL's official Python binding to the `libdggal` C library; `py4dggs` is the
pure-Python counterpart to it, no C library required.

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

`tokyo.centroid` is the cell's *own* center (not your input lat/lon) - a
`GeoPoint(lat, lon)` in WGS84 degrees. `tokyo.vertices` gives you the cell's
boundary as a list of `GeoPoint`s (6 for a hexagon, 5 for one of the 12
pentagons every DGGS grid has).

## Converting between representations

A zone has three interchangeable forms, a `Zone` object, its canonical
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

`Zone` is a plain immutable `(grid, value)` pair, cheap to construct, safe to
put in a `set`/use as a `dict` key (it's hashable), and comparing two `Zone`s
only makes sense if they came from the *same* grid.

## Choosing a grid

`py4dggs` registers **six** grids, every combination of 3 projections ×
2 topology/indexing families:

| | **Z7** (aperture-7, digit-addressed) | **I3H** (aperture-3, rhombic-addressed) |
|---|---|---|
| **ISEA** projection | `IGEO7` | `ISEA3H` |
| **IVEA** projection | `IVEA7H` | `IVEA3H` |
| **RTEA** projection | `RTEA7H` | `RTEA3H` |

**Which one should you use?**

- **Start with `IGEO7`** unless you have a specific reason not to  it's the
  most tested grid (the original one this library was built to reproduce
  bit-for-bit) and its aperture-7 hexagons are the shape most people picture
  when they think "hex grid".
- **Pick a projection variant (IVEA/RTEA) only if you need it**, they exist
  because DGGAL supports all three `VGCRadialVertex` orientations; if nothing
  in your project already commits you to IVEA or RTEA, `ISEA`-based grids
  (`IGEO7`/`ISEA3H`) are the default choice.
- **Pick aperture-3 (`*3H`) if you need one of these things `*7H` doesn't
  have**: sub-zones/tiling (see below - I3H-only), or you specifically want
  the rhombic addressing scheme (text ids like `"D8-1B8-A"` rather than
  Z7's digit strings like `"05460005"`).
- Every grid exposes the **exact same `Zone`/`Grid` API** shown throughout
  this tutorial - switching grids never changes your code shape, only which
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

`get_grid("IGEO7")` is the same object as `IGEO7`, use `get_grid(name)` when
the grid is chosen at runtime (e.g. from a config value or CLI flag), and the
plain import (`from py4dggs import IGEO7`) when it's a compile-time choice.

## Neighbours (the k-ring)

Every zone knows its immediate neighbours, 6 for an ordinary hexagon, 5 for
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

**Congruent vs. non-congruent hierarchy a real gotcha to know about.** The
two topology families answer "how many parents does a cell have?" differently:

- **Z7 grids** (`IGEO7`/`IVEA7H`/`RTEA7H`) are *congruent*: every cell has
  **exactly one** parent (`.parent` is never ambiguous) and 7 children (6 for
  a pentagon) - a parent/child relationship is just appending/dropping one
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
assuming `.parent` is the *only* parent - it's always safe to call, but on an
I3H grid it's only ever showing you one of up to three true parents.

## Sub-zones - generating many cells at once (I3H grids only)

This is the closest thing in this library to "generate a grid": instead of
computing one zone at a time, **sub-zones** give you *every* descendant of a
zone `relative_depth` levels down, as one ordered, indexable batch - the OGC
"descendants-at-depth" operation. It is currently implemented on the I3H grids
(`ISEA3H`/`IVEA3H`/`RTEA3H`); on the Z7 grids (`IGEO7`/`IVEA7H`/`RTEA7H`) these
five methods raise `NotImplementedError` for now.

Z7 sub-zones are **planned, not a non-goal**, a roadmap item, not a scope
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

`sub_zone_index`/`sub_zone_at_index` currently build the full `sub_zones()`
list internally (`O(count_sub_zones)`, not `O(1)`) - fine for the tile sizes
this pattern is meant for, but worth knowing if you're reaching for a single
slot out of a very deep tile in a hot loop.

## Exporting cells as GeoJSON

`zone.vertices` gives you a cell's boundary exactly as DGGAL reports it. That is
the library's correctness contract, so it stays raw: longitudes in [-180, 180],
in the engine's own order, with no repeated closing point. The `py4dggs.geojson`
module sits on top of that and emits RFC 7946 geometry you can hand straight to
a map. It is exported from the top-level package, so `from py4dggs import
geojson` finds it; the examples below import the three functions directly:

```python
from py4dggs import IGEO7
from py4dggs.geojson import zone_geometry

lisbon = IGEO7.zone_from_geo(lat=38.7223, lon=-9.1393, res=8)
geom = zone_geometry(lisbon)
ring = geom["coordinates"][0]

print(geom["type"], len(ring))
print(ring[0] == ring[-1])
```
```
Polygon 7
True
```

A hexagon has 6 vertices but the ring has 7, because RFC 7946 wants the first
point repeated at the end. Every emitted ring is also wound counterclockwise, so
consumers that care about orientation (PostGIS, Shapely) get what they expect.

`zone_feature` wraps a single cell and `feature_collection` wraps any iterable of
them, which is usually what you want. By default each feature carries its zone's
text id; pass `properties` to replace that, either as one dict applied to every
feature or as a callable taking a zone:

```python
from py4dggs.geojson import zone_feature, feature_collection

print(zone_feature(lisbon)["properties"])

patch = feature_collection([lisbon, *lisbon.neighbors])
print(patch["type"], len(patch["features"]))

labelled = feature_collection(
    [lisbon, *lisbon.neighbors],
    lambda z: {"id": z.text_id, "res": z.resolution},
)
print(labelled["features"][0]["properties"])
```
```
{'zone': '0064156546'}
FeatureCollection 7
{'id': '0064156546', 'res': 8}
```

The result is plain dicts, so writing it out is just
`json.dump(patch, open("patch.geojson", "w"))` and the file drops straight into
QGIS or geojson.io. Any iterable works as the input, so
`feature_collection(tile.sub_zones(2))` exports a whole tile from the sub-zones
section above.

## GeoJSON's two hard cases: the antimeridian and the poles

Two kinds of cell need more than copying coordinates across, and they are the
reason to call this module instead of writing the loop yourself.

A cell sitting on +/-180 has vertices on both sides of it. Each longitude is
individually in range, so nothing looks wrong until a planar reader joins them
up and draws a sliver the long way round the globe:

```python
from py4dggs import IGEO7
from py4dggs.geojson import zone_geometry

straddler = IGEO7.zone_from_geo(lat=-20.0, lon=179.99, res=5)
lons = [v.lon for v in straddler.vertices]
print(round(max(lons) - min(lons), 2))

geom = zone_geometry(straddler)
print(geom["type"], len(geom["coordinates"]))
```
```
359.89
MultiPolygon 2
```

That cell is barely a degree across, but its raw longitude span is 359.89
degrees. `zone_geometry` cuts it at the antimeridian and returns the two halves
as a `MultiPolygon`, per RFC 7946 section 3.1.9.

The other case is a cell containing a pole. It winds a full 360 degrees in
longitude yet carries no vertex at latitude +/-90, so it cannot be closed as a
planar ring at all. `zone_geometry` re-anchors that boundary to run from -180 to
+180 and seals it across the pole itself, giving back a `Polygon` with the two
pole corners added. There is no example here because picking the cell that
encloses a pole is fiddly (quantizing the exact pole is a boundary tie-break, so
the cell you get back is not always the one that contains it); see
`pole_cell()` in `tests/test_geojson.py` for a reliable way to find it.

Finally, a zone with no geometry (DGGAL's nullZone) raises `InvalidZoneError`
instead of returning a degenerate shape, which brings us to how the rest of the
library reports bad input.

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
uv run python examples/01_first_zone.py   # this tutorial's first sections
uv run python examples/05_tile_store.py   # the sub-zones/tile-store walkthrough
```

See `examples/README.md` for the full list. Both those scripts and every code
block in this file are executed by `tests/test_docs_examples.py`, so if something
here stops working, the test suite says so.

## Where to go next

- **`README.md`** the "Verification against pydggal" section explains
  exactly what's proven correct (and how) for each of the 6 grids, and the
  "Grids"/"Sub-zones" sections are the terse reference version of what this
  tutorial walked through with examples.
- **`ARCHITECTURE.md`** start here once you want to *read* the source: how
  a `Grid` is composed from a `Projection`+`Topology`+`Indexing`, and how to
  find the eC/DGGAL source line a given Python function ports.
