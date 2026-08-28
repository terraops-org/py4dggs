"""Sub-zones: using a DGGS as a fixed-size tile store.

Pick a coarse zone as a "tile" and a relative_depth, and sub-zones give you a
stable, ordered mapping between an array slot and a fine-resolution cell -- the
raster-band model, but DGGS-native (equal-area cells, no reprojection).

Sub-zones are currently implemented on the I3H grids only; the Z7 grids raise
NotImplementedError. That is a roadmap item, not a permanent limitation.

Run:  python examples/05_tile_store.py
"""
from py4dggs import ISEA3H, IGEO7, InvalidZoneError

tile = ISEA3H.zone_from_geo(lat=35.6762, lon=139.6503, res=6)
DEPTH = 3

# ALWAYS ask the count first: it is a cheap closed form, while materialising the
# sequence is not. This is the gate that keeps a caller-supplied depth safe.
n = tile.count_sub_zones(DEPTH)
print(f"tile {tile.text_id} at depth {DEPTH}: {n} sub-zones")

cells = tile.sub_zones(DEPTH)
print(f"materialised {len(cells)}; first={cells[0].text_id} last={cells[-1].text_id}")

# The mapping is stable in both directions: slot -> cell, and cell -> slot.
slot = 17
cell = tile.sub_zone_at_index(DEPTH, slot)
print(f"\nslot {slot} -> {cell.text_id} -> slot {tile.sub_zone_index(cell)}")
assert tile.sub_zone_index(cell) == slot

# Index 0 is the centroid descendant, available without building the list.
print(f"first_sub_zone (no materialisation): {tile.first_sub_zone(DEPTH).text_id}")

# A zone is its own sub-zone at depth 0.
print(f"depth 0 -> {tile.sub_zones(0)[0].text_id}, index {tile.sub_zone_index(tile)}")

# Depths that would not fit in memory are refused rather than attempted.
coarse = ISEA3H.zone_from_geo(lat=0.0, lon=0.0, res=0)
huge = coarse.count_sub_zones(33)
print(f"\ncount at depth 33 is answerable: {huge:,}")
try:
    coarse.sub_zones(33)
except InvalidZoneError as e:
    print(f"...but building it is refused: {e}")

# Z7 grids do not have a sub-zone order yet.
try:
    IGEO7.zone_from_geo(lat=0.0, lon=0.0, res=4).sub_zones(2)
except NotImplementedError as e:
    print(f"\nIGEO7 sub-zones: NotImplementedError({e}) -- planned, see TUTORIAL.md")
