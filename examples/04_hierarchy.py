"""Hierarchy: parents and children -- and why the two grid families differ.

The Z7 grids use a CONGRUENT digit hierarchy: a child's id is its parent's id
plus one digit, so every cell has exactly 1 parent and 7 children.

The I3H grids use DGGAL's GEOMETRIC hierarchy, which is non-congruent: a cell can
have 3 parents, and children are found by quantizing the cell's own boundary
vertices at the next resolution. Both are correct; they are different relations.

Run:  python examples/04_hierarchy.py
"""
from py4dggs import IGEO7, ISEA3H

# --- Z7: congruent, id-based ------------------------------------------------ #
z = IGEO7.zone_from_geo(lat=-33.8688, lon=151.2093, res=6)  # Sydney
parent = z.parents[0]
print("Z7 (IGEO7) -- congruent digit hierarchy")
print(f"  zone      {z.text_id}")
print(f"  parent    {parent.text_id}   (drop the last digit)")
print(f"  parents   {len(z.parents)}")
print(f"  children  {len(parent.children)}  -> {[c.text_id for c in parent.children]}")
assert z in parent.children
assert z.text_id.startswith(parent.text_id), "a Z7 child id extends its parent's"

# --- I3H: geometric, non-congruent ------------------------------------------ #
a = ISEA3H.zone_from_geo(lat=-33.8688, lon=151.2093, res=6)
print("\nI3H (ISEA3H) -- geometric hierarchy")
print(f"  zone      {a.text_id}")
print(f"  parents   {len(a.parents)}  -> {[p.text_id for p in a.parents]}")
print(f"  children  {len(a.children)}")
print(f"  centroid parent: {a.centroid_parent.text_id if a.centroid_parent else None}")
print(f"  is centroid child: {a.is_centroid_child}")

# Walking up: ancestors all the way to resolution 0.
print("\nwalking IGEO7 up to the base cell:")
cur = z
while cur.parents:
    cur = cur.parents[0]
    print(f"   res {cur.resolution}: {cur.text_id}")
