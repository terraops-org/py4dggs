"""Neighbours (the k-ring), including the pentagon case.

Run:  python examples/03_neighbours.py
"""
from py4dggs import IGEO7, ISEA3H

zone = IGEO7.zone_from_geo(lat=52.3676, lon=4.9041, res=6)  # Amsterdam
neighbours = zone.neighbors

print(f"{zone.text_id} has {len(neighbours)} neighbours:")
for n in neighbours:
    print(f"   {n.text_id}  {n.centroid}")

# Adjacency is symmetric: each neighbour lists us back.
assert all(zone in n.neighbors for n in neighbours), "adjacency must be mutual"
print("\nsymmetric:", True)

# Every icosahedral DGGS has exactly 12 pentagon cells per resolution, and a
# pentagon has 5 neighbours, not 6. Code that assumes 6 will be wrong 12 times.
pentagon = IGEO7.zone_from_text("060")
print(f"\npentagon {pentagon.text_id}: is_pentagon={pentagon.is_pentagon}, "
      f"{len(pentagon.neighbors)} neighbours")

# The aperture-3 grids use an exact topological neighbour rule, not a geometric
# approximation -- they agree with DGGAL cell-for-cell.
a3 = ISEA3H.zone_from_geo(lat=52.3676, lon=4.9041, res=6)
print(f"\n{a3.text_id} (ISEA3H): {len(a3.neighbors)} neighbours")
