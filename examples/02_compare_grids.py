"""The same point on all six registered grids.

py4dggs registers six grids: three aperture-7 (Z7 digit addressing) and three
aperture-3 (I3H rhombic addressing), each in an ISEA / IVEA / RTEA projection
variant. They are genuinely different tessellations -- the ids are not
interchangeable between them, even when they look alike.

Run:  python examples/02_compare_grids.py
"""
from py4dggs import get_grid

LAT, LON, RES = 35.6762, 139.6503, 5  # Tokyo

print(f"{'grid':8}  {'text id':14}  {'centroid':40}  nbrs  kids")
print("-" * 78)
GRIDS = ["IGEO7", "IVEA7H", "RTEA7H", "ISEA3H", "IVEA3H", "RTEA3H"]
for name in GRIDS:
    g = get_grid(name)
    z = g.zone_from_geo(lat=LAT, lon=LON, res=RES)
    c = z.centroid
    print(f"{name:8}  {z.text_id:14}  "
          f"lat={c.lat:11.6f} lon={c.lon:12.6f}          "
          f"{len(z.neighbors):4}  {len(z.children):4}")

# The aperture-3 grids share an identical packing, so the SAME text id exists on
# all three -- but it denotes a different cell on each. Grid identity matters.
a = get_grid("ISEA3H").zone_from_geo(lat=LAT, lon=LON, res=RES)
b = get_grid("IVEA3H").zone_from_text(a.text_id)
print(f"\nsame text id on ISEA3H and IVEA3H: {a.text_id} == {b.text_id}")
print(f"  ...but the same zone?  {a == b}")
print(f"  their centroids differ by "
      f"{abs(a.centroid.lat - b.centroid.lat):.6f} deg latitude")
