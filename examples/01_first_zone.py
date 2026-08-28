"""A single cell: from a coordinate to a zone and back.

Run:  python examples/01_first_zone.py
"""
from py4dggs import IGEO7, Zone

# Lisbon, at resolution 7. Every grid's entry point is the same: zone_from_geo.
zone = IGEO7.zone_from_geo(lat=38.7223, lon=-9.1393, res=7)

print("text id      ", zone.text_id)      # the canonical, shareable identifier
print("value        ", zone.value)        # the packed int -- equals DGGAL's own DGGRSZone
print("resolution   ", zone.resolution)
print("centroid     ", zone.centroid)     # NOT the point you asked for: the CELL's centre
print("is pentagon  ", zone.is_pentagon)  # 12 cells per resolution are pentagons
print("vertices     ", len(zone.vertices), "corners")

# Three representations of the same cell, all interchangeable.
assert IGEO7.zone_from_text(zone.text_id) == zone
assert Zone(IGEO7, zone.value) == zone

# Quantization is many-to-one: every point inside the cell gives the same zone.
nearby = IGEO7.zone_from_geo(lat=38.7224, lon=-9.1394, res=7)
print("\nnearby point -> same cell:", nearby == zone)

# Resolution controls cell size. Coarser = fewer, bigger cells.
for res in (3, 7, 11):
    z = IGEO7.zone_from_geo(lat=38.7223, lon=-9.1393, res=res)
    print(f"  res {res:2d}: {z.text_id}")
