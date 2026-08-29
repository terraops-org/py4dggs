"""Turning cells into RFC 7946 GeoJSON, including the antimeridian case.

Run:  python examples/07_geojson.py
"""
import json

from py4dggs import IGEO7, ISEA3H
from py4dggs.geojson import zone_geometry, zone_feature, feature_collection

lisbon = IGEO7.zone_from_geo(lat=38.7223, lon=-9.1393, res=8)
geom = zone_geometry(lisbon)
ring = geom["coordinates"][0]

print("one cell:")
print(f"  zone            {lisbon.text_id}")
print(f"  geometry        {geom['type']}")
print(f"  ring points     {len(ring)} (6 vertices plus the repeated closing point)")
print(f"  ring is closed  {ring[0] == ring[-1]}")

print("\none feature (properties default to the zone's text id):")
print(f"  {zone_feature(lisbon)['properties']}")

print("\na cell and its neighbours as a FeatureCollection:")
patch = feature_collection([lisbon, *lisbon.neighbors])
print(f"  features        {len(patch['features'])}")
print(f"  serialised      {len(json.dumps(patch))} bytes")

print("\nany iterable of zones works, so a whole I3H tile exports in one call:")
tile = ISEA3H.zone_from_geo(lat=35.6762, lon=139.6503, res=6)
tile_fc = feature_collection(
    tile.sub_zones(2),
    lambda z: {"id": z.text_id, "res": z.resolution},
)
print(f"  tile            {tile.text_id}")
print(f"  sub-zones       {len(tile_fc['features'])}")
print(f"  first feature   {tile_fc['features'][0]['properties']}")

print("\na cell on the antimeridian is cut into two parts:")
straddler = IGEO7.zone_from_geo(lat=-20.0, lon=179.99, res=5)
lons = [v.lon for v in straddler.vertices]
straddler_geom = zone_geometry(straddler)
print(f"  zone            {straddler.text_id}")
print(f"  raw lon span    {max(lons) - min(lons):.2f} degrees, for a cell ~1 degree across")
print(f"  geometry        {straddler_geom['type']} with {len(straddler_geom['coordinates'])} parts")

print("\nto view any of this, write it out and open it in QGIS or geojson.io:")
print('  json.dump(patch, open("patch.geojson", "w"))')
