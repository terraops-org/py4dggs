#!/usr/bin/env python3
"""Generate IGEO7/Z7 golden tables from pydggal.

pydggal return shapes (probed 2026-06-16, pydggal v0.0.6 / ecrt 0.0.6):
  zone                                  -> int (DGGRSZone); getZoneTextID accepts int or DGGRSZone
  getZoneFromWGS84Centroid(level, GeoPoint) -> zone
  getZoneTextID(zone)                   -> str like "0064156"
  getZoneFromTextID(textId)             -> zone
  getZoneWGS84Centroid(zone)            -> GeoPoint with .lat/.lon in DEGREES
  getZoneWGS84Vertices(zone)            -> ecrt Array of GeoPoint (.lat/.lon degrees)
  getZoneNeighbors(zone)                -> ecrt Array of DGGRSZone (6 for hex, 5 for pentagon)
  getZoneChildren/getZoneParents        -> NON-CONGRUENT geometric hierarchy (13 / 2) -- NOT used.

Hierarchy here is the Z7 INDEX hierarchy (what the explorer + JS port use):
  parent  = drop the last direction digit (1 parent)
  children = append digit 0..6, keep those pydggal validates as real cells
             (7 for hexagons, 6 for pentagons)

Run: python generate.py
"""
import sys
import json
import datetime
from pathlib import Path

from dggal import *

app = Application(appGlobals=globals())
pydggal_setup(app)

# Grid to generate for: default ISEA7H_Z7 (IGEO7, tables/); pass another DGGRS
# class name (e.g. IVEA7H_Z7) to generate its variant tables under tables/<name>/.
_GRID_NAME = sys.argv[1] if len(sys.argv) > 1 else "ISEA7H_Z7"
grid = globals()[_GRID_NAME]()

# VENDORED DIVERGENCE (the only edit to this file vs igeo7-spec's copy -- see
# PROVENANCE.md): upstream writes to `<script dir>/tables`, but here the tables
# ARE the script's own directory (tests/golden/), so there is no extra level.
TABLES = Path(__file__).parent
if _GRID_NAME != "ISEA7H_Z7":
    TABLES = TABLES / _GRID_NAME.replace("_Z7", "").lower()  # e.g. tests/golden/ivea7h
TABLES.mkdir(parents=True, exist_ok=True)
ROUND = 7

# --- test point set: cities PLUS hard cases (see RUNNER_CONTRACT / spec) ---
CITIES = [
    ("Lisbon", 38.7223, -9.1393), ("Tokyo", 35.6762, 139.6503),
    ("NewYork", 40.7128, -74.0060), ("SaoPaulo", -23.5505, -46.6333),
    ("Sydney", -33.8688, 151.2093), ("London", 51.5074, -0.1278),
    ("Nairobi", -1.2921, 36.8219), ("Reykjavik", 64.1466, -21.9426),
]
HARD = [
    ("NorthPole", 90.0, 0.0), ("SouthPole", -90.0, 0.0),
    ("Antimeridian", 0.0, 180.0), ("AntimeridianNeg", 0.0, -179.999),
    ("Equator0", 0.0, 0.0), ("Vert0Lon", 0.0, 11.2),
]
POINTS = [{"name": n, "lat": la, "lon": lo} for (n, la, lo) in CITIES + HARD]


def _provenance():
    return {
        "source": "pydggal",
        "grid": _GRID_NAME,
        "dggal_version": "0.0.6",
        "binding_commit": "v0.0.6",
        "generator": "verify/generate.py",
        "generated": datetime.datetime.now(datetime.timezone.utc)
                          .replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "status": "DGGAL-ecosystem reference — pending lab confirmation",
    }


def _write(name, cases):
    path = TABLES / f"{name}.json"
    path.write_text(json.dumps({"provenance": _provenance(), "cases": cases}, indent=2))
    print(f"wrote {path.name}  ({len(cases)} cases)")


def _zone_at(pt, res):
    gpt = GeoPoint(); gpt.lat = pt["lat"]; gpt.lon = pt["lon"]
    return grid.getZoneFromWGS84Centroid(res, gpt)


def _textid(zone):
    return grid.getZoneTextID(zone)


def _valid_textid(textId):
    """True iff textId round-trips through pydggal (a real cell)."""
    try:
        z = grid.getZoneFromTextID(textId)
        return _textid(z) == textId
    except Exception:
        return False


def gen_forward():
    cases = []
    for pt in POINTS:
        for res in range(11):  # 0..10
            cases.append({
                "name": pt["name"],
                "lat": round(pt["lat"], ROUND),
                "lon": round(pt["lon"], ROUND),
                "res": res,
                "textId": _textid(_zone_at(pt, res)),
            })
    _write("forward", cases)


def gen_inverse():
    cases, seen = [], set()
    for pt in POINTS:
        for res in range(11):
            zone = _zone_at(pt, res)
            tid = _textid(zone)
            if tid in seen:
                continue
            seen.add(tid)
            c = grid.getZoneWGS84Centroid(zone)
            cases.append({
                "textId": tid, "res": res,
                "centroidLat": round(float(c.lat), ROUND),
                "centroidLon": round(float(c.lon), ROUND),
            })
    _write("inverse", cases)


def gen_hierarchy():
    cases, seen = [], set()
    for pt in POINTS:
        for res in range(1, 9):  # need a parent (res>=1)
            tid = _textid(_zone_at(pt, res))
            if tid in seen:
                continue
            seen.add(tid)
            parent = tid[:-1]  # Z7 index parent: drop last direction digit
            children = [tid + d for d in "0123456" if _valid_textid(tid + d)]
            cases.append({
                "textId": tid, "res": res,
                "parent": parent,
                "children": sorted(children),
            })
    _write("hierarchy", cases)


def gen_vertices():
    cases, seen = [], set()
    for pt in POINTS:
        for res in range(0, 8):
            zone = _zone_at(pt, res)
            tid = _textid(zone)
            if tid in seen:
                continue
            seen.add(tid)
            verts = [[round(float(v.lat), ROUND), round(float(v.lon), ROUND)]
                     for v in grid.getZoneWGS84Vertices(zone)]
            cases.append({"textId": tid, "res": res, "vertices": verts})
    _write("vertices", cases)


def gen_kring():
    cases, seen, skipped = [], set(), 0
    for pt in POINTS:
        for res in range(2, 9):
            zone = _zone_at(pt, res)
            tid = _textid(zone)
            if tid in seen:
                continue
            seen.add(tid)
            try:
                disk = sorted({tid, *(_textid(n) for n in grid.getZoneNeighbors(zone))})
            except OverflowError:
                # pydggal v0.0.6 marshals neighbour zones as int64; uint64 high-bit
                # zones overflow. Record visibly rather than silently dropping.
                skipped += 1
                cases.append({"textId": tid, "res": res, "k": 1,
                              "skipped": "pydggal-int64-overflow"})
                continue
            cases.append({"textId": tid, "res": res, "k": 1, "disk": disk})
    _write("kring", cases)
    if skipped:
        print(f"  WARNING: {skipped}/{len(cases)} kring cases skipped (pydggal int64 overflow)")


def gen_isea3h_centroid():
    """Point-keyed centroid table for ISEA3H (no text-id — deferred for our lib)."""
    cases = []
    for pt in POINTS:
        for res in range(11):  # 0..10
            zone = _zone_at(pt, res)
            c = grid.getZoneWGS84Centroid(zone)
            cases.append({
                "name": pt["name"],
                "lat": round(pt["lat"], ROUND),
                "lon": round(pt["lon"], ROUND),
                "res": res,
                "centroidLat": round(float(c.lat), ROUND),
                "centroidLon": round(float(c.lon), ROUND),
            })
    _write("centroid", cases)


def gen_isea3h_vertices():
    """Point-keyed vertices table for ISEA3H."""
    cases = []
    for pt in POINTS:
        for res in range(11):  # 0..10
            zone = _zone_at(pt, res)
            verts = [[round(float(v.lat), ROUND), round(float(v.lon), ROUND)]
                     for v in grid.getZoneWGS84Vertices(zone)]
            cases.append({
                "name": pt["name"],
                "lat": round(pt["lat"], ROUND),
                "lon": round(pt["lon"], ROUND),
                "res": res,
                "vertices": verts,
            })
    _write("vertices", cases)


def gen_isea3h_neighbours():
    """Point-keyed neighbour-centroid table for ISEA3H.

    I3H zones do NOT overflow pydggal's int64 marshalling the way some ISEA7H_Z7
    zones do (see gen_kring's OverflowError note), but we guard defensively
    anyway and record "skipped" rather than crash if a future pydggal version
    regresses this.
    """
    cases = []
    for pt in POINTS:
        for res in range(11):  # 0..10
            zone = _zone_at(pt, res)
            entry = {
                "name": pt["name"],
                "lat": round(pt["lat"], ROUND),
                "lon": round(pt["lon"], ROUND),
                "res": res,
            }
            try:
                neighbour_centroids = sorted(
                    [round(float(c.lat), ROUND), round(float(c.lon), ROUND)]
                    for c in (grid.getZoneWGS84Centroid(n) for n in grid.getZoneNeighbors(zone))
                )
                entry["neighbourCentroids"] = neighbour_centroids
            except OverflowError:
                entry["neighbourCentroids"] = "skipped"
            cases.append(entry)
    _write("neighbours", cases)


def gen_isea3h_textid():
    """Value-keyed text-id table for the I3H family: pydggal zone int -> canonical
    text-id. Value-keyed (not point-keyed like the geometry tables) so it tests
    text-id purely, decoupled from the quantization boundary tie-breaks already
    frozen in centroid/vertices."""
    cases = []
    for pt in POINTS:
        for res in range(11):  # 0..10
            zone = _zone_at(pt, res)
            cases.append({
                "name": pt["name"],
                "res": res,
                "value": str(int(zone)),
                "textId": _textid(zone),
            })
    _write("textid", cases)


def gen_isea3h_hierarchy():
    """Value-keyed hierarchy table for the I3H family (A2): pydggal zone int ->
    parents / children / centroidParent / isCentroidChild. Value-keyed (like
    text-id) so it verifies the hierarchy purely, immune to the quantization
    boundary tie-break frozen in the geometry tables."""
    cases = []
    for pt in POINTS:
        for res in range(11):  # 0..10
            zone = _zone_at(pt, res)
            cpv = int(grid.getZoneCentroidParent(zone))
            cases.append({
                "name": pt["name"],
                "res": res,
                "value": str(int(zone)),
                "parents": sorted(str(int(p)) for p in grid.getZoneParents(zone)),
                "children": sorted(str(int(c)) for c in grid.getZoneChildren(zone)),
                "centroidParent": None if cpv == (1 << 64) - 1 else str(cpv),
                "isCentroidChild": bool(grid.isZoneCentroidChild(zone)),
            })
    _write("hierarchy", cases)


def gen_isea3h_subzones():
    """Value-keyed sub-zone table for the I3H family (A3): pydggal zone int ->
    ordered sub-zone int lists at relative depths 1-3. Value-keyed (like
    text-id/hierarchy) so it verifies Grid.sub_zones purely, immune to the
    quantization boundary tie-break frozen in the geometry tables."""
    cases = []
    for pt in POINTS:
        for res in range(11):  # 0..10
            zone = _zone_at(pt, res)
            by_depth = {}
            for depth in (1, 2, 3):
                subs = grid.getSubZones(zone, depth)
                by_depth[str(depth)] = [str(int(subs[i])) for i in range(subs.count)]
            cases.append({
                "name": pt["name"],
                "res": res,
                "value": str(int(zone)),
                "subZonesByDepth": by_depth,
            })
    _write("subzones", cases)


if __name__ == "__main__":
    if _GRID_NAME in ("ISEA3H", "IVEA3H", "RTEA3H"):
        # Point-keyed geometry tables + value-keyed text-id (A1) + hierarchy (A2)
        # + sub-zones (A3).
        gen_isea3h_centroid()
        gen_isea3h_vertices()
        gen_isea3h_neighbours()
        gen_isea3h_textid()
        gen_isea3h_hierarchy()
        gen_isea3h_subzones()
    else:
        gen_forward()
        gen_inverse()
        gen_hierarchy()
        gen_vertices()
        gen_kring()
