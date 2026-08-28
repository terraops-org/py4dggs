"""Frozen golden-table conformance for py4dggs.ISEA3H vs the pydggal baseline
(point-keyed geometry; regenerate with `python generate.py ISEA3H`)."""
import json, os, pathlib, pytest
from py4dggs import ISEA3H

_DEFAULT = pathlib.Path(__file__).resolve().parent / "golden" / "isea3h"  # vendored, see tests/golden/PROVENANCE.md
_DIR = pathlib.Path(os.environ.get("DGGS_GOLDEN_TABLES_ISEA3H", _DEFAULT))
_HAVE = _DIR.exists()
pytestmark = pytest.mark.skipif(not _HAVE, reason=f"ISEA3H golden tables not at {_DIR}; run generate.py ISEA3H")
TOL = 0.0001

def _load(name):
    return json.loads((_DIR / f"{name}.json").read_text())["cases"] if _HAVE else []

def _srt(pairs):
    return sorted((round(a, 4), round(b, 4)) for a, b in pairs)

@pytest.mark.parametrize("c", _load("centroid"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_centroid(c):
    z = ISEA3H.zone_from_geo(c["lat"], c["lon"], c["res"])
    assert abs(z.centroid.lat - c["centroidLat"]) < TOL
    assert abs(z.centroid.lon - c["centroidLon"]) < TOL

@pytest.mark.parametrize("c", _load("vertices"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_vertices(c):
    z = ISEA3H.zone_from_geo(c["lat"], c["lon"], c["res"])
    dv = [(p.lat, p.lon) for p in z.vertices]
    assert len(dv) == len(c["vertices"])
    assert _srt(dv) == _srt([(a, b) for a, b in c["vertices"]])

@pytest.mark.parametrize("c", _load("textid"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_textid(c):
    val = int(c["value"])
    assert ISEA3H.indexing.to_text(val) == c["textId"]
    assert ISEA3H.indexing.from_text(c["textId"]) == val

@pytest.mark.parametrize("c", _load("hierarchy"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_hierarchy(c):
    val = int(c["value"])
    assert sorted(str(p) for p in ISEA3H.parents(val)) == c["parents"]
    assert sorted(str(cc) for cc in ISEA3H.children(val)) == c["children"]
    cp = ISEA3H.centroid_parent(val)
    assert (str(cp) if cp is not None else None) == c["centroidParent"]
    assert ISEA3H.is_centroid_child(val) == c["isCentroidChild"]

@pytest.mark.parametrize("c", _load("subzones"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_subzones(c):
    # value-keyed (A3) -> exact for every cell, no boundary xfail
    val = int(c["value"])
    for depth_str, expected in c["subZonesByDepth"].items():
        assert ISEA3H.sub_zones(val, int(depth_str)) == tuple(int(s) for s in expected)
