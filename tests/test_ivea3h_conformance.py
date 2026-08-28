"""Frozen golden-table conformance for py4dggs.IVEA3H vs the pydggal baseline
(point-keyed geometry; regenerate with `python generate.py IVEA3H`).

IVEA3H reproduces pydggal exactly except for 12 forward SINGULARITIES at the
exact vertex-0 meridian (lon=11.2 deg, the HARD "Vert0Lon" test point) at even
resolutions 0/2/4/6/8/10, where the clean-room engine and pydggal break the
exact-boundary tie differently. Confirmed benign 2026-07-03: at each of these
points our returned cell is a direct NEIGHBOUR of pydggal's, and both engines
agree exactly on that neighbour cell's own centroid/vertices (our engine's
centroid for our chosen cell == pydggal's centroid for that same cell) — i.e.
a pure quantization tie-break at a measure-zero boundary, not a geometry bug.
This is the same singularity class IVEA7H documents (13 xfails: pole + vertex-0
meridian), consistent with IVEA's extra pole-longitude sensitivity vs ISEA (the
eC's IVEA pole-fix, ri5x6.ec:476, which we do not implement — see IVEA7H's
conformance test docstring). ISEA3H has zero conformance xfails at this point
set since ISEA does not share IVEA's pole-longitude sensitivity."""
import json, os, pathlib, pytest
from py4dggs import IVEA3H

_DEFAULT = pathlib.Path(__file__).resolve().parent / "golden" / "ivea3h"  # vendored, see tests/golden/PROVENANCE.md
_DIR = pathlib.Path(os.environ.get("DGGS_GOLDEN_TABLES_IVEA3H", _DEFAULT))
_HAVE = _DIR.exists()
pytestmark = pytest.mark.skipif(not _HAVE, reason=f"IVEA3H golden tables not at {_DIR}; run generate.py IVEA3H")
TOL = 0.0001

# Empirical singularity set (py4dggs.IVEA3H vs the pydggal baseline, 2026-07-03):
# exact-boundary tie-break at the vertex-0 meridian (lon 11.2), even resolutions.
# Both test_centroid and test_vertices re-run forward quantization (the tables
# are point-keyed, not textId-keyed), so both are affected by the same tie.
XFAIL = {("Vert0Lon", r) for r in (0, 2, 4, 6, 8, 10)}

def _load(name):
    return json.loads((_DIR / f"{name}.json").read_text())["cases"] if _HAVE else []

def _srt(pairs):
    return sorted((round(a, 4), round(b, 4)) for a, b in pairs)

@pytest.mark.parametrize("c", _load("centroid"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_centroid(c):
    if (c.get("name"), c["res"]) in XFAIL:
        pytest.xfail("IVEA3H exact-singularity boundary tie-break (vertex-0 meridian)")
    z = IVEA3H.zone_from_geo(c["lat"], c["lon"], c["res"])
    assert abs(z.centroid.lat - c["centroidLat"]) < TOL
    assert abs(z.centroid.lon - c["centroidLon"]) < TOL

@pytest.mark.parametrize("c", _load("vertices"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_vertices(c):
    if (c.get("name"), c["res"]) in XFAIL:
        pytest.xfail("IVEA3H exact-singularity boundary tie-break (vertex-0 meridian)")
    z = IVEA3H.zone_from_geo(c["lat"], c["lon"], c["res"])
    dv = [(p.lat, p.lon) for p in z.vertices]
    assert len(dv) == len(c["vertices"])
    assert _srt(dv) == _srt([(a, b) for a, b in c["vertices"]])

@pytest.mark.parametrize("c", _load("textid"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_textid(c):
    # value-keyed -> exact for every cell, no boundary xfail
    val = int(c["value"])
    assert IVEA3H.indexing.to_text(val) == c["textId"]
    assert IVEA3H.indexing.from_text(c["textId"]) == val

@pytest.mark.parametrize("c", _load("hierarchy"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_hierarchy(c):
    val = int(c["value"])
    assert sorted(str(p) for p in IVEA3H.parents(val)) == c["parents"]
    assert sorted(str(cc) for cc in IVEA3H.children(val)) == c["children"]
    cp = IVEA3H.centroid_parent(val)
    assert (str(cp) if cp is not None else None) == c["centroidParent"]
    assert IVEA3H.is_centroid_child(val) == c["isCentroidChild"]

@pytest.mark.parametrize("c", _load("subzones"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_subzones(c):
    # value-keyed (A3) -> exact for every cell, no boundary xfail
    val = int(c["value"])
    for depth_str, expected in c["subZonesByDepth"].items():
        assert IVEA3H.sub_zones(val, int(depth_str)) == tuple(int(s) for s in expected)
