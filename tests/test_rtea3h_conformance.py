"""Frozen golden-table conformance for py4dggs.RTEA3H vs the pydggal baseline
(point-keyed geometry; regenerate with `python generate.py RTEA3H`).

Mirrors test_ivea3h_conformance.py / test_isea3h_conformance.py. This table
checks forward-quantize + centroid, and forward-quantize + vertices -- it does
NOT check neighbours, so the k-ring geometric-approximation tolerance
documented in test_rtea3h_fuzz.py's module docstring does not apply here.

RTEA3H reproduces pydggal exactly except for 22 forward SINGULARITIES at the
exact vertex-0 meridian (lon=11.2 deg, the HARD "Vert0Lon" test point) at
EVERY resolution 0-10, where the clean-room engine and pydggal break the
exact-boundary tie differently. Confirmed benign 2026-07-03: at each of these
points our returned cell is a direct NEIGHBOUR of pydggal's, and both engines
agree exactly on that neighbour cell's own centroid/vertices (our engine's
centroid/vertices for pydggal's chosen cell, queried directly, match pydggal's
exactly) -- i.e. a pure quantization tie-break at a measure-zero boundary, not
a geometry bug. This is the same singularity class IVEA3H documents (12
xfails: Vert0Lon at EVEN resolutions only, 0/2/4/6/8/10) and IVEA7H documents
(13 xfails: pole + vertex-0 meridian) -- RTEA3H's set is broader (ALL 11
resolutions, not just even) because the RTEA vertex-assignment permutation
(the eC's `rtea` case of `VGCRadialVertex`) warps the cells straddling this
exact meridian differently than IVEA's does; RTEA7H's own forward xfail set
documents the same "RTEA hits it at every resolution where IVEA/ISEA only hit
a subset" pattern for its NorthPole singularity. ISEA3H has zero conformance
xfails at this point set since ISEA does not share IVEA's/RTEA's
pole-longitude sensitivity."""
import json, os, pathlib, pytest
from py4dggs import RTEA3H

_DEFAULT = pathlib.Path(__file__).resolve().parent / "golden" / "rtea3h"  # vendored, see tests/golden/PROVENANCE.md
_DIR = pathlib.Path(os.environ.get("DGGS_GOLDEN_TABLES_RTEA3H", _DEFAULT))
_HAVE = _DIR.exists()
pytestmark = pytest.mark.skipif(not _HAVE, reason=f"RTEA3H golden tables not at {_DIR}; run generate.py RTEA3H")
TOL = 0.0001

# Empirical singularity set (py4dggs.RTEA3H vs the pydggal baseline, 2026-07-03):
# exact-boundary tie-break at the vertex-0 meridian (lon 11.2), ALL resolutions
# 0-10 (broader than IVEA3H's even-only set -- see module docstring). Both
# test_centroid and test_vertices re-run forward quantization (the tables are
# point-keyed, not textId-keyed), so both are affected by the same tie.
XFAIL = {("Vert0Lon", r) for r in range(11)}

def _load(name):
    return json.loads((_DIR / f"{name}.json").read_text())["cases"] if _HAVE else []

def _srt(pairs):
    return sorted((round(a, 4), round(b, 4)) for a, b in pairs)

@pytest.mark.parametrize("c", _load("centroid"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_centroid(c):
    if (c.get("name"), c["res"]) in XFAIL:
        pytest.xfail("RTEA3H exact-singularity boundary tie-break (vertex-0 meridian)")
    z = RTEA3H.zone_from_geo(c["lat"], c["lon"], c["res"])
    assert abs(z.centroid.lat - c["centroidLat"]) < TOL
    assert abs(z.centroid.lon - c["centroidLon"]) < TOL

@pytest.mark.parametrize("c", _load("vertices"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_vertices(c):
    if (c.get("name"), c["res"]) in XFAIL:
        pytest.xfail("RTEA3H exact-singularity boundary tie-break (vertex-0 meridian)")
    z = RTEA3H.zone_from_geo(c["lat"], c["lon"], c["res"])
    dv = [(p.lat, p.lon) for p in z.vertices]
    assert len(dv) == len(c["vertices"])
    assert _srt(dv) == _srt([(a, b) for a, b in c["vertices"]])

@pytest.mark.parametrize("c", _load("textid"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_textid(c):
    # value-keyed -> exact for every cell, no boundary xfail
    val = int(c["value"])
    assert RTEA3H.indexing.to_text(val) == c["textId"]
    assert RTEA3H.indexing.from_text(c["textId"]) == val

@pytest.mark.parametrize("c", _load("hierarchy"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_hierarchy(c):
    val = int(c["value"])
    assert sorted(str(p) for p in RTEA3H.parents(val)) == c["parents"]
    assert sorted(str(cc) for cc in RTEA3H.children(val)) == c["children"]
    cp = RTEA3H.centroid_parent(val)
    assert (str(cp) if cp is not None else None) == c["centroidParent"]
    assert RTEA3H.is_centroid_child(val) == c["isCentroidChild"]

@pytest.mark.parametrize("c", _load("subzones"), ids=lambda c: f"{c.get('name','?')}-r{c['res']}")
def test_subzones(c):
    # value-keyed (A3) -> exact for every cell, no boundary xfail
    val = int(c["value"])
    for depth_str, expected in c["subZonesByDepth"].items():
        assert RTEA3H.sub_zones(val, int(depth_str)) == tuple(int(s) for s in expected)
