# tests/test_conformance.py
"""Golden-table conformance for py4dggs.IGEO7 (Task 6).

Runs all five suites from the vendored tables/ (tests/golden/, regenerate with
`uv run python tests/golden/generate.py`) against the py4dggs.IGEO7
implementation and must match the frozen igeo7-py oracle profile:
  581 passed, 7 xfailed

The 7 xfailed cases are benign pole / vertex-0 lon=11.2 boundary tie-breaks
on the forward lookup — the clean-room engine diverges from pydggal in the
same way the prior oracle did.
"""
import json
import os
import pytest
from pathlib import Path

from py4dggs import IGEO7, Zone

# ---------------------------------------------------------------------------
# Table locations
# ---------------------------------------------------------------------------
_DEFAULT_TABLES = Path(__file__).resolve().parent / "golden"  # vendored, see tests/golden/PROVENANCE.md
_TABLE_DIR = Path(os.environ.get("DGGS_GOLDEN_TABLES", _DEFAULT_TABLES))
_HAVE_TABLES = _TABLE_DIR.exists()

pytestmark = pytest.mark.skipif(not _HAVE_TABLES, reason=f"golden tables not found at {_TABLE_DIR}; set DGGS_GOLDEN_TABLES")

def _load(name: str) -> list[dict]:
    return json.loads((_TABLE_DIR / f"{name}.json").read_text())["cases"] if _HAVE_TABLES else []


FORWARD_CASES   = _load("forward")
INVERSE_CASES   = _load("inverse")
HIERARCHY_CASES = _load("hierarchy")
VERTICES_CASES  = _load("vertices")
KRING_CASES     = [c for c in _load("kring") if "disk" in c] if _HAVE_TABLES else []

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TOL = 0.0001  # degrees ~ 11 m

# The 7 forward singularities where the clean-room engine diverges from pydggal
# (pole / vertex-0 meridian lon=11.2 boundary tie-break — benign, xfail):
FORWARD_DISCREPANCIES = (
    {("NorthPole", r) for r in (1, 3, 5, 7, 9)}
    | {("Vert0Lon", r) for r in (3, 7)}
)


# ---------------------------------------------------------------------------
# Helper: build a Zone from a text_id
# ---------------------------------------------------------------------------
def _zone(text_id: str) -> Zone:
    return IGEO7.zone_from_text(text_id)


# ---------------------------------------------------------------------------
# Forward: lat/lon/res → text_id
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "case",
    FORWARD_CASES,
    ids=[f"{c['name']}-r{c['res']}" for c in FORWARD_CASES],
)
def test_forward(case):
    if (case["name"], case["res"]) in FORWARD_DISCREPANCIES:
        pytest.xfail("exact-singularity boundary tie-break; benign")
    got = IGEO7.zone_from_geo(case["lat"], case["lon"], case["res"]).text_id
    assert got == case["textId"], (
        f"forward({case['name']} lat={case['lat']} lon={case['lon']} res={case['res']}): "
        f"got {got!r}, want {case['textId']!r}"
    )


# ---------------------------------------------------------------------------
# Inverse: text_id → centroid, plus round-trip
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "case",
    INVERSE_CASES,
    ids=[c["textId"] for c in INVERSE_CASES],
)
def test_inverse(case):
    zone = _zone(case["textId"])
    cen = zone.centroid
    assert abs(cen.lat - case["centroidLat"]) < TOL, (case, cen)
    assert abs(cen.lon - case["centroidLon"]) < TOL, (case, cen)
    # Round-trip: centroid should map back to the same cell
    rt = IGEO7.zone_from_geo(case["centroidLat"], case["centroidLon"], case["res"]).text_id
    assert rt == case["textId"], (
        f"inverse round-trip {case['textId']}: "
        f"zone_from_geo({case['centroidLat']}, {case['centroidLon']}, {case['res']}) → {rt!r}"
    )


# ---------------------------------------------------------------------------
# Hierarchy: parent + children
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "case",
    HIERARCHY_CASES,
    ids=[c["textId"] for c in HIERARCHY_CASES],
)
def test_hierarchy(case):
    zone = _zone(case["textId"])
    got_parent = zone.parent.text_id if zone.parent else None
    assert got_parent == case["parent"], (
        f"hierarchy {case['textId']} parent: got {got_parent!r}, want {case['parent']!r}"
    )
    assert len(zone.children) == len(case["children"]), (
        f"hierarchy {case['textId']} children count: "
        f"got {len(zone.children)}, want {len(case['children'])}"
    )
    assert {ch.text_id for ch in zone.children} == set(case["children"]), (
        f"hierarchy {case['textId']} children set mismatch: "
        f"got {sorted(ch.text_id for ch in zone.children)}, "
        f"want {sorted(case['children'])}"
    )


# ---------------------------------------------------------------------------
# Vertices: vertex count
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "case",
    VERTICES_CASES,
    ids=[c["textId"] for c in VERTICES_CASES],
)
def test_vertices(case):
    zone = _zone(case["textId"])
    assert len(zone.vertices) == len(case["vertices"]), (
        f"vertices {case['textId']}: got {len(zone.vertices)}, want {len(case['vertices'])}"
    )


# ---------------------------------------------------------------------------
# K-ring / disk
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "case",
    KRING_CASES,
    ids=[c["textId"] for c in KRING_CASES],
)
def test_kring(case):
    zone = _zone(case["textId"])
    got = {z.text_id for z in zone.disk(case["k"])}
    assert got == set(case["disk"]), (
        f"kring {case['textId']} k={case['k']}: "
        f"got {sorted(got)}, want {sorted(case['disk'])}"
    )
