"""Guard: the vendored golden tables must actually be there.

Every `test_*_conformance.py` module carries a module-level `skipif` on "tables
found?", and its cases are *parametrized from the table files* — so if the
tables go missing those suites generate ZERO cases and vanish without even
showing up as skips. That is exactly how CI once reported a green 247 tests
where a full run is 4209.

This module therefore has NO skipif and reads no tables: it just asserts the
fixtures exist, so losing them is a loud failure instead of a quiet 95% drop in
coverage. See tests/golden/PROVENANCE.md.
"""
import json
import pathlib

import pytest

GOLDEN = pathlib.Path(__file__).resolve().parent / "golden"

# The IGEO7 (ISEA7H_Z7) tables live at the top level; every other grid gets a
# subdirectory. Keep in step with tests/golden/PROVENANCE.md.
IGEO7_TABLES = ["forward.json", "inverse.json", "vertices.json", "kring.json", "hierarchy.json"]
A3_TABLES = ["centroid.json", "vertices.json", "neighbours.json", "textid.json",
             "hierarchy.json", "subzones.json"]
A7_TABLES = IGEO7_TABLES  # the A7 variants carry the same five tables
A7_VARIANTS = ["ivea7h", "rtea7h"]
A3_VARIANTS = ["isea3h", "ivea3h", "rtea3h"]


def test_golden_directory_exists():
    assert GOLDEN.is_dir(), (
        f"vendored golden tables missing at {GOLDEN} — the conformance suites "
        f"would silently collect zero cases. See tests/golden/PROVENANCE.md."
    )


def test_generator_is_vendored_too():
    """The repo must be able to regenerate its own fixtures without the spec repo."""
    assert (GOLDEN / "generate.py").is_file()
    assert (GOLDEN / "PROVENANCE.md").is_file()


def _assert_table_has_cases(p: pathlib.Path, label: str) -> None:
    """A table is only useful if it has CASES.

    Every conformance module loads `data["cases"]` and parametrizes from it, so a
    file that is present, valid JSON and truthy -- but whose `cases` list is
    empty -- still collects zero tests. A bare `assert data` waved those through,
    reopening the same silent-degradation hole this module exists to close.
    """
    assert p.is_file(), f"missing {label} golden table {p}"
    data = json.loads(p.read_text())
    assert data, f"empty {label} golden table {p}"
    assert isinstance(data, dict) and "cases" in data, (
        f"{label} golden table {p} has no 'cases' key — the conformance suites "
        f"read data['cases'] and would collect zero tests"
    )
    assert data["cases"], (
        f"{label} golden table {p} has an EMPTY 'cases' list — the conformance "
        f"suite parametrized from it would silently vanish"
    )


@pytest.mark.parametrize("name", IGEO7_TABLES)
def test_igeo7_tables_present_and_non_empty(name):
    _assert_table_has_cases(GOLDEN / name, "IGEO7")


@pytest.mark.parametrize("grid", A7_VARIANTS + A3_VARIANTS)
def test_variant_grid_directories_present(grid):
    d = GOLDEN / grid
    assert d.is_dir(), f"missing golden tables for {grid} at {d}"
    assert any(d.glob("*.json")), f"no tables in {d}"


@pytest.mark.parametrize("grid", A7_VARIANTS)
@pytest.mark.parametrize("name", A7_TABLES)
def test_a7_variant_tables_present_and_non_empty(grid, name):
    """Per-FILE, like the A3 variants below.

    `test_variant_grid_directories_present` only asks whether the directory holds
    *any* .json, so four of the five tables could vanish and the guard would stay
    green -- while four conformance suites quietly collected nothing.
    """
    _assert_table_has_cases(GOLDEN / grid / name, grid)


@pytest.mark.parametrize("grid", A3_VARIANTS)
@pytest.mark.parametrize("name", A3_TABLES)
def test_a3_tables_present_and_non_empty(grid, name):
    _assert_table_has_cases(GOLDEN / grid / name, grid)


def test_no_dependency_on_the_sibling_spec_repo():
    """Regression: the tables used to be read from ../igeo7-spec, a PRIVATE repo.

    A fresh clone must be self-sufficient, so no test module may reach outside
    the repo for fixtures by default.
    """
    tests_dir = pathlib.Path(__file__).resolve().parent
    me = pathlib.Path(__file__).resolve()
    needle = "igeo7-spec"
    offenders = [
        p.name
        for p in sorted(tests_dir.glob("test_*.py"))
        if p.resolve() != me                      # this file names the path it forbids
        and f'parents[2] / "{needle}"' in p.read_text()
    ]
    assert not offenders, (
        f"these modules still default to the sibling igeo7-spec repo: {offenders}"
    )
