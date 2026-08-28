"""The documentation's code must actually run.

`README.md` and `TUTORIAL.md` teach the library through worked examples, and
`examples/` ships them as standalone scripts -- but until now nothing executed
any of it, so an API change could silently invalidate every example while the
suite stayed green. This repo already learned that lesson once from the golden
tables (see test_golden_tables_present.py): a doc example that no longer runs is
the same silent degradation, just aimed at readers instead of CI.

The doc blocks are executed per FILE, in order, in one shared namespace, because
the tutorial builds on itself (e.g. its pentagon example reuses an earlier
import). A block that raises fails the file's test and names the block.
"""
import pathlib
import re
import runpy
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"

_BLOCK_RE = re.compile(r"^```python\n(.*?)^```", re.S | re.M)


def _python_blocks(md: pathlib.Path) -> list[str]:
    return _BLOCK_RE.findall(md.read_text())


@pytest.mark.parametrize("doc", ["README.md", "TUTORIAL.md"])
def test_documentation_examples_run(doc):
    path = ROOT / doc
    assert path.is_file(), f"missing {path}"
    blocks = _python_blocks(path)
    assert blocks, f"{doc} has no ```python blocks — did the fences change?"
    ns: dict = {}
    for i, block in enumerate(blocks, 1):
        try:
            exec(compile(block, f"<{doc} block {i}>", "exec"), ns)
        except Exception as e:  # noqa: BLE001 -- we want the block's identity in the message
            first = block.strip().splitlines()[0]
            pytest.fail(
                f"{doc} block {i} failed ({type(e).__name__}: {e})\n"
                f"  first line: {first}"
            )


def test_examples_directory_exists():
    assert EXAMPLES.is_dir(), (
        f"missing {EXAMPLES} — the runnable examples are part of the docs "
        f"surface, not optional extras"
    )
    scripts = sorted(EXAMPLES.glob("*.py"))
    assert scripts, f"no example scripts in {EXAMPLES}"


def _example_scripts():
    return sorted(p.name for p in EXAMPLES.glob("*.py")) if EXAMPLES.is_dir() else []


@pytest.mark.parametrize("script", _example_scripts())
def test_example_script_runs(script):
    """Run each example as a real script, in a subprocess, exactly as a reader
    would (`python examples/foo.py`) -- so a missing `__main__` guard or an
    import that only works from the repo root is caught."""
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / script)],
        capture_output=True, text=True, timeout=120, cwd=ROOT,
    )
    assert proc.returncode == 0, (
        f"examples/{script} exited {proc.returncode}\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )
    assert proc.stdout.strip(), f"examples/{script} printed nothing — an example should show its result"
