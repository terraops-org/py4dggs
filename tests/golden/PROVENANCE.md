# Golden conformance tables — provenance

Frozen fixtures that pin every grid's behaviour to **DGGAL** (via `pydggal`).
They are the second half of this project's verification strategy: the
`test_*_fuzz.py` suites check us against a *live* pydggal, these tables check us
against a *recorded* pydggal that cannot drift under us.

## Where these came from

Copied (not moved) from the multi-implementation spec repo:

| | |
|---|---|
| Source repo | `terraops-org/igeo7-spec` (private) |
| Source path | `verify/tables/` and `verify/generate.py` |
| Source commit | `82caa2a235add3125253d96d451c6a18df679b83` |
| Commit date | 2026-07-06 |
| Vendored on | 2026-07-29 |
| Oracle used | `dggal==0.0.6` (`ecrt==0.0.6`) |

`igeo7-spec` remains the generator of record and serves other implementations
(there is a JS runner there verifying against these same tables). It was **not**
modified by this vendoring.

## Why they are vendored

`igeo7-spec` is private, so CI could not read it without a credential — and the
failure was silent: the conformance tests are *parametrized from these files*,
so with the tables absent they generate **zero** cases and never even appear as
skips. CI reported ~247 tests instead of 4209 and still showed green.

Vendoring also finishes what `CLAUDE.md` already asks for — "`to.py4dggs-py`
itself has zero sibling-repo dependencies" — of which these tables were the last
exception. A fresh `git clone` now runs the entire suite with no sibling repo,
no token, and no configuration.

## Regenerating

Needs the `dggal` dev dependency (already in `uv sync`). Run from the repo root:

```bash
uv run python tests/golden/generate.py            # IGEO7  -> tests/golden/*.json
uv run python tests/golden/generate.py IVEA7H_Z7  # -> tests/golden/ivea7h/
uv run python tests/golden/generate.py RTEA7H_Z7  # -> tests/golden/rtea7h/
uv run python tests/golden/generate.py ISEA3H     # -> tests/golden/isea3h/
uv run python tests/golden/generate.py IVEA3H     # -> tests/golden/ivea3h/
uv run python tests/golden/generate.py RTEA3H     # -> tests/golden/rtea3h/
```

Regenerating is a **deliberate** act — it re-freezes the baseline against
whatever `dggal` is installed. Do it only when adding a grid or intentionally
moving to a new DGGAL version, and say which in the commit message. A diff in
these files is otherwise a red flag, not a routine update.

## The one divergence from upstream

`generate.py` is byte-identical to `igeo7-spec`'s copy except for the `TABLES`
path constant: upstream writes to `<script dir>/tables`, but here the tables ARE
the script's directory, so the extra level is dropped. The change is marked
in-file with a `VENDORED DIVERGENCE` comment. Everything else — the test point
set, rounding, table schemas — is untouched, so the two copies stay diffable.

## Layout

```
tests/golden/
  generate.py         the generator (only needs dggal + stdlib)
  PROVENANCE.md       this file
  forward.json        IGEO7 (ISEA7H_Z7) — the default grid, tables at top level
  inverse.json
  vertices.json
  kring.json
  hierarchy.json
  ivea7h/  rtea7h/    aperture-7 projection variants
  isea3h/  ivea3h/  rtea3h/    aperture-3 family (centroid/vertices/neighbours/
                               textid/hierarchy/subzones)
```

## Overriding the path

Each conformance module defaults to the directory above and still honours an
env override, so an out-of-tree table set can be pointed at without editing code:

```
DGGS_GOLDEN_TABLES            IGEO7 (tables root)
DGGS_GOLDEN_TABLES_IVEA7H     DGGS_GOLDEN_TABLES_RTEA7H
DGGS_GOLDEN_TABLES_ISEA3H     DGGS_GOLDEN_TABLES_IVEA3H     DGGS_GOLDEN_TABLES_RTEA3H
```

If a path does not exist the affected module skips rather than fails — which is
why a *silent* drop in test count is the thing to watch for.
