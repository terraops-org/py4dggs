# Runnable examples

Standalone scripts, ordered roughly by difficulty. Each one runs on its own with
no arguments and prints what it computed:

```bash
python examples/01_first_zone.py
```

| Script | Shows |
|---|---|
| `01_first_zone.py` | coordinate -> zone, the three representations, resolution |
| `02_compare_grids.py` | the same point on all six grids, and why ids are not portable between them |
| `03_neighbours.py` | the k-ring, mutual adjacency, and the pentagon case |
| `04_hierarchy.py` | Z7's congruent digit hierarchy vs I3H's geometric one |
| `05_tile_store.py` | sub-zones as a fixed-size tile store, and the materialisation bound |
| `06_error_handling.py` | what raises, what returns empty, and what to catch |
| `07_geojson.py` | RFC 7946 export: one cell, a FeatureCollection, and the antimeridian split |

These are executed by the test suite (`tests/test_docs_examples.py`), together
with every ```python block in `README.md` and `TUTORIAL.md` — so an API change
that breaks an example fails CI instead of quietly misleading a reader.

For prose with more context, read `TUTORIAL.md` first; for how the library is put
together, `ARCHITECTURE.md`.
