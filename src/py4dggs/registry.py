# src/py4dggs/registry.py
"""IGEO7 singleton and grid registry."""
from py4dggs.grid import Grid
from py4dggs.types import GridConfig
from py4dggs.projections.isea import ISEAProjection
from py4dggs.projections.ivea import IVEAProjection
from py4dggs.projections.rtea import RTEAProjection
from py4dggs.topologies.hex_a7 import HexAperture7Topology
from py4dggs.topologies.hex_a3 import HexAperture3Topology
from py4dggs.indexings.z7 import Z7Indexing
from py4dggs.indexings.i3h import I3HIndexing

IGEO7 = Grid(ISEAProjection(), HexAperture7Topology(), Z7Indexing(), GridConfig(), "IGEO7")
IVEA7H = Grid(IVEAProjection(), HexAperture7Topology(), Z7Indexing(), GridConfig(), "IVEA7H")
RTEA7H = Grid(RTEAProjection(), HexAperture7Topology(), Z7Indexing(), GridConfig(), "RTEA7H")
ISEA3H = Grid(ISEAProjection(), HexAperture3Topology(), I3HIndexing(), GridConfig(), "ISEA3H")
IVEA3H = Grid(IVEAProjection(), HexAperture3Topology(), I3HIndexing(), GridConfig(), "IVEA3H")
RTEA3H = Grid(RTEAProjection(), HexAperture3Topology(), I3HIndexing(), GridConfig(), "RTEA3H")
_GRIDS = {"IGEO7": IGEO7, "IVEA7H": IVEA7H, "RTEA7H": RTEA7H, "ISEA3H": ISEA3H, "IVEA3H": IVEA3H, "RTEA3H": RTEA3H}

def get_grid(name: str) -> Grid:
    if name not in _GRIDS:
        raise KeyError(f"unknown grid {name!r}; known: {sorted(_GRIDS)}")
    return _GRIDS[name]
