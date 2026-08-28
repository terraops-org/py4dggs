"""ISEA projection variant. The shared icosahedral equal-area kernel lives in
:mod:`py4dggs.projections.icovertex`; this file binds it to the ISEA vertex assignment
(the eC ``isea`` case of ``VGCRadialVertex`` — ``icoVertexGreatCircle.ec:127``)."""
from py4dggs.projections.icovertex import _IcoVertexProjection


class ISEAProjection(_IcoVertexProjection):
    """ISEA (Icosahedral Snyder Equal-Area) = the shared kernel with ``radial_vertex="isea"``."""

    radial_vertex = "isea"
