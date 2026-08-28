"""RTEA projection variant. The shared icosahedral equal-area kernel lives in
:mod:`py4dggs.projections.icovertex`; this file binds it to the RTEA vertex
assignment (the eC ``rtea`` case of ``VGCRadialVertex`` — a permutation of the
isea/ivea ``AB``/``AC``/``BC`` edge angles, ``icoVertexGreatCircle.ec``)."""
from py4dggs.projections.icovertex import _IcoVertexProjection


class RTEAProjection(_IcoVertexProjection):
    radial_vertex = "rtea"
