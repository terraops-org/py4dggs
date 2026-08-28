"""IVEA projection variant. The shared icosahedral equal-area kernel lives in
:mod:`py4dggs.projections.icovertex`; this file binds it to the IVEA vertex
assignment (the eC ``ivea`` case of ``VGCRadialVertex`` — ``B``/``C``,
``beta``/``gamma`` and ``AB``/``AC`` swapped vs ISEA; ``icoVertexGreatCircle.ec``)."""
from py4dggs.projections.icovertex import _IcoVertexProjection


class IVEAProjection(_IcoVertexProjection):
    radial_vertex = "ivea"
