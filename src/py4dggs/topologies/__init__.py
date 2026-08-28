# src/py4dggs/topologies/__init__.py
"""Topology implementations (cell tiling + quantization) for the DGGS reference.

A Topology owns the aperture and cell shape: it turns a projection's planar point
into a (base cell, direction-digit) address and back into planar cell geometry.
Each concrete grid plugs in one Topology (see ``py4dggs.interfaces.Topology``)."""
