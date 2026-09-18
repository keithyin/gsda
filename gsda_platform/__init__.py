"""gsda_platform — the control-plane (Platform Server) for the analysis platform.

Named ``gsda_platform`` (NOT ``platform``) on purpose: a top-level ``platform``
package shadows the Python stdlib ``platform`` module, which breaks ``import fastapi``
(pydantic's import chain dies). The module layout matches §47 of the design doc.
"""

__version__ = "0.1.0"
