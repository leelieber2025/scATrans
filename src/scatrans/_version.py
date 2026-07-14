"""Package version — single source of truth.

Bump ``__version__`` here when preparing a release. Packaging
(``pyproject.toml`` dynamic version), ``scatrans.__version__``, docs
``release``, and the release zip helper all read this module.

Do not regenerate this file from setuptools_scm; it is intentionally
committed.
"""

from __future__ import annotations

__all__ = ["__version__", "version", "version_tuple"]

__version__ = "0.10.3"
version = __version__

# Parsed once for callers that need a comparable tuple (major, minor, patch, …).
_parts: list[int | str] = []
for _seg in __version__.split("."):
    try:
        _parts.append(int(_seg))
    except ValueError:
        _parts.append(_seg)
version_tuple: tuple[int | str, ...] = tuple(_parts)
