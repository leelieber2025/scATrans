#!/usr/bin/env python3
"""Build a single clean GitHub release zip (source + docs + tests, no build artifacts).

Version is taken from ``src/scatrans/_version.py`` (single source of truth).
Optionally syncs ``CITATION.cff`` and ``packaging/ecosystem-packages/meta.yaml``
to match before packaging.
"""

from __future__ import annotations

import ast
import fnmatch
import re
import zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAMP = date.today().isoformat()

# Patterns applied only to paths collected via RELEASE_INCLUDES.
RELEASE_EXCLUDES = [
    "backup/**",
    "build/**",
    "dist/**",
    "dist-github-release/**",
    "scatrans.egg-info/**",
    "src/scatrans.egg-info/**",
    "**/__pycache__/**",
    "**/*.pyc",
    "**/.pytest_cache/**",
    "**/.ruff_cache/**",
    "**/.mypy_cache/**",
    "**/.ipynb_checkpoints/**",
    "docs/_build/**",
    "docs/api/generated/**",
    ".coverage",
    "coverage.xml",
    "coverage-*.xml",
    "htmlcov/**",
    ".claude/**",
    ".git/**",
    "**/*.h5ad",
    "**/*.zip",
]

RELEASE_INCLUDES = [
    ".github/**",
    ".gitignore",
    ".pre-commit-config.yaml",
    ".readthedocs.yaml",
    "CHANGELOG.md",
    "CITATION.cff",
    "LICENSE",
    "MANIFEST.in",
    "README.md",
    "SECURITY.md",
    "conftest.py",
    "conftest_fixtures.py",
    "docs/**",
    "examples/**",
    "packaging/**",
    "pyproject.toml",
    "setup.cfg",
    "scripts/make_release_zips.py",
    "src/scatrans/**",
    "tests/**",
]


def _read_version() -> str:
    path = ROOT / "src" / "scatrans" / "_version.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if (
                isinstance(target, ast.Name)
                and target.id == "__version__"
                and isinstance(node.value, ast.Constant)
                and isinstance(node.value.value, str)
            ):
                return node.value.value
    raise SystemExit(f"Could not parse __version__ from {path}")


def _sync_metadata(version: str) -> None:
    """Keep auxiliary packaging metadata aligned with _version.py."""
    citation = ROOT / "CITATION.cff"
    if citation.is_file():
        text = citation.read_text(encoding="utf-8")
        new = re.sub(r"(?m)^version:\s*.*$", f"version: {version}", text, count=1)
        if new != text:
            citation.write_text(new, encoding="utf-8")
            print(f"Synced CITATION.cff → version: {version}")

    eco = ROOT / "packaging" / "ecosystem-packages" / "meta.yaml"
    if eco.is_file():
        text = eco.read_text(encoding="utf-8")
        new = re.sub(r"(?m)^version:\s*.*$", f"version: {version}", text, count=1)
        if new != text:
            eco.write_text(new, encoding="utf-8")
            print(f"Synced packaging/ecosystem-packages/meta.yaml → version: {version}")


def _matches(path: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(path, pat) for pat in patterns)


def _expand_include(pattern: str) -> list[Path]:
    if pattern.endswith("/**"):
        base = ROOT / pattern[:-3]
        if not base.exists():
            return []
        return [p for p in base.rglob("*") if p.is_file()]
    path = ROOT / pattern
    return [path] if path.is_file() else []


def _collect_files(include_patterns: list[str], exclude_patterns: list[str]) -> list[Path]:
    seen: set[Path] = set()
    for pattern in include_patterns:
        for path in _expand_include(pattern):
            rel = path.relative_to(ROOT).as_posix()
            if _matches(rel, exclude_patterns):
                continue
            seen.add(path)
    return sorted(seen)


def _write_zip(paths: list[Path], out_path: Path, *, arc_root: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()
    with zipfile.ZipFile(out_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in paths:
            rel = path.relative_to(ROOT).as_posix()
            zf.write(path, f"{arc_root}/{rel}")
    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"Wrote {out_path.name}: {len(paths)} files, {size_mb:.2f} MiB")


def main() -> None:
    version = _read_version()
    _sync_metadata(version)

    files = _collect_files(RELEASE_INCLUDES, RELEASE_EXCLUDES)
    required = {
        "pyproject.toml",
        "README.md",
        "CHANGELOG.md",
        "LICENSE",
        "src/scatrans/__init__.py",
        "src/scatrans/_version.py",
        "src/scatrans/pl.py",
        "docs/conf.py",
        "docs/user_guide/plotting.md",
        "docs/api/index.md",
        ".readthedocs.yaml",
        "tests/test_pl_smoke.py",
    }
    names = {p.relative_to(ROOT).as_posix() for p in files}
    missing = sorted(required - names)
    if missing:
        raise SystemExit(f"Release zip missing required files: {missing}")

    arc_root = f"scatrans-{version}"
    out_zip = ROOT / f"scatrans-{version}-github-release-{STAMP}.zip"
    _write_zip(files, out_zip, arc_root=arc_root)

    print(f"Version: {version} (from src/scatrans/_version.py)")
    print(f"Single clean release package ready: {out_zip.name}")
    print("Upload to: https://github.com/leelieber2025/scATrans/releases")


if __name__ == "__main__":
    main()
