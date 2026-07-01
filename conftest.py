"""Root-level conftest.py.

Ensures project root and src/ are on sys.path (CI, --cov, editable installs,
unpacked sdists) and registers shared test fixtures via pytest_plugins.
"""

import sys
from pathlib import Path

root = Path(__file__).resolve().parent
for p in (str(root), str(root / "src")):
    if p not in sys.path:
        sys.path.insert(0, p)

pytest_plugins = ["conftest_fixtures"]