"""Shared pytest fixtures for scATrans test suite.

Fixtures live in ``conftest_fixtures.py`` at the project root and are loaded
via pytest_plugins so CI reliably registers them.
"""

pytest_plugins = ["conftest_fixtures"]
