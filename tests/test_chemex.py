"""Unit and regression test for the chemex package."""

import chemex
import sys


def test_chemex_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "chemex" in sys.modules
