from __future__ import annotations

import importlib.util
import types

from tests.conftest import REPO_DIR


def load_module_from_repo_path(*, rel_path: str, module_name: str) -> types.ModuleType:
    """Load a Python module from a repository-relative path.

    Useful for top-level scripts (e.g. ``dicast.py``) that are not importable as
    part of the ``dicast_lib`` package. Library code under ``dicast_lib/`` should
    be imported directly (``from dicast_lib import utils``).
    """
    script_path = REPO_DIR / rel_path
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
