"""Start-up helpers for the scripts/ and gem/ Python tools: locating
prad2py, the database directory and PyQt6's bundled Qt libraries.

Qt-free on purpose: fix_qt_lib_path() has to run before PyQt6 is imported,
and the command-line tools use the rest without Qt.  gem/ tools import
this module through their ../scripts sys.path entry, which works in the
source tree and under <prefix>/share/prad2evviewer alike.
"""
from __future__ import annotations

import os
import site
import sys
from pathlib import Path
from typing import Optional, Tuple

# The repository root in a checkout, <prefix>/share/prad2evviewer when
# installed.
_ROOT_DIR = Path(__file__).resolve().parent.parent

# Probed in this order; each hit is put first on sys.path, so the last
# existing one wins.
_BUILD_SUBDIRS = ("build/python", "build-release/python",
                  "build/Release/python")

PRAD2PY_HINT = ("Build it with:\n"
                "    cmake -DBUILD_PYTHON=ON -S . -B build && "
                "cmake --build build\n"
                "and add build/python/ to PYTHONPATH.")


def import_prad2py(build_first: bool = True) -> Tuple[Optional[object], str]:
    """Import prad2py and return ``(module, "")``, or ``(None, error)``
    with ``error`` as ``"<ExceptionType>: <message>"``.

    The Python build directories of this checkout take precedence over
    PYTHONPATH, or with ``build_first`` False are tried only when a plain
    import fails.  An installed tree has no build directories, so there
    the wrapper's PYTHONPATH decides.
    """
    if not build_first:
        try:
            import prad2py
            return prad2py, ""
        except ImportError:
            pass
    for sub in _BUILD_SUBDIRS:
        cand = _ROOT_DIR / sub
        if cand.is_dir() and str(cand) not in sys.path:
            sys.path.insert(0, str(cand))
    try:
        import prad2py
        return prad2py, ""
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def database_dir(use_env: bool = True) -> Path:
    """$PRAD2_DATABASE_DIR when set (and ``use_env``), else the database/
    directory next to scripts/."""
    env = os.environ.get("PRAD2_DATABASE_DIR") if use_env else None
    return Path(env).resolve() if env else (_ROOT_DIR / "database").resolve()


def find_database_file(name: str, use_env: bool = True) -> Optional[Path]:
    """First existing file among $PRAD2_DATABASE_DIR/<name> (with
    ``use_env``), <root>/database/<name>, ./database/<name> and ./<name>,
    or None.  ``name`` may contain subdirectories (``runinfo/general.json``).
    """
    dirs = []
    env = os.environ.get("PRAD2_DATABASE_DIR") if use_env else None
    if env:
        dirs.append(Path(env))
    dirs += [_ROOT_DIR / "database", Path.cwd() / "database", Path.cwd()]
    for d in dirs:
        p = d / name
        if p.is_file():
            return p.resolve()
    return None


def fix_qt_lib_path() -> None:
    """Re-exec the interpreter with PyQt6's bundled Qt6/lib first on
    LD_LIBRARY_PATH.  Call it before importing PyQt6.

    On some systems the system libQt6DBus.so.6 is built against another Qt
    than PyQt6's bundled Qt6Core, and importing PyQt6 fails with an
    undefined Qt_6_PRIVATE_API symbol.  With the bundled directory first
    the dynamic linker takes every Qt library from one build.
    LD_LIBRARY_PATH is only read at process start, hence the re-exec.  A
    no-op when there is no bundled Qt or it is already on the path.
    """
    dirs = []
    try:
        dirs += site.getsitepackages()
    except AttributeError:
        pass
    try:
        dirs.append(site.getusersitepackages())
    except Exception:  # noqa: BLE001
        pass
    for sp in dirs:
        qt6_lib = os.path.join(sp, "PyQt6", "Qt6", "lib")
        if os.path.isdir(qt6_lib):
            cur = os.environ.get("LD_LIBRARY_PATH", "")
            if qt6_lib not in cur.split(":"):
                env = dict(os.environ,
                           LD_LIBRARY_PATH=qt6_lib + (":" + cur if cur else ""))
                os.execve(sys.executable, [sys.executable] + sys.argv, env)
            return
