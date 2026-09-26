#pragma once
//=============================================================================
// InstallPaths.h — runtime data-directory resolution for PRad-II binaries.
//
// Every shipped executable (prad2_server, prad2evviewer, …) and the
// prad2py Python module need to locate ``database/`` and ``resources/``
// at run time.  The preference order is:
//
//     1. The ``PRAD2_DATABASE_DIR`` / ``PRAD2_RESOURCE_DIR`` env var,
//        set verbatim (what prad2_setup.sh + the bin/ wrappers configure).
//     2. A path relative to *this module*, resolved via dladdr on POSIX
//        and ``GetModuleHandleExW`` on Windows.  Makes installed binaries
//        relocatable — move the install tree and they still find data.
//     3. The build-time ``DATABASE_DIR`` / ``RESOURCE_DIR`` constants
//        (last-resort fallback for dev-in-tree use).
//
// The "module" for a static-linked executable is the executable itself;
// for the prad2py extension module it is the .so / .pyd file.  Using
// dladdr (and its Windows equivalent) rather than ``/proc/self/exe``
// directly is what makes the helper work for both — /proc/self/exe for
// a Python process running ``import prad2py`` would point at the Python
// interpreter, not at our module.
//=============================================================================

#include <initializer_list>
#include <string>

namespace prad2 {

/// Absolute path to the directory containing the DLL / shared object /
/// executable that defines ``module_dir``.  Returns an empty string
/// if runtime resolution fails (very rare — exhausted platform APIs).
std::string module_dir();

/// Resolve a PRad-II data directory.
///
/// @param env_name        Override env var (e.g. ``PRAD2_DATABASE_DIR``).
///                        If set, its value is returned unchanged.
/// @param rel_candidates  Relative paths from the module directory to try
///                        in order.  The first one that exists wins.  For
///                        example, an executable installed to
///                        ``<prefix>/bin/`` reaches the database via
///                        ``"../share/prad2evviewer/database"``; the
///                        prad2py extension at ``<prefix>/lib/prad2py/``
///                        reaches it via ``"../../share/prad2evviewer/database"``.
/// @param compile_default Last-resort fallback (typically the ``DATABASE_DIR``
///                        preprocessor constant).  May be null — in that
///                        case an empty string is returned if lookup fails.
std::string resolve_data_dir(const char *env_name,
                             std::initializer_list<const char *> rel_candidates,
                             const char *compile_default);

/// The database directory: ``resolve_data_dir("PRAD2_DATABASE_DIR",
/// {"../share/prad2evviewer/database", "../../share/prad2evviewer/database"},
/// DATABASE_DIR)`` with the ``DATABASE_DIR`` prad2dec was built with, so it
/// resolves from both ``<prefix>/bin/`` and ``<prefix>/lib/prad2py/``.
std::string database_dir();

/// A database file for the command-line tools: ``database_dir()/name`` if it
/// can be opened, else the first of ``name``, ``database/name`` and
/// ``../database/name`` (working-directory relative, returned as written)
/// that can be opened.  Empty if none can.
std::string find_database_file(const std::string &name);

/// True for ``/…``, ``\\…`` and Windows drive paths (``C:…``).
bool is_absolute_path(const std::string &path);

/// A database-relative config path joined to ``db_dir`` (``db_dir/path``).
/// Empty or absolute paths, and any path when ``db_dir`` is empty, are
/// returned unchanged.
std::string resolve_db_path(const std::string &path, const std::string &db_dir);

} // namespace prad2
