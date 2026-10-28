// InstallPaths.cpp — implementation of the run-time data-directory resolver.
//
// See InstallPaths.h for the lookup policy.  module_dir() uses dladdr() on
// our own symbol on Linux/macOS and GetModuleHandleExW(..._FROM_ADDRESS) +
// GetModuleFileNameW on Windows.  (Not /proc/self/exe: inside the prad2py
// module that names the Python interpreter, not the .so.)
//
// Needs C++17 <filesystem>.  Link ${CMAKE_DL_LIBS} on Linux (empty on
// glibc ≥ 2.34, -ldl on older systems) — handled by prad2dec's
// CMakeLists.
//=============================================================================

#include "InstallPaths.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <system_error>

#if defined(_WIN32)
  #ifndef WIN32_LEAN_AND_MEAN
  #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#else
  #include <dlfcn.h>
#endif

// Set by the top-level CMake project; "." for a standalone prad2dec build.
#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

namespace prad2 {

std::string module_dir()
{
#if defined(_WIN32)
    HMODULE h = nullptr;
    // Use the address of this function as the anchor — whichever DLL /
    // exe embeds the prad2dec static library will own it.
    if (!GetModuleHandleExW(
            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            reinterpret_cast<LPCWSTR>(&module_dir), &h)) {
        return {};
    }
    wchar_t buf[MAX_PATH];
    DWORD n = GetModuleFileNameW(h, buf, MAX_PATH);
    if (n == 0 || n >= MAX_PATH) return {};
    std::wstring ws(buf, n);
    std::error_code ec;
    fs::path p = fs::weakly_canonical(fs::path(ws), ec);
    if (ec) p = fs::path(ws);
    return p.parent_path().string();
#else
    Dl_info info{};
    if (!dladdr(reinterpret_cast<const void *>(&module_dir), &info) ||
        !info.dli_fname) {
        return {};
    }
    std::error_code ec;
    fs::path p = fs::weakly_canonical(fs::path(info.dli_fname), ec);
    if (ec) p = fs::path(info.dli_fname);
    return p.parent_path().string();
#endif
}

std::string resolve_data_dir(const char *env_name,
                             std::initializer_list<const char *> rel_candidates,
                             const char *compile_default)
{
    if (env_name) {
        if (const char *env = std::getenv(env_name); env && *env) {
            return env;
        }
    }

    std::string base = module_dir();
    if (!base.empty()) {
        std::error_code ec;
        for (const char *rel : rel_candidates) {
            if (!rel) continue;
            fs::path cand = fs::path(base) / rel;
            fs::path norm = fs::weakly_canonical(cand, ec);
            if (ec) { ec.clear(); norm = cand; }
            if (fs::is_directory(norm, ec)) {
                return norm.string();
            }
            ec.clear();
        }
    }

    return compile_default ? std::string(compile_default) : std::string();
}

std::string database_dir()
{
    return resolve_data_dir("PRAD2_DATABASE_DIR",
                            {"../share/prad2evviewer/database",
                             "../../share/prad2evviewer/database"},
                            DATABASE_DIR);
}

std::string find_database_file(const std::string &name)
{
    auto readable = [](const std::string &p) { return std::ifstream(p).good(); };
    const std::string db = database_dir();
    if (!db.empty() && readable(db + "/" + name)) return db + "/" + name;
    for (const std::string &p : {name, "database/" + name, "../database/" + name})
        if (readable(p)) return p;
    return {};
}

bool is_absolute_path(const std::string &path)
{
    if (path.empty()) return false;
    if (path[0] == '/' || path[0] == '\\') return true;
    return path.size() >= 2 && path[1] == ':';   // Windows drive letter
}

std::string resolve_db_path(const std::string &path, const std::string &db_dir)
{
    if (path.empty() || db_dir.empty() || is_absolute_path(path)) return path;
    return db_dir + "/" + path;
}

} // namespace prad2
