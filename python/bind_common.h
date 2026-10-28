#pragma once
// bind_common.h — helpers shared by the prad2py binding translation units.
//
// Bulk arrays are returned as fresh numpy copies: allocate the array, then
// copy into it.  Never wrap a C++ buffer with the 3-argument
// py::array_t(shape, strides, ptr) constructor: without a base handle its
// ownership is ambiguous, and it was observed to corrupt unrelated buffers
// once the owning C++ object was reused or freed.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace py = pybind11;

// Release the GIL for the duration of a bound C++ call; arguments are
// converted before and the result is cast after, with the GIL held.
using release_gil = py::call_guard<py::gil_scoped_release>;

// 1-D numpy copy of p[0..n).
template <class T>
py::array_t<T> to_numpy(const T *p, std::size_t n)
{
    py::array_t<T> arr(static_cast<py::ssize_t>(n));
    std::copy_n(p, n, arr.mutable_data());
    return arr;
}

template <class T>
py::array_t<T> to_numpy(const std::vector<T> &v)
{
    return to_numpy(v.data(), v.size());
}

// Row-major (rows, cols) numpy copy of p[0..rows*cols).
template <class T>
py::array_t<T> to_numpy2d(const T *p, py::ssize_t rows, py::ssize_t cols)
{
    py::array_t<T> arr({rows, cols});
    std::copy_n(p, rows * cols, arr.mutable_data());
    return arr;
}

// Python list of copies of p[0..n).
template <class T>
py::list list_of(const T *p, int n)
{
    py::list out;
    for (int i = 0; i < n; ++i) out.append(p[i]);
    return out;
}

// i if 0 <= i < n, else IndexError(msg).
inline int checked_index(int i, int n, const char *msg)
{
    if (i < 0 || i >= n) throw py::index_error(msg);
    return i;
}

// prad2::database_dir()/daq_config.json.
std::string default_daq_config_path();
