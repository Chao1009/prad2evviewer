// prad2py.cpp — main module glue.
//
// The actual binding code lives in per-area translation units:
//   bind_dec.cpp   → prad2py.dec (evio reader, event data, waveform
//                                   analyzers, EPICS / HV archive)
//   bind_det.cpp   → prad2py.det (GEM + HyCal systems, reconstruction,
//                                   transforms, PipelineBuilder)
//
// Each of those files defines a ``register_XXX(py::module_ &m)`` entry
// point that adds a submodule to the top-level module.
//
// No "do everything" helpers at module root — analyses should drive the
// per-event loop themselves via ``dec.EvChannel.select_event()`` plus
// ``info()``/``fadc()``/``gem()``/``tdc()``/``vtp()`` and accumulate into
// numpy / Python on their own terms.

#include "bind_common.h"

#include "InstallPaths.h"

#include <string>

std::string default_daq_config_path()
{
    return prad2::database_dir() + "/daq_config.json";
}

void register_dec(py::module_ &m);
void register_det(py::module_ &m);

PYBIND11_MODULE(prad2py, m)
{
    m.doc() = "PRad-II (prad2dec + prad2det) Python bindings.";

    m.attr("__version__")    = "0.4.0";

    // Resolved at module-load time: env → module-relative → compile-time.
    // Exposes the same path `default_daq_config()` will use so analyses
    // can point their own lookups at the right place.
    m.attr("DATABASE_DIR")   = prad2::database_dir();

    m.def("default_daq_config", &default_daq_config_path,
          "Return the default daq_config.json path used by analyses.");

    register_dec(m);    // prad2py.dec
    register_det(m);    // prad2py.det
}
