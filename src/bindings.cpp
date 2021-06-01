// ggCaller header
#include "bindings.h"

PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    py::class_<Graph>(m, "Graph")
            .def("build", py::overload_cast<(const std::string&,
                                           const int,
                                           const std::vector<std::string>&,
                                           const std::vector<std::string>&,
                                           size_t,
                                           bool,
                                           const bool,
                                           const std::string&)> (&Graph::build), "Build graph from list.")
            .def("build", py::overload_cast<(const std::string&,
                                            const std::string&,
                                            const std::vector<std::string>&,
                                            const std::vector<std::string>&,
                                            size_t,
                                            const bool)> (&Graph::build), "Build graph from existing GFA and colours file.")
            .def("findORFs", &Graph::findORFs);

}