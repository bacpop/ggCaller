// ggCaller header
#include "bindings.h"

PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    py::class_<Graph, std::shared_ptr<Graph>>(m, "Graph")
            .def(py::init<>())
            .def("read", &Graph::read)
            .def("build", &Graph::build)
            .def("findORFs", &Graph::findORFs)
            .def("generate_sequence", &Graph::generate_sequence)
            .def("connect_ORFs", &Graph::connect_ORFs)
            .def("generate_clusters", &Graph::generate_clusters)
            .def("refind_gene", &Graph::refind_gene);
}