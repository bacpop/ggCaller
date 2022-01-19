// ggCaller header
#include "bindings.h"

PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    py::class_<Graph, std::shared_ptr<Graph>>(m, "Graph")
            .def("read", &Graph::read)
            .def("build", &Graph::build)
            .def("data_in", &Graph::in)
            .def("data_out", &Graph::out)
            .def("findORFs", &Graph::findORFs)
            .def("generate_sequence", &Graph::generate_sequence)
            .def("connect_ORFs", &Graph::connect_ORFs)
            .def("generate_clusters", &Graph::generate_clusters)
            .def("refind_gene", &Graph::refind_gene)
            .def("search_graph", &Graph::search_graph)
            .def("node_size", &Graph::node_size);

    m.def("create_graph", []() { return std::make_shared<Graph>(); });

    m.def("get_distances_align", &get_distances_align, "Get distances based on alignment.",
        py::arg("matrix_in"),
        py::arg("no_threads"));

    m.def("get_distances_pa", &get_distances_pa, "Get distances based on feature presence/absence.",
        py::arg("matrix_in"),
        py::arg("no_threads"));
}