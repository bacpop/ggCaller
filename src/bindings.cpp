// ggCaller header
#include "bindings.h"

PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    py::class_<Graph, std::unique_ptr<Graph>>(m, "Graph")
            .def(py::init<>())
            .def("read", &Graph::read)
            .def("build", &Graph::build)
            .def("data_in", &Graph::in)
            .def("data_out", &Graph::out)
            .def("findGenes", &Graph::findGenes)
            .def("generate_sequence", &Graph::generate_sequence)
            .def("refind_gene", &Graph::refind_gene)
            .def("search_graph", &Graph::search_graph)
            .def("node_size", &Graph::node_size)
            .def("rb_hash", &Graph::rb_hash)
            .def("ORF_location", &Graph::ORF_location);

    m.def("get_distances_align", &get_distances_align, "Get distances based on alignment.",
        py::arg("matrix_in"),
        py::arg("no_threads"));

    m.def("get_distances_pa", &get_distances_pa, "Get distances based on feature presence/absence.",
        py::arg("matrix_in"),
        py::arg("no_threads"));

    m.def("read_cluster_file", &read_cluster_file, "Read cluster map file",
        py::arg("cluster_file"));

    m.def("clear_graph", &clear_graph, "Call graph destructor",
        py::arg("graph"));

//    m.def("use_count", &use_count, "Get pointer count",
//        py::arg("graph_ptr"));
}