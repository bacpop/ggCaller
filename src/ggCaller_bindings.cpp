// ggCaller header
#include "ggCaller_classes.h"

std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyunitigMap, size_t> py_ggCaller_graphexists (const std::string& graphfile,
                                                                                       const std::string& coloursfile,
                                                                                       const std::vector<std::string>& start_codons,
                                                                                       const std::vector<std::string>& stop_codons_for,
                                                                                       const std::vector<std::string>& stop_codons_rev,
                                                                                       size_t num_threads,
                                                                                       const bool is_ref,
                                                                                       const bool write_idx,
                                                                                       const bool repeat,
                                                                                       const bool no_filter,
                                                                                       const size_t max_path_length,
                                                                                       const size_t min_ORF_length,
                                                                                       const size_t max_ORF_overlap) {
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Reading coloured compacted DBG..." << endl;

    // read in graph
    ColoredCDBG<> ccdbg;
    ccdbg.read(graphfile, coloursfile, num_threads);

    // set local variables
    const int kmer = ccdbg.getK();
    const int overlap = kmer - 1;
    omp_set_num_threads(num_threads);

    // iterate over unitigs, search for presence of sub-sequence
    cout << "Constructing unitig matrices for graph traversal..." << endl;

    // get the number of colours, generate empty colour vector
    const size_t nb_colours = ccdbg.getNbColors();
    const std::vector<bool> empty_colour_arr(nb_colours, 0);

    // get colour names
    std::vector<std::string> input_colours = ccdbg.getColorNames();

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    auto graph_tuple = std::move(index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));

    // clear ccdbg to free memory
    ccdbg.clear();

    // generate complete paths
    cout << "Generating complete stop-stop paths..." << endl;
    auto path_pair = std::move(traverse_graph(graph_tuple, repeat, empty_colour_arr, max_path_length));

    // generate ORF sequences
    cout << "Generating ORF sequences from complete paths..." << endl;
    auto ORF_pair = std::move(call_ORFs(path_pair, std::get<0>(graph_tuple), stop_codons_for, start_codons, overlap, min_ORF_length));

    // clear path_pair to free memory
    path_pair.first.clear();
    path_pair.second.clear();

    // generate fmindices and check for artificial sequences
    std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>> ORF_colours_tuple;
    if (is_ref)
    {
        cout << "Checking for artificial sequences..." << endl;
        ORF_colours_tuple = std::move(filter_artificial_ORFS(ORF_pair.first, input_colours, write_idx));
    } else{
        ORF_colours_tuple = std::move(sort_ORF_colours(ORF_pair.first));
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
        cout << "Calculating gene overlap" << endl;
        ORF_overlap_map = std::move(calculate_overlaps(std::get<0>(graph_tuple), ORF_pair.second, ORF_colours_tuple, overlap, 90));
    }

    // convert robin_hood maps to unordered_maps for return. Mapped by ID, so just iterate over numbers
    PyORFColoursMap ORF_colours_map;
    for (size_t i = 0; i < std::get<0>(ORF_colours_tuple).size(); i++)
    {
        ORF_colours_map[i] = std::move(std::get<0>(ORF_colours_tuple)[i]);
    }
    std::get<0>(ORF_colours_tuple).clear();

    PyORFIDMap ORF_ID_Map;
    for (size_t i = 0; i < std::get<1>(ORF_colours_tuple).size(); i++)
    {
        ORF_ID_Map[i] = std::move(std::get<1>(ORF_colours_tuple)[i]);
    }
    std::get<1>(ORF_colours_tuple).clear();

    PyunitigMap unitig_map;
    for (size_t i = 0; i < std::get<0>(graph_tuple).size(); i++)
    {
        unitig_map[i] = std::move(std::get<0>(graph_tuple)[i]);
    }
    std::get<0>(graph_tuple).clear();

    std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyunitigMap, size_t> return_tuple = std::make_tuple(ORF_overlap_map, ORF_colours_map, ORF_ID_Map, unitig_map, nb_colours);

    return return_tuple;
}

std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyunitigMap, size_t> py_ggCaller_graphbuild (const std::string& infile1,
                                                                                      const int kmer,
                                                                                      const std::vector<std::string>& start_codons,
                                                                                      const std::vector<std::string>& stop_codons_for,
                                                                                      const std::vector<std::string>& stop_codons_rev,
                                                                                      size_t num_threads,
                                                                                      bool is_ref,
                                                                                      const bool write_idx,
                                                                                      const bool repeat,
                                                                                      const bool write_graph,
                                                                                      const bool no_filter,
                                                                                      const size_t max_path_length,
                                                                                      const size_t min_ORF_length,
                                                                                      const size_t max_ORF_overlap,
                                                                                      const std::string& infile2) {
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Building coloured compacted DBG..." << endl;

    // set local variables
    const int overlap = kmer - 1;
    omp_set_num_threads(num_threads);

    if (infile2 != "NA") {
        is_ref = 0;
    }

    // generate graph, writing if write_graph == true
    size_t lastindex = infile1.find_last_of(".");
    std::string outgraph = infile1.substr(0, lastindex);
    ColoredCDBG<> ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, write_graph, outgraph);

    // get the number of colours, generate empty colour vector
    const size_t nb_colours = ccdbg.getNbColors();
    const std::vector<bool> empty_colour_arr(nb_colours, 0);

    // get colour names
    std::vector<std::string> input_colours = ccdbg.getColorNames();

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    auto graph_tuple = std::move(index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));

    // clear ccdbg to free memory
    ccdbg.clear();

    // generate complete paths
    cout << "Generating complete stop-stop paths..." << endl;
    auto path_pair = std::move(traverse_graph(graph_tuple, repeat, empty_colour_arr, max_path_length));

    // generate ORF sequences
    cout << "Generating ORF sequences from complete paths..." << endl;
    auto ORF_pair = std::move(call_ORFs(path_pair, std::get<0>(graph_tuple), stop_codons_for, start_codons, overlap, min_ORF_length));

    // clear path_pair to free memory
    path_pair.first.clear();
    path_pair.second.clear();

    // generate fmindices and check for artificial sequences
    std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>> ORF_colours_tuple;
    if (is_ref)
    {
        cout << "Checking for artificial sequences..." << endl;
        ORF_colours_tuple = std::move(filter_artificial_ORFS(ORF_pair.first, input_colours, write_idx));
    } else{
        ORF_colours_tuple = std::move(sort_ORF_colours(ORF_pair.first));
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
        cout << "Calculating gene overlap" << endl;
        ORF_overlap_map = std::move(calculate_overlaps(std::get<0>(graph_tuple), ORF_pair.second, ORF_colours_tuple, overlap, 90));
    }

    // convert robin_hood maps to unordered_maps for return. Mapped by ID, so just iterate over numbers
    PyORFColoursMap ORF_colours_map;
    for (size_t i = 0; i < std::get<0>(ORF_colours_tuple).size(); i++)
    {
        ORF_colours_map[i] = std::move(std::get<0>(ORF_colours_tuple)[i]);
    }
    std::get<0>(ORF_colours_tuple).clear();

    PyORFIDMap ORF_ID_Map;
    for (size_t i = 0; i < std::get<1>(ORF_colours_tuple).size(); i++)
    {
        ORF_ID_Map[i] = std::move(std::get<1>(ORF_colours_tuple)[i]);
    }
    std::get<1>(ORF_colours_tuple).clear();

    PyunitigMap unitig_map;
    for (size_t i = 0; i < std::get<0>(graph_tuple).size(); i++)
    {
        unitig_map[i] = std::move(std::get<0>(graph_tuple)[i]);
    }
    std::get<0>(graph_tuple).clear();

    std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyunitigMap, size_t> return_tuple = std::make_tuple(ORF_overlap_map, ORF_colours_map, ORF_ID_Map, unitig_map, nb_colours);

    return return_tuple;
}

PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    py::class_<unitigDict>(m, "UnitigMap")
            .def_readwrite("head_colour", &unitigDict::unitig_head_colour)
            .def_readwrite("tail_colour", &unitigDict::unitig_tail_colour)
            .def_readwrite("colours_equal", &unitigDict::head_tail_colours_equal)
            .def_readwrite("neighbours", &unitigDict::neighbours)
            .def_readwrite("end_contig", &unitigDict::end_contig);

    m.def("call_genes_existing", &py_ggCaller_graphexists, "Traverses pre-existing Bifrost graph, calling open reading frames.",
    py::arg("graphfile"),
    py::arg("coloursfile"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 0,
    py::arg("write_idx") = 1,
    py::arg("repeat") = 0,
    py::arg("no_repeat") = 0,
    py::arg("max_path_length") = 10000,
    py::arg("min_ORF_length") = 90,
    py::arg("max_ORF_overlap") = 60);

    m.def("call_genes_build", &py_ggCaller_graphbuild, "Builds and then traverses Bifrost graph, calling open reading frames.",
    py::arg("infile1"),
    py::arg("kmer"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 1,
    py::arg("write_idx") = 1,
    py::arg("repeat") = 0,
    py::arg("write_graph") = 1,
    py::arg("no_repeat") = 0,
    py::arg("max_path_length") = 10000,
    py::arg("min_ORF_length") = 90,
    py::arg("max_ORF_overlap") = 60,
    py::arg("infile2") = "NA");
}