// ggCaller header
#include "ggCaller_classes.h"

std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyUnitigMap, size_t, size_t> py_ggCaller_graphexists (const std::string& graphfile,
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

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // initialise persistent variables
    GraphTuple graph_tuple;
    int kmer;
    int overlap;
    size_t nb_colours;
    std::vector<std::string> input_colours;

     // scope for ccdbg
    {
        // read in graph
        ColoredCDBG<> ccdbg;
        ccdbg.read(graphfile, coloursfile, num_threads);

        //set local variables
        kmer = ccdbg.getK();
        overlap = kmer - 1;

        // get the number of colours
        nb_colours = ccdbg.getNbColors();

        // get colour names
        input_colours = ccdbg.getColorNames();

        // generate codon index for graph
        cout << "Generating graph stop codon index..." << endl;
        graph_tuple = std::move(index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));
    }

    // set persistant variable
    std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>, ColourNodeStrandMap> ORF_tuple;

    // scope for path_pair and seq_idx
    {
        // generate empty colour vector
        std::vector<bool> empty_colour_arr(nb_colours, 0);

        // generate complete paths
        cout << "Generating complete stop-stop paths..." << endl;
        auto path_pair = std::move(traverse_graph(graph_tuple, repeat, empty_colour_arr, max_path_length));

        // generate FMIndexes if is_ref
        std::vector<fm_index_coll> seq_idx;
        if (is_ref)
        {
            seq_idx = generate_fmindex(input_colours, write_idx);
        }

        // generate ORF sequences
        cout << "Generating ORF sequences from complete paths..." << endl;
        ORF_tuple = std::move(call_ORFs(path_pair, std::get<0>(graph_tuple), stop_codons_for, start_codons, overlap, min_ORF_length, is_ref, seq_idx, nb_colours));
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
        cout << "Calculating gene overlap" << endl;
        ORF_overlap_map = std::move(calculate_overlaps(std::get<0>(graph_tuple), ORF_tuple, overlap, 90));
    }

    // convert robin_hood maps to unordered_maps for return. Mapped by ID, so just iterate over numbers
    PyORFColoursMap ORF_colours_map;
    for (size_t i = 0; i < std::get<0>(ORF_tuple).size(); i++)
    {
        ORF_colours_map[i] = std::move(std::get<0>(ORF_tuple)[i]);
    }

    PyORFIDMap ORF_ID_Map;
    for (size_t i = 0; i < std::get<1>(ORF_tuple).size(); i++)
    {
        ORF_ID_Map[i] = std::move(std::get<1>(ORF_tuple)[i]);
    }

    PyUnitigMap unitig_map;
    for (size_t i = 0; i < std::get<0>(graph_tuple).size(); i++)
    {
        unitig_map[i] = std::move(std::get<0>(graph_tuple)[i]);
    }

    std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyUnitigMap, size_t, size_t> return_tuple = std::make_tuple(ORF_overlap_map, ORF_colours_map, ORF_ID_Map, unitig_map, nb_colours, overlap);

    return return_tuple;
}

std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyUnitigMap, size_t, size_t> py_ggCaller_graphbuild (const std::string& infile1,
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

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // initialise persistent variables
    GraphTuple graph_tuple;
    const int overlap = kmer - 1;
    size_t nb_colours;
    std::vector<std::string> input_colours;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    // scope for ccdbg
    {
        // generate graph, writing if write_graph == true
        size_t lastindex = infile1.find_last_of(".");
        std::string outgraph = infile1.substr(0, lastindex);
        ColoredCDBG<> ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, write_graph, outgraph);

        // get the number of colours
        nb_colours = ccdbg.getNbColors();

        // get colour names
        input_colours = ccdbg.getColorNames();

        // generate codon index for graph
        cout << "Generating graph stop codon index..." << endl;
        graph_tuple = std::move(index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));
    }

    // set persistant variable
    std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>, ColourNodeStrandMap> ORF_tuple;

    // scope for path_pair and seq_idx
    {
        // generate empty colour vector
        std::vector<bool> empty_colour_arr(nb_colours, 0);

        // generate complete paths
        cout << "Generating complete stop-stop paths..." << endl;
        auto path_pair = std::move(traverse_graph(graph_tuple, repeat, empty_colour_arr, max_path_length));

        // generate FMIndexes if is_ref
        std::vector<fm_index_coll> seq_idx;
        if (is_ref)
        {
            seq_idx = generate_fmindex(input_colours, write_idx);
        }

        // generate ORF sequences
        cout << "Generating ORF sequences from complete paths..." << endl;
        ORF_tuple = std::move(call_ORFs(path_pair, std::get<0>(graph_tuple), stop_codons_for, start_codons, overlap, min_ORF_length, is_ref, seq_idx, nb_colours));
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
        cout << "Calculating gene overlap" << endl;
        ORF_overlap_map = std::move(calculate_overlaps(std::get<0>(graph_tuple), ORF_tuple, overlap, 90));
    }

    // convert robin_hood maps to unordered_maps for return. Mapped by ID, so just iterate over numbers
    PyORFColoursMap ORF_colours_map;
    for (size_t i = 0; i < std::get<0>(ORF_tuple).size(); i++)
    {
        ORF_colours_map[i] = std::move(std::get<0>(ORF_tuple)[i]);
    }

    PyORFIDMap ORF_ID_Map;
    for (size_t i = 0; i < std::get<1>(ORF_tuple).size(); i++)
    {
        ORF_ID_Map[i] = std::move(std::get<1>(ORF_tuple)[i]);
    }

    PyUnitigMap unitig_map;
    for (size_t i = 0; i < std::get<0>(graph_tuple).size(); i++)
    {
        unitig_map[i] = std::move(std::get<0>(graph_tuple)[i]);
    }

    std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyUnitigMap, size_t, size_t> return_tuple = std::make_tuple(ORF_overlap_map, ORF_colours_map, ORF_ID_Map, unitig_map, nb_colours, overlap);

    return return_tuple;
}

GraphTuple py_index_graph_exists(const std::string& graphfile,
                               const std::string& coloursfile,
                               const std::vector<std::string>& start_codons,
                               const std::vector<std::string>& stop_codons_for,
                               const std::vector<std::string>& stop_codons_rev,
                               size_t num_threads,
                               const bool is_ref)
{
        // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Reading coloured compacted DBG..." << endl;

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // initialise persistent variables
    GraphPair graph_pair;
    int kmer;
    int overlap;
    size_t nb_colours;
    std::vector<std::string> input_colours;

     // scope for ccdbg
    {
        // read in graph
        ColoredCDBG<> ccdbg;
        ccdbg.read(graphfile, coloursfile, num_threads);

        //set local variables
        kmer = ccdbg.getK();
        overlap = kmer - 1;

        // get the number of colours
        nb_colours = ccdbg.getNbColors();

        // get colour names
        input_colours = ccdbg.getColorNames();

        // generate codon index for graph
        cout << "Generating graph stop codon index..." << endl;
        graph_pair = std::move(index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(graph_pair.first, graph_pair.second, nb_colours, overlap);

    return graph_tuple;
}

GraphTuple py_index_graph_build(const std::string& infile1,
                              const int kmer,
                              const std::vector<std::string>& start_codons,
                              const std::vector<std::string>& stop_codons_for,
                              const std::vector<std::string>& stop_codons_rev,
                              size_t num_threads,
                              bool is_ref,
                              const bool write_graph,
                              const std::string& infile2)
{
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Building coloured compacted DBG..." << endl;

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // initialise persistent variables
    GraphPair graph_tuple;
    int overlap = kmer - 1;
    size_t nb_colours;
    std::vector<std::string> input_colours;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    // scope for ccdbg
    {
        // generate graph, writing if write_graph == true
        size_t lastindex = infile1.find_last_of(".");
        std::string outgraph = infile1.substr(0, lastindex);
        ColoredCDBG<> ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, write_graph, outgraph);

        // get the number of colours
        nb_colours = ccdbg.getNbColors();

        // get colour names
        input_colours = ccdbg.getColorNames();

        // generate codon index for graph
        cout << "Generating graph stop codon index..." << endl;
        graph_pair = std::move(index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(graph_pair.first, graph_pair.second, nb_colours, overlap);

    return graph_tuple;
}

std::pair<ORFOverlapMap, PyORFIDMap> calculate_ORFs ()


PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    py::class_<unitigDict>(m, "UnitigMap")
            .def_readwrite("head_colour", &unitigDict::unitig_head_colour)
            .def_readwrite("tail_colour", &unitigDict::unitig_tail_colour)
            .def_readwrite("colours_equal", &unitigDict::head_tail_colours_equal)
            .def_readwrite("neighbours", &unitigDict::neighbours)
            .def_readwrite("end_contig", &unitigDict::end_contig)
            .def_readwrite("seq", &unitigDict::unitig_seq)
            .def_readwrite("forward_stop" &unitigDict::forward_stop)
            .def_readwrite("reverse_stop" &unitigDict::reverse_stop)
            .def_readwrite("full_codon" &unitigDict::full_codon)
            .def_readwrite("part_codon" &unitigDict::part_codon)
            .def_readwrite("unitig_size" &unitigDict::unitig_size);

//    m.def("call_genes_existing", &py_ggCaller_graphexists, "Traverses pre-existing Bifrost graph, calling open reading frames.",
//    py::arg("graphfile"),
//    py::arg("coloursfile"),
//    py::arg("start_codons"),
//    py::arg("stop_codons_for"),
//    py::arg("stop_codons_rev"),
//    py::arg("num_threads") = 1,
//    py::arg("is_ref") = 0,
//    py::arg("write_idx") = 1,
//    py::arg("repeat") = 0,
//    py::arg("no_repeat") = 0,
//    py::arg("max_path_length") = 10000,
//    py::arg("min_ORF_length") = 90,
//    py::arg("max_ORF_overlap") = 60);
//
//    m.def("call_genes_build", &py_ggCaller_graphbuild, "Builds and then traverses Bifrost graph, calling open reading frames.",
//    py::arg("infile1"),
//    py::arg("kmer"),
//    py::arg("start_codons"),
//    py::arg("stop_codons_for"),
//    py::arg("stop_codons_rev"),
//    py::arg("num_threads") = 1,
//    py::arg("is_ref") = 1,
//    py::arg("write_idx") = 1,
//    py::arg("repeat") = 0,
//    py::arg("write_graph") = 1,
//    py::arg("no_repeat") = 0,
//    py::arg("max_path_length") = 10000,
//    py::arg("min_ORF_length") = 90,
//    py::arg("max_ORF_overlap") = 60,
//    py::arg("infile2") = "NA");

    m.def("index_existing", &py_index_graph_exists, "Traverses pre-existing Bifrost graph, calling open reading frames.",
    py::arg("graphfile"),
    py::arg("coloursfile"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 0);

    m.def("index_build", &py_ggCaller_graphbuild, "Builds and then traverses Bifrost graph, calling open reading frames.",
    py::arg("infile1"),
    py::arg("kmer"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 1,
    py::arg("write_graph") = 1,
    py::arg("infile2") = "NA");
}