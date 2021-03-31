// ggCaller header
#include "ggCaller_classes.h"

GraphTuple py_index_graph_exists(const std::string& graphfile,
                               const std::string& coloursfile,
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
    GraphTuple graph_tuple = std::make_tuple(graph_pair.first, graph_pair.second, input_colours, nb_colours, overlap);

    return graph_tuple;
}

GraphTuple py_index_graph_build(const std::string& infile1,
                              const int kmer,
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
    GraphPair graph_pair;
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
    GraphTuple graph_tuple = std::make_tuple(graph_pair.first, graph_pair.second, input_colours, nb_colours, overlap);

    return graph_tuple;
}

std::pair<ORFOverlapMap, ORFVector> py_calculate_ORFs (const UnitigVector& graph_vector,
                                                     const size_t& colour_ID,
                                                     const std::unordered_set<size_t>& node_ids,
                                                     const bool& repeat,
                                                     const size_t& overlap,
                                                     const size_t& max_path_length,
                                                     bool& is_ref,
                                                     const bool& no_filter,
                                                     const std::vector<std::string>& stop_codons_for,
                                                     const std::vector<std::string>& start_codons_for,
                                                     const size_t min_ORF_length,
                                                     const size_t max_overlap,
                                                     const bool write_idx,
                                                     const std::string& FM_fasta_file)
{
    std::pair<ORFVector, NodeStrandMap> ORF_pair;
    // traverse graph, set scope for all_paths and fm_idx
    {
        // recursive traversal
        AllPaths all_paths = traverse_graph(graph_vector, colour_ID, node_ids, repeat, max_path_length);

        // if no FM_fasta_file specified, cannot generate FM Index
        if (FM_fasta_file == "NA")
        {
            is_ref = false;
        }

        // generate FM_index if is_ref
        fm_index_coll fm_idx;
        if (is_ref)
        {
            fm_idx = index_fasta(FM_fasta_file, write_idx);
        }

        // generate ORF calls
        ORF_pair = call_ORFs(all_paths, graph_vector, stop_codons_for, start_codons_for, overlap, min_ORF_length, is_ref, fm_idx);
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
        ORF_overlap_map = std::move(calculate_overlaps(graph_vector, ORF_pair, overlap, max_overlap));
    }

    std::pair<ORFOverlapMap, ORFVector> return_pair = std::make_pair(ORF_overlap_map, ORF_pair.first);

    return return_pair;
}

PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    py::class_<unitigDict>(m, "UnitigMap")
            .def_readwrite("full_colour", &unitigDict::unitig_full_colour)
            .def_readwrite("colours_equal", &unitigDict::head_tail_colours_equal)
            .def_readwrite("neighbours", &unitigDict::neighbours)
            .def_readwrite("end_contig", &unitigDict::end_contig)
            .def_readwrite("seq", &unitigDict::unitig_seq)
            .def_readwrite("forward_stop", &unitigDict::forward_stop)
            .def_readwrite("reverse_stop", &unitigDict::reverse_stop)
            .def_readwrite("full_codon", &unitigDict::full_codon)
            .def_readwrite("part_codon", &unitigDict::part_codon)
            .def_readwrite("unitig_size", &unitigDict::unitig_size);


    m.def("index_existing", &py_index_graph_exists, "Indexes pre-existing Bifrost graph.",
    py::arg("graphfile"),
    py::arg("coloursfile"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 0);

    m.def("index_build", &py_index_graph_build, "Builds and then indexes Bifrost graph.",
    py::arg("infile1"),
    py::arg("kmer"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 1,
    py::arg("write_graph") = 1,
    py::arg("infile2") = "NA");

    m.def("calculate_ORFs", &py_calculate_ORFs, "Calculates ORFs for a single colour in a graph.",
    py::arg("graph_vector"),
    py::arg("colour_ID"),
    py::arg("node_ids"),
    py::arg("repeat"),
    py::arg("overlap"),
    py::arg("max_path_length"),
    py::arg("is_ref"),
    py::arg("no_filter"),
    py::arg("stop_codons_for"),
    py::arg("start_codons_for"),
    py::arg("min_ORF_length"),
    py::arg("max_overlap"),
    py::arg("write_idx"),
    py::arg("FM_fasta_file"));
}