// ggCaller header
#include "ggCaller_classes.h"

int py_ggCaller_graphexists (const std::string& graphfile,
                             const std::string& coloursfile,
                             const std::string& outfile,
                             const std::vector<std::string>& start_codons,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             size_t num_threads,
                             const bool is_ref,
                             const bool write_idx,
                             const bool repeat,
                             const size_t& max_path_length,
                             const size_t& min_ORF_length) {
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
    const auto graph_tuple = index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours);

    // generate complete paths
    cout << "Generating complete stop-stop paths..." << endl;
    const auto path_tuple = traverse_graph(ccdbg, graph_tuple, repeat, empty_colour_arr, max_path_length);

    // generate ORF sequences
    cout << "Generating ORF sequences from complete paths..." << endl;
    auto ORF_tuple = call_ORFs(ccdbg, path_tuple, stop_codons_for, start_codons, overlap, min_ORF_length);

    // generate fmindices and check for artificial sequences if sequences supplied are reference files
    if (is_ref)
    {
        cout << "Checking for artificial sequences..." << endl;
        auto ORF_colours_map = filter_artificial_ORFS(std::get<0>(ORF_tuple), std::get<1>(ORF_tuple), input_colours, write_idx);
    }

    // write fasta files to file
    cout << "Writing gene calls to file..." << endl;
    write_to_file(outfile, std::get<0>(ORF_tuple));



    cout << "Done." << endl;

    return 0;
}

int py_ggCaller_graphbuild (const std::string& infile1,
                            const int& kmer,
                            const std::string& outfile,
                            const std::vector<std::string>& start_codons,
                            const std::vector<std::string>& stop_codons_for,
                            const std::vector<std::string>& stop_codons_rev,
                            size_t num_threads,
                            bool is_ref,
                            const bool write_idx,
                            const bool repeat,
                            const bool write_graph,
                            const size_t& max_path_length,
                            const size_t& min_ORF_length,
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

    ColoredCDBG<> ccdbg;

    if (write_graph) {
        // build and write graph
        size_t lastindex = infile1.find_last_of(".");
        std::string outgraph = infile1.substr(0, lastindex);
        ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, true, outgraph);
    } else {
        // build graph only
        std::string outgraph = "NA";
        ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, false, outgraph);
    }

    // iterate over unitigs, search for presence of sub-sequence
    cout << "Constructing unitig matrices for graph traversal..." << endl;

    // get the number of colours, generate empty colour vector
    const size_t nb_colours = ccdbg.getNbColors();
    const std::vector<bool> empty_colour_arr(nb_colours, 0);

    // get colour names
    std::vector<std::string> input_colours = ccdbg.getColorNames();

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    const auto graph_tuple = index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours);

    // generate complete paths
    cout << "Generating complete stop-stop paths..." << endl;
    const auto path_tuple = traverse_graph(ccdbg, graph_tuple, repeat, empty_colour_arr, max_path_length);

    // generate ORF sequences
    cout << "Generating ORF sequences from complete paths..." << endl;
    auto ORF_tuple = call_ORFs(ccdbg, path_tuple, stop_codons_for, start_codons, overlap, min_ORF_length);

    // generate fmindices and check for artificial sequences if sequences supplied are reference files
    if (is_ref)
    {
        cout << "Checking for artificial sequences..." << endl;
        auto ORF_colours_map = filter_artificial_ORFS(std::get<0>(ORF_tuple), std::get<1>(ORF_tuple), input_colours, write_idx);
    }

    // write fasta files to file
    cout << "Writing gene calls to file..." << endl;
    write_to_file(outfile, std::get<0>(ORF_tuple));

    cout << "Done." << endl;

    return 0;
}

PYBIND11_MODULE(ggCaller_cpp, m)
{
    m.doc() = "Call ORFs in Bifrost graph.";

    m.def("call_genes_existing", &py_ggCaller_graphexists, "Traverses pre-existing Bifrost graph, calling open reading frames.",
    py::arg("graphfile"),
    py::arg("coloursfile"),
    py::arg("outfile"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 0,
    py::arg("write_idx") = 1,
    py::arg("repeat") = 0,
    py::arg("max_path_length") = 10000,
    py::arg("min_ORF_length") = 90);

    m.def("call_genes_build", &py_ggCaller_graphbuild, "Builds and then traverses Bifrost graph, calling open reading frames.",
    py::arg("infile1"),
    py::arg("kmer"),
    py::arg("outfile"),
    py::arg("start_codons"),
    py::arg("stop_codons_for"),
    py::arg("stop_codons_rev"),
    py::arg("num_threads") = 1,
    py::arg("is_ref") = 1,
    py::arg("write_idx") = 1,
    py::arg("repeat") = 0,
    py::arg("write_graph") = 1,
    py::arg("max_path_length") = 10000,
    py::arg("min_ORF_length") = 90,
    py::arg("infile2") = "NA");
}