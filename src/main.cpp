//
// Created by sth19 on 15/12/2020.
//

# include "ggCaller_classes.h"

int main(int argc, char *argv[]) {

    int num_threads = 1;
    bool is_ref = true;
    const std::string outfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/calls.fasta";
    omp_set_num_threads(num_threads);
    const bool write_graph = true;
    const bool write_idx = true;
    const bool repeat = false;
    const int max_path_length = 10000;
    const int min_ORF_length = 90;

    const std::vector<std::string> stop_codons_for = {"TAA", "TGA", "TAG"};
    const std::vector<std::string> stop_codons_rev = {"TTA", "TCA", "CTA"};

    const std::vector<std::string> start_codons_for = {"ATG", "GTG", "TTG"};


    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact coloured DBG
//    cout << "Building coloured compacted DBG..." << endl;
//
//    const std::string infile1 = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group2_capsular_fa_list.txt";
//    const std::string infile2 = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group3_capsular_fa_list.txt";
//
//    const int kmer = 31;
//
//
//    if (infile2 != "NA") {
//        is_ref = 0;
//    }
//
//    ColoredCDBG<> ccdbg;
//
//
//    if (write_graph) {
//        // build and write graph
//        size_t lastindex = infile1.find_last_of(".");
//        std::string outgraph = infile1.substr(0, lastindex);
//        ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, true, outgraph);
//    } else {
//        // build graph only
//        std::string outgraph = "NA";
//        ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, false, outgraph);
//    }


    cout << "Reading coloured compacted DBG..." << endl;

    const std::string graphfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group3_capsular_fa_list.gfa";
    const std::string coloursfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group3_capsular_fa_list.bfg_colors";

    // read in graph
    ColoredCDBG<> ccdbg;
    ccdbg.read(graphfile, coloursfile, num_threads);

    // set local variables
    const int kmer = ccdbg.getK();
    const int overlap = kmer - 1;

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

    // generate ORF sequences - get this bit to work!
    cout << "Generating ORF sequences from complete paths..." << endl;
    auto ORF_tuple = call_ORFs(ccdbg, path_tuple, stop_codons_for, start_codons_for, overlap, min_ORF_length);

    int ORF_num = std::get<0>(ORF_tuple).at("+").size();
    int path_num = std::get<1>(ORF_tuple).at("+").size();

    // generate fmindices and check for artificial sequences
    cout << "Checking for artificial sequences..." << endl;
    auto ORF_colours_tuple = filter_artificial_ORFS(std::get<0>(ORF_tuple), std::get<1>(ORF_tuple), input_colours, write_idx);

    ORF_num = std::get<0>(ORF_tuple).at("+").size();
    path_num = std::get<1>(ORF_tuple).at("+").size();

    cout << "Calculating gene overlap" << endl;
    auto ORF_overlap_map = calculate_overlaps(std::get<0>(graph_tuple), std::get<1>(ORF_tuple), ORF_colours_tuple, overlap, 90);

    // write fasta files to file
    cout << "Writing gene calls to file..." << endl;
    write_to_file(outfile, std::get<0>(ORF_tuple));

    cout << "Done." << endl;

    return 0;
}
