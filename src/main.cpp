//
// Created by sth19 on 15/12/2020.
//

# include "ggCaller_classes.h"

void write_to_file (const std::string& outfile_name,
                    const ORFColoursMap& ORF_colours_map,
                    const ORFIDMap& ORF_ID_map,
                    const size_t nb_colours)
{
    int gene_id = 1;
    ofstream outfile;
    outfile.open(outfile_name);

    robin_hood::unordered_map<size_t, std::vector<bool>> ORF_print_map;
    std::vector<bool> empty_vector(nb_colours, 0);
    for (size_t i = 0; i < ORF_ID_map.size(); i++)
    {
        ORF_print_map[i] = empty_vector;
    }

    // iterate over all colours in colours_map
    for (const auto& colour : ORF_colours_map)
    {
        // iterate over all ORF_IDs that have that colour and add in correct positive into vector
        for (const auto& ORF : colour.second)
        {
            ORF_print_map[ORF][colour.first] = 1;
        }
    }

    // iterate over ORFs in all_ORFs
    for (const auto& ORF : ORF_print_map)
    {
        // generate string for colours
        std::string colours;
        for (const auto& i : ORF.second)
        {
            colours += std::to_string(i);
        }
        // append to file
        outfile << ">" << std::to_string(gene_id) << "_" << colours << "\n" << ORF_ID_map.at(ORF.first).first << "\n";
        gene_id++;
    }

    outfile.close();
}

int main(int argc, char *argv[]) {

    int num_threads = 4;
    bool is_ref = true;
    const std::string outfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/group3_capsular_fa_list_profiling_4threads.fasta";
    //omp_set_num_threads(num_threads);
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
//    const std::string infile1 = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/plasmid_clique_556_list.txt";
//    const std::string infile2 = "NA";
//
//    const int kmer = 31;
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

    //set local variables
    const int kmer = ccdbg.getK();
    const int overlap = kmer - 1;

    // get the number of colours, generate empty colour vector
    const size_t nb_colours = ccdbg.getNbColors();
    const std::vector<bool> empty_colour_arr(nb_colours, 0);

    // get colour names
    std::vector<std::string> input_colours = ccdbg.getColorNames();

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    const auto graph_tuple = std::move(index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));

    // clear ccdbg to free memory
    ccdbg.clear();

    // generate complete paths
    cout << "Generating complete stop-stop paths..." << endl;
    auto path_pair = std::move(traverse_graph(graph_tuple, repeat, empty_colour_arr, max_path_length));

    // generate ORF sequences - get this bit to work!
    cout << "Generating ORF sequences from complete paths..." << endl;
    auto ORF_pair = std::move(call_ORFs(path_pair, std::get<0>(graph_tuple), stop_codons_for, start_codons_for, overlap, min_ORF_length));

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

    cout << "Calculating gene overlap" << endl;
    auto overlap_map = std::move(calculate_overlaps(std::get<0>(graph_tuple), ORF_pair.second, ORF_colours_tuple, overlap, 90));

    // write fasta files to file
    write_to_file (outfile, std::get<0>(ORF_colours_tuple), std::get<1>(ORF_colours_tuple), nb_colours);

    cout << "Done." << endl;

    return 0;
}
