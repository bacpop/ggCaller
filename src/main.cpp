//
// Created by sth19 on 15/12/2020.
//

# include "ggCaller_classes.h"

std::string generate_sequence(const unitigMap& graph_map,
                              const std::vector<int>& nodelist,
                              const std::vector<indexTriplet>& node_coords,
                              const size_t& overlap)
{
    std::string sequence;
    for (size_t i = 0; i < nodelist.size(); i++)
    {
        // initialise sequence items
        std::string unitig_seq;
        std::string substring;

        // parse information
        const auto& id = nodelist[i];
        const auto& coords = node_coords[i];
        bool strand = (id >= 0) ? true : false;

        if (strand)
        {
            unitig_seq = graph_map.at(abs(id)).unitig_seq;
        } else {
            unitig_seq = reverse_complement(graph_map.at(abs(id)).unitig_seq);
        }

        if (sequence.empty())
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - std::get<0>(coords)) + 1;
            substring = unitig_seq.substr(std::get<0>(coords), node_seq_len);
        } else
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - overlap) + 1;
            // need to account for overlap, if overlap is greater than the end of the node, sequence already accounted for
            if (node_seq_len > 0)
            {
                substring = unitig_seq.substr(overlap, node_seq_len);
            }
        }
        sequence += substring;
    }
    return sequence;
}

void write_to_file (const unitigMap& graph_map,
                    const std::string& outfile_name,
                    const ORFColoursMap& ORF_colours_map,
                    const ORFIDMap& ORF_ID_map,
                    const size_t nb_colours,
                    const size_t overlap)
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

        // generate ORF from node sequence
        const auto ORF_sequence = std::move(generate_sequence(graph_map, std::get<0>(ORF_ID_map.at(ORF.first)), std::get<1>(ORF_ID_map.at(ORF.first)), overlap));
        const auto TIS_sequence = std::move(generate_sequence(graph_map, std::get<3>(ORF_ID_map.at(ORF.first)), std::get<4>(ORF_ID_map.at(ORF.first)), overlap));

        std::string total_seq;

        total_seq += TIS_sequence;
        total_seq += ORF_sequence;

        // append to file
        outfile << ">" << std::to_string(gene_id) << "_" << colours << "\n" << total_seq << "\n";
        gene_id++;
    }

    outfile.close();
}

int main(int argc, char *argv[]) {

    int num_threads = 4;
    bool is_ref = true;
    const std::string outfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/plasmid_clique_556_list_integer_paths.fasta";
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

    const std::string graphfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/plasmid_clique_556_list.gfa";
    const std::string coloursfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/plasmid_clique_556_list.bfg_colors";

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

    std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>, ColourNodeStrandMap> ORF_tuple;
    // scope for path_pair
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

        // generate ORF sequences - get this bit to work!
        cout << "Generating ORF sequences from complete paths..." << endl;
        ORF_tuple = std::move(call_ORFs(path_pair, std::get<0>(graph_tuple), stop_codons_for, start_codons_for, overlap, min_ORF_length, is_ref, seq_idx, nb_colours));
    }

    cout << "Calculating gene overlap" << endl;
    auto overlap_map = std::move(calculate_overlaps(std::get<0>(graph_tuple), ORF_tuple, overlap, 90));

    // write fasta files to file
    write_to_file(std::get<0>(graph_tuple), outfile, std::get<0>(ORF_tuple), std::get<1>(ORF_tuple), nb_colours, overlap);

    cout << "Done." << endl;

    return 0;
}
