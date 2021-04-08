//
// Created by sth19 on 15/12/2020.
//

# include "ggCaller_classes.h"

std::string generate_sequence(const UnitigVector& graph_vector,
                              const std::vector<int>& nodelist,
                              const std::vector<indexPair>& node_coords,
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
            unitig_seq = graph_vector.at(abs(id) - 1).unitig_seq;
        } else {
            unitig_seq = reverse_complement(graph_vector.at(abs(id) - 1).unitig_seq);
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

void write_to_file (const robin_hood::unordered_map<std::string, std::vector<bool>>& ORF_print_map,
                    const std::string& outfile_name)
{
    int gene_id = 1;
    ofstream outfile;
    outfile.open(outfile_name);

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
        outfile << ">" << std::to_string(gene_id) << "_" << colours << "\n" << ORF.first << "\n";
        gene_id++;
    }

    outfile.close();
}

int main(int argc, char *argv[]) {

    int num_threads = 4;
    bool is_ref = true;
    const std::string outfile = "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/group3_capsular_fa_list.fasta";
    omp_set_num_threads(num_threads);
    const bool write_graph = true;
    const bool write_idx = true;
    const bool repeat = false;
    const size_t max_path_length = 10000;
    const size_t min_ORF_length = 90;
    const size_t max_ORF_overlap = 60;
    const bool no_filter = false;

    const std::vector<std::string> stop_codons_for = {"TAA", "TGA", "TAG"};
    const std::vector<std::string> stop_codons_rev = {"TTA", "TCA", "CTA"};

    const std::vector<std::string> start_codons_for = {"ATG", "GTG", "TTG"};


    // Set number of threads
    if (num_threads < 1) {
        num_threads = 1;
    }

//    GraphTuple graph_tuple = py_index_graph_exists(
//        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.gfa",
//        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.bfg_colors",
//        stop_codons_for, stop_codons_rev, num_threads, is_ref);

    GraphTuple graph_tuple = py_index_graph_build(
            "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group3_capsular_fa_list.txt", 31,
            stop_codons_for, stop_codons_rev, num_threads, is_ref, write_graph, "NA");

    const auto& graph_vector = std::get<0>(graph_tuple);
    const auto& node_colour_vector = std::get<1>(graph_tuple);
    const auto& input_colours = std::get<2>(graph_tuple);
    const auto& nb_colours = std::get<3>(graph_tuple);
    const auto& overlap = std::get<4>(graph_tuple);

    // initialise print map
    robin_hood::unordered_map<std::string, std::vector<bool>> ORF_print_map;

    for (int colour_ID = 0; colour_ID < node_colour_vector.size(); colour_ID++)
    {
        const auto& node_set = node_colour_vector.at(colour_ID);

        std::pair<ORFOverlapMap, ORFVector> ORF_pair = py_calculate_ORFs(graph_vector, colour_ID, node_set, repeat,
                                                               overlap, max_path_length, is_ref, no_filter,
                                                               stop_codons_for, start_codons_for, min_ORF_length,
                                                               max_ORF_overlap, write_idx, input_colours.at(colour_ID));

        auto& ORF_vector = ORF_pair.second;

        // add to ORF_print_map for writing to file
        std::vector<bool> empty_colours_vector(nb_colours, 0);
        for (const auto ORF : ORF_vector)
        {
            // generate ORF seq
            const auto ORF_sequence = std::move(generate_sequence(graph_vector, std::get<0>(ORF), std::get<1>(ORF), overlap));
            const auto TIS_sequence = std::move(generate_sequence(graph_vector, std::get<3>(ORF), std::get<4>(ORF), overlap));

            std::string total_seq;

            total_seq += TIS_sequence;
            total_seq += ORF_sequence;

            // add colour index to the colours vector
            if (ORF_print_map.find(total_seq) == ORF_print_map.end())
            {
                ORF_print_map[total_seq] = empty_colours_vector;
            }

            ORF_print_map[total_seq][colour_ID] = 1;
        }
    }

    // print to file
    write_to_file(ORF_print_map, outfile);

    return 0;
}