//
// Created by sth19 on 15/12/2020.
//

#include "unitigDict.h"
#include "traversal.h"
#include "call_ORFs.h"
#include "match_string.h"
#include "gene_overlap.h"
#include "graph.h"

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
    const std::string outfile = "/home/shorsfield/software/ggCaller/all_test_capsular_loci_list_v1.2.4.fasta";
    //omp_set_num_threads(num_threads);
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

    // initialise and build graph
    Graph unitig_graph = Graph();

    GraphTuple graph_tuple = unitig_graph.build(
            "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.txt", 31,
            stop_codons_for, stop_codons_rev, num_threads, is_ref, write_graph, "NA");

//    GraphTuple graph_tuple = unitig_graph.read(
//            "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group3_capsular_fa_list.gfa",
//            "/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group3_capsular_fa_list.bfg_colors",
//            stop_codons_for, stop_codons_rev, num_threads, is_ref);

    const auto& node_colour_vector = std::get<0>(graph_tuple);
    const auto& input_colours = std::get<1>(graph_tuple);
    const auto& nb_colours = std::get<2>(graph_tuple);
    const auto& overlap = std::get<3>(graph_tuple);

    // initialise print map
    robin_hood::unordered_map<std::string, std::vector<bool>> ORF_print_map;

    //#pragma omp parallel for
    for (size_t colour_ID = 0; colour_ID < node_colour_vector.size(); colour_ID++)
    {

        std::pair<ORFOverlapMap, ORFVector> ORF_pair = unitig_graph.findORFs(colour_ID, node_colour_vector[colour_ID], repeat,
                                                                             overlap, max_path_length, is_ref, no_filter,
                                                                             stop_codons_for, start_codons_for, min_ORF_length,
                                                                             max_ORF_overlap, write_idx, input_colours[colour_ID]);

        auto& ORF_vector = ORF_pair.second;

        // test out ORF pairing
        // generate random numbers of ORFs to add (10% of size of total ORFs)
        std::vector<std::pair<size_t,size_t>> random_ORFs;
        for (int i = 0; i < round(ORF_vector.size()/10); i++)
        {
            size_t a = rand() % ORF_vector.size() + 1;
            size_t b = rand() % ORF_vector.size() + 1;
            random_ORFs.push_back(std::make_pair(a, b));
        }

        // iterate over random numbers, inserting and then finding nearest ORFs
        unitig_graph.add_ORF_info(colour_ID, random_ORFs, ORF_vector);

        auto neighbours = unitig_graph.get_neighbouring_ORFs(colour_ID, random_ORFs, ORF_vector);

        // add to ORF_print_map for writing to file
        std::vector<bool> empty_colours_vector(nb_colours, 0);
        for (const auto ORF : ORF_vector)
        {
            // generate ORF seq
            const auto ORF_sequence = std::move(unitig_graph.generate_sequence(std::get<0>(ORF), std::get<1>(ORF), overlap));
            const auto TIS_sequence = std::move(unitig_graph.generate_sequence(std::get<3>(ORF), std::get<4>(ORF), overlap));

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