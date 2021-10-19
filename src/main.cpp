//
// Created by sth19 on 15/12/2020.
//

#include "unitigDict.h"
#include "traversal.h"
#include "call_ORFs.h"
#include "match_string.h"
#include "gene_overlap.h"
#include "graph.h"
#include "ORF_clustering.h"

// parse a fasta and return sequences
std::pair<std::vector<std::pair<int, int>>, std::vector<std::string>> parse_fasta (const std::string& fasta)
{

    std::vector<std::pair<int, int>> id_list;
    std::vector<std::string> seq_list;
    std::string line;
    std::ifstream infile(fasta);

    // create entry pair
    std::pair<int, int> entry_pair;
    std::string entry_seq;

    std::string delimiter = "_";

    // parse fasta
    while (std::getline(infile, line)) {
        // remove new line characters
        line.erase(std::remove(line.begin(), line.end(), '\n'),
                   line.end());

        // line may be empty so you *must* ignore blank lines
        if (line.empty())
        {
            continue;
        }

        if (line[0] == '>') {
            int genome = std::stoi(line.substr(1, line.find(delimiter)));
            int ORF_ID = std::stoi(line.substr(line.find(delimiter) + 1));
            std::get<0>(entry_pair) = genome;
            std::get<1>(entry_pair) = ORF_ID;
        } else
        {
            entry_seq = line;
            id_list.push_back(std::move(entry_pair));
            seq_list.push_back(std::move(entry_seq));
        }
    }

    auto return_pair = std::make_pair(id_list, seq_list);
    return return_pair;
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

//    GraphTuple graph_tuple = unitig_graph.build(
//            "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.txt", 31,
//            stop_codons_for, stop_codons_rev, num_threads, is_ref, write_graph, "NA");

    GraphTuple graph_tuple = unitig_graph.read(
            "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/plasmid_clique_119_230_372_list.gfa",
            "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/plasmid_clique_119_230_372_list.bfg_colors",
            stop_codons_for, stop_codons_rev, num_threads, is_ref);

    const auto& node_colour_vector = std::get<0>(graph_tuple);
    const auto& input_colours = std::get<1>(graph_tuple);
    const auto& nb_colours = std::get<2>(graph_tuple);
    const auto& overlap = std::get<3>(graph_tuple);

    // initialise print map
    robin_hood::unordered_map<std::string, std::vector<bool>> ORF_print_map;

    // initialise colour_ORF_map
    std::unordered_map<size_t, ORFNodeMap> colour_ORF_map;

    // parse ORFs known to be genes
    const auto known_genes_pair = parse_fasta("/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/group3_test_ORFs_for_panaroo.fasta");
    const auto& known_genes_entries = known_genes_pair.first;
    const auto& known_genes = known_genes_pair.second;


    //#pragma omp parallel for
//    for (size_t colour_ID = 0; colour_ID < node_colour_vector.size(); colour_ID++)
    for (size_t colour_ID = 10; colour_ID < 11; colour_ID++)
    {

        std::pair<ORFOverlapMap, ORFVector> ORF_pair = unitig_graph.findORFs(colour_ID, node_colour_vector[colour_ID], repeat,
                                                                             overlap, max_path_length, is_ref, no_filter,
                                                                             stop_codons_for, start_codons_for, min_ORF_length,
                                                                             max_ORF_overlap, write_idx, input_colours[colour_ID]);

        auto& ORF_vector = ORF_pair.second;

        // test out ORF pairing
        // generate random numbers of ORFs to add (10% of size of total ORFs)
//        std::vector<std::pair<size_t,size_t>> random_ORFs;
//        for (int i = 0; i < round(ORF_vector.size()/10); i++)
//        {
//            size_t a = rand() % ORF_vector.size() + 1;
//            size_t b = rand() % ORF_vector.size() + 1;
//            random_ORFs.push_back(std::make_pair(a, b));
//        }

        // sequence vector for testing
        std::unordered_map<size_t, size_t> known_gene_ids;
        size_t ORF_count = 0;
        std::vector<std::string> ORF_sequence_list;
        for (const auto ORF : ORF_vector)
        {
            const auto ORF_sequence = std::move(unitig_graph.generate_sequence(std::get<0>(ORF), std::get<1>(ORF), overlap));
            ORF_sequence_list.push_back(ORF_sequence);
            auto index_iter = std::find(known_genes.begin(), known_genes.end(), ORF_sequence);
            if (index_iter != known_genes.end())
            {
                bool colour_found = false;
                while (!colour_found && index_iter != known_genes.end())
                {
                    auto index = std::distance(known_genes.begin(), index_iter);
                    // check if gene is called for colour
                    if (known_genes_entries.at(index).first == colour_ID)
                    {
                        known_gene_ids[index] = ORF_count;
                        colour_found = true;
                    }
                    index_iter = std::find(index_iter + 1, known_genes.end(), ORF_sequence);
                }
            }
            ORF_count++;
        }

        // generate a map that holds only high scoring ORFs
        ORFNodeMap ORF_node_map;
        for (const auto& gene : known_gene_ids)
        {
            // get true entry number
            const auto& ORF_ID = known_genes_entries.at(gene.first).second;

            ORF_node_map[ORF_ID] = ORF_vector.at(gene.second);
        }

        // move to colour_ORF_map
        colour_ORF_map[colour_ID] = std::move(ORF_node_map);


//        // generate a map to hold path IDs and associated ORFs
//        std::vector<std::pair<size_t, size_t>> original_to_new_IDs;
//        std::vector<size_t> target_ORFs(original_targets.size());
//        for (int i = 0; i < known_gene_ids.size(); i++)
//        {
//            const auto& ORF_info = ORF_vector.at(known_gene_ids.at(i));
//            auto it =  std::find(original_targets.begin(), original_targets.end(), fasta_order.at(i));
//
//            if (it != original_targets.end())
//            {
//                size_t index = it - original_targets.begin();
//                target_ORFs[index] = known_gene_ids.at(i);
//            }
//
//            std::pair<size_t, size_t> new_pair;
//            new_pair.first = fasta_order.at(i);
//            new_pair.second = known_gene_ids.at(i);
//
//            original_to_new_IDs.push_back(new_pair);
//        }


        //auto neighbours = unitig_graph.connect_ORFs(colour_ID, ORF_vector, target_ORFs, 10000);


//        // add to ORF_print_map for writing to file
//        std::vector<bool> empty_colours_vector(nb_colours, 0);
//        for (const auto ORF : ORF_vector)
//        {
//            // generate ORF seq
//            const auto ORF_sequence = std::move(unitig_graph.generate_sequence(std::get<0>(ORF), std::get<1>(ORF), overlap));
//            const auto TIS_sequence = std::move(unitig_graph.generate_sequence(std::get<3>(ORF), std::get<4>(ORF), overlap));
//
//            std::string total_seq;
//
//            total_seq += TIS_sequence;
//            total_seq += ORF_sequence;
//
//            // add colour index to the colours vector
//            if (ORF_print_map.find(total_seq) == ORF_print_map.end())
//            {
//                ORF_print_map[total_seq] = empty_colours_vector;
//            }
//
//            ORF_print_map[total_seq][colour_ID] = 1;
//        }
    }

    // calculate identity between ORFs in colour_ORF_map
    const auto return_pair = unitig_graph.generate_clusters(colour_ORF_map, overlap, 0.7, 0.98);


    // print to file
    //write_to_file(ORF_print_map, outfile);

    return 0;
}