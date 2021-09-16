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
std::vector<std::string> parse_fasta (const std::string& fasta)
{
    std::vector<std::string> seq_list;
    std::string line, DNA_sequence;
    std::ifstream infile(fasta);

    // ignore header
    std::getline(infile, line);

    // parse fasta
    while (std::getline(infile, line)) {
        // remove new line characters
        line.erase(std::remove(line.begin(), line.end(), '\n'),
                   line.end());

        // erase DNA_sequence
        DNA_sequence.clear();

        // line may be empty so you *must* ignore blank lines
        if (line.empty())
        {
            continue;
        }

        if (line[0] != '>') {
            DNA_sequence = line;
        }

        // add to seq list
        if (DNA_sequence.size() != 0)
        {
            seq_list.push_back(DNA_sequence);
        }
    }

    return seq_list;
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
            "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.gfa",
            "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.bfg_colors",
            stop_codons_for, stop_codons_rev, num_threads, is_ref);

    const auto& node_colour_vector = std::get<0>(graph_tuple);
    const auto& input_colours = std::get<1>(graph_tuple);
    const auto& nb_colours = std::get<2>(graph_tuple);
    const auto& overlap = std::get<3>(graph_tuple);

    // initialise print map
    robin_hood::unordered_map<std::string, std::vector<bool>> ORF_print_map;

    // initialise colour_ORF_map
    std::unordered_map<size_t, ORFNodeMap> colour_ORF_map;


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
//        std::vector<std::pair<size_t,size_t>> random_ORFs;
//        for (int i = 0; i < round(ORF_vector.size()/10); i++)
//        {
//            size_t a = rand() % ORF_vector.size() + 1;
//            size_t b = rand() % ORF_vector.size() + 1;
//            random_ORFs.push_back(std::make_pair(a, b));
//        }

        // parse ORFs known to be genes
        auto known_genes = parse_fasta("/mnt/c/Users/sth19/CLionProjects/Bifrost_API/data/plasmid_clique_119_230_372_test_ORFs_for_panaroo_new.fasta");
        const std::vector<size_t> fasta_order = {1176, 327, 1219, 826, 233, 990, 222, 606, 871, 1127, 713, 563, 1146, 484, 490, 493, 752, 665, 496, 981, 852, 618, 345, 598, 836, 792, 115, 476, 471, 474, 1113, 136, 520, 99, 932, 697, 923, 1215, 744, 971, 846, 938, 1057, 1118, 571, 867, 687, 953, 975, 645, 117, 1061, 53, 65, 605, 1145, 556, 373, 652, 1088, 384, 276, 30, 928, 1233, 992, 1175, 1209, 219, 1207, 439, 1104, 551, 104, 899, 946, 1016, 180, 845, 302, 1208, 607, 770, 107, 336, 632, 121, 514, 38, 466, 39, 91, 375, 221, 1144, 339, 492, 644, 74, 468, 1131, 379, 828, 743, 802, 945, 1070, 1106, 661, 21, 358, 879};
        const std::vector<size_t> original_targets = {514, 520, 21, 30, 1057, 1061, 38, 551, 39, 556, 1070, 563, 53, 571, 1088, 65, 1106, 91, 605, 606, 1118, 607, 99, 1127, 104, 618, 107, 115, 117, 1144, 121, 1145, 1146, 645, 136, 652, 1175, 1176, 665, 687, 180, 1207, 1208, 697, 1209, 1215, 1233, 219, 221, 222, 744, 752, 770, 276, 792, 802, 302, 836, 845, 846, 336, 339, 852, 345, 358, 871, 879, 373, 375, 384, 899, 923, 928, 932, 938, 945, 946, 953, 975, 466, 981, 471, 474, 476, 992, 484, 490, 492, 493, 496, 1016};


        // sequence vector for testing
        std::vector<size_t> known_gene_ids(fasta_order.size());
        size_t ORF_count = 0;
        for (const auto ORF : ORF_vector)
        {
            const auto ORF_sequence = std::move(unitig_graph.generate_sequence(std::get<0>(ORF), std::get<1>(ORF), overlap));
            auto index_iter = std::find(known_genes.begin(), known_genes.end(), ORF_sequence);
            if (index_iter != known_genes.end())
            {
                auto index = std::distance(known_genes.begin(), index_iter);
                known_gene_ids[index] = ORF_count;
            }
            ORF_count++;
        }

        // generate a map to hold path IDs and associated ORFs
        std::vector<std::pair<size_t, size_t>> original_to_new_IDs;
        std::vector<size_t> target_ORFs(original_targets.size());
        for (int i = 0; i < known_gene_ids.size(); i++)
        {
            const auto& ORF_info = ORF_vector.at(known_gene_ids.at(i));
            auto it =  std::find(original_targets.begin(), original_targets.end(), fasta_order.at(i));

            if (it != original_targets.end())
            {
                size_t index = it - original_targets.begin();
                target_ORFs[index] = known_gene_ids.at(i);
            }

            std::pair<size_t, size_t> new_pair;
            new_pair.first = fasta_order.at(i);
            new_pair.second = known_gene_ids.at(i);

            original_to_new_IDs.push_back(new_pair);
        }

        // generate a map that holds only high scoring ORFs
        ORFNodeMap ORF_node_map;
        for (size_t i = 0; i < ORF_vector.size(); i++)
        {
            ORF_node_map[i] = ORF_vector.at(i);
        }

        // move to colour_ORF_map
        colour_ORF_map[colour_ID] = std::move(ORF_node_map);

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
    unitig_graph.generate_clusters(colour_ORF_map, overlap, 0.7, 0.98);


    // print to file
    //write_to_file(ORF_print_map, outfile);

    return 0;
}