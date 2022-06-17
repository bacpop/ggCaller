#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "call_ORFs.h"
#include "unitigDict.h"
#include "indexing.h"
#include "gene_overlap.h"
#include "ORF_connection.h"

PathVector iter_nodes_binary (const ColoredCDBG<MyUnitigMap>& ccdbg,
                              const std::vector<Kmer>& head_kmer_arr,
                              const NodeTuple& head_node_tuple,
                              const size_t current_colour,
                              const size_t length_max,
                              const size_t overlap,
                              const bool repeat,
                              const bool is_ref,
                              const boost::dynamic_bitset<>& ref_set,
                              const fm_index_coll& fm_idx);

PathVector iter_nodes_path (const ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::vector<Kmer>& head_kmer_arr,
                            const NodeTuple& head_node_tuple,
                            const size_t current_colour,
                            const size_t length_max,
                            const size_t overlap,
                            const fm_index_coll& fm_idx,
                            const robin_hood::unordered_map<int, robin_hood::unordered_set<int>>& edge_map,
                            const bool ref_strand,
                            const std::pair<int, int>& start_end_nodes);

ORFNodeRobMap calculate_genome_paths(const std::vector<Kmer>& head_kmer_arr,
                                     ColoredCDBG<MyUnitigMap>& ccdbg,
                                     const std::string& fasta_file,
                                     const int kmer,
                                     const int colour_ID,
                                     const size_t nb_colours,
                                     const size_t& max_path_length,
                                     const std::vector<std::string>& stop_codons_for,
                                     const std::vector<std::string>& start_codons_for,
                                     const size_t min_ORF_length,
                                     torch::jit::script::Module& TIS_model,
                                     const double& minimum_ORF_score,
                                     const bool no_filter,
                                     tbb::concurrent_unordered_map<size_t, float>& all_TIS_scores);

ORFNodeRobMap traverse_graph(const ColoredCDBG<MyUnitigMap>& ccdbg,
                             const std::vector<Kmer>& head_kmer_arr,
                             const size_t colour_ID,
                             const std::vector<size_t>& node_ids,
                             const bool repeat,
                             const size_t max_path_length,
                             const size_t overlap,
                             const bool is_ref,
                             const boost::dynamic_bitset<>& ref_set,
                             const fm_index_coll& fm_idx,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& start_codons_for,
                             const size_t min_ORF_length,
                             torch::jit::script::Module& TIS_model,
                             const double& minimum_ORF_score,
                             const bool no_filter,
                             tbb::concurrent_unordered_map<size_t, float>& all_TIS_scores);

#endif //TRAVERSAL_H
