#ifndef GENE_REFINDING_H
#define GENE_REFINDING_H

#include "unitigDict.h"
#include "match_string.h"
#include "ORF_clustering.h"

std::pair<std::vector<std::vector<size_t>>, std::string> calculate_node_ranges(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                               const std::vector<Kmer>& head_kmer_arr,
                                                                               const int& overlap,
                                                                               const std::vector<int>& full_nodelist);

std::pair<std::vector<int>, bool> assign_seq(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                             const std::vector<Kmer>& head_kmer_arr,
                                             RefindPathVector& unitig_complete_paths,
                                             const int kmer,
                                             const bool is_ref,
                                             const fm_index_coll& fm_idx,
                                             const std::vector<int>& curr_nodelist,
                                             const bool downstream);

RefindPathVector iter_nodes_length (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                  const std::vector<Kmer>& head_kmer_arr,
                                  const NodeTuple& head_node_tuple,
                                  const size_t& current_colour,
                                  const size_t& radius,
                                  const bool& repeat,
                                  const bool& is_ref,
                                  const fm_index_coll& fm_idx,
                                  const int overlap,
                                  const std::unordered_set<int>& to_avoid);

RefindTuple traverse_outward(const ColoredCDBG<MyUnitigMap>& ccdbg,
                             const std::vector<Kmer>& head_kmer_arr,
                             const size_t& colour_ID,
                             const std::pair<std::vector<int>, std::vector<indexPair>>& ORF_info_const,
                             const size_t& radius,
                             const bool is_ref,
                             const int kmer,
                             const fm_index_coll& fm_idx,
                             const bool repeat,
                             const std::unordered_set<int>& to_avoid);

RefindMap refind_in_nodes(const ColoredCDBG<MyUnitigMap>& ccdbg,
                          const std::vector<Kmer>& head_kmer_arr,
                          const size_t colour_ID,
                          const NodeSearchDict& node_search_dict,
                          const size_t radius,
                          const bool is_ref,
                          const int kmer,
                          const fm_index_coll& fm_idx,
                          const bool repeat,
                          const std::unordered_set<int>& to_avoid);

#endif //GENE_REFINDING_H