#ifndef GENE_REFINDING_H
#define GENE_REFINDING_H

#include "unitigDict.h"
#include "match_string.h"
#include "ORF_clustering.h"

std::vector<std::vector<size_t>> calculate_node_ranges(const GraphVector& graph_vector,
                                                       const int& overlap,
                                                       const std::vector<int>& full_nodelist);

std::pair<std::vector<int>, std::pair<ContigLoc, bool>> assign_seq(const size_t& colour_ID,
                                                                   const GraphVector& graph_vector,
                                                                   const PathVector& unitig_complete_paths,
                                                                   const int kmer,
                                                                   const bool is_ref,
                                                                   const fm_index_coll& fm_idx,
                                                                   std::string& stream_seq,
                                                                   const size_t& ORF_end,
                                                                   const std::string& ORF_seq);

PathVector iter_nodes_length (const GraphVector& graph_vector,
                              const NodeTuple& head_node_tuple,
                              const size_t& current_colour,
                              const size_t& radius,
                              const bool& repeat,
                              const bool& is_ref,
                              const fm_index_coll& fm_idx);

RefindTuple traverse_outward(const GraphVector& graph_vector,
                             const size_t& colour_ID,
                             const ORFNodeVector& ORF_info,
                             const size_t& radius,
                             const bool is_ref,
                             const int kmer,
                             const fm_index_coll& fm_idx,
                             const std::vector<size_t>& contig_locs,
                             const bool repeat);

RefindMap refind_in_nodes(const GraphVector& graph_vector,
                          const size_t& colour_ID,
                          const std::unordered_map<int, std::unordered_map<std::string, ORFNodeVector>>& node_search_dict,
                          const size_t& radius,
                          const bool is_ref,
                          const int kmer,
                          const fm_index_coll& fm_idx,
                          const bool repeat);

#endif //GENE_REFINDING_H