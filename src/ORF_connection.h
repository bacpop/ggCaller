#ifndef ORF_CONNECTION_H
#define ORF_CONNECTION_H

#include "unitigDict.h"
#include "traversal.h"

//PathOverlapMap overlapping_paths(const GraphVector& graph_vector,
//                                 const PathVector& all_paths);
//
//std::vector<std::pair<size_t, size_t>> pair_ORF_paths (const PathOverlapMap& path_overlap_map,
//                                                       const std::unordered_map<size_t, std::vector<size_t>>& ORF_path_map,
//                                                       const ORFVector& ORF_vector,
//                                                       const std::unordered_set<size_t>& target_ORFs,
//                                                       std::unordered_set<size_t>& unpaired_downstream,
//                                                       std::unordered_set<size_t> unpaired_upstream);
//
//void order_ORFs_in_path(std::vector<size_t>& ORF_list,
//                        const ORFVector& ORF_vector);
//
//std::vector<size_t> check_next_path (const PathOverlapMap& path_overlap_map,
//                                     const size_t& head_node,
//                                     const std::unordered_map<size_t, std::vector<size_t>>& ORF_path_map,
//                                     const ORFVector& ORF_vector,
//                                     const int& stream);

void add_ORF_info (GraphVector& graph_vector,
                  const size_t& colour_ID,
                  const std::vector<size_t>& target_ORFs,
                  const ORFVector& ORF_vector);

std::vector<size_t> order_ORFs_in_node(const GraphVector& graph_vector,
                                       const std::unordered_set<size_t>& node_ORFs,
                                       const int& node_id,
                                       const ORFVector& ORF_vector);

std::vector<std::pair<size_t, size_t>> pair_ORF_nodes (const GraphVector& graph_vector,
                                                      const size_t& colour_ID,
                                                      const std::vector<size_t>& target_ORFs,
                                                      const ORFVector& ORF_vector,
                                                      const size_t& max_ORF_path_length,
                                                      const int& stream,
                                                      std::unordered_set<int>& prev_node_set);

std::vector<std::pair<size_t, size_t>> check_next_ORFs (const GraphVector& graph_vector,
                                                        const int& head_node,
                                                        const size_t& stream_source,
                                                        const size_t& current_colour,
                                                        const int& stream,
                                                        const ORFVector& ORF_vector,
                                                        const size_t& max_ORF_path_length,
                                                        std::unordered_set<int>& prev_node_set);

#endif //ORF_CONNECTION_H
