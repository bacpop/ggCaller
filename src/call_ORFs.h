#ifndef CALL_ORFS_H
#define CALL_ORFS_H

#include "unitigDict.h"
#include "concurrent_hash_map.h"
#include "match_string.h"
#include "indexing.h"

// call_ORFs
void generate_ORFs(std::vector<ORFNodeMap>& total_ORF_node_map,
                         const GraphVector& graph_vector,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const size_t& colour_ID,
                         const std::vector<bool>& colours,
                         const std::vector<bool>& colour_complete,
                         const std::vector<int>& path,
                         const int& overlap,
                         const size_t min_len,
                         const bool is_ref,
                         const std::vector<fm_index_coll>& fm_idx_vector);

std::tuple<std::string, std::vector<int>, std::vector<indexPair>> calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                                                                                    const std::vector<int>& nodelist,
                                                                                    const std::vector<std::vector<size_t>>& node_ranges,
                                                                                    const int& overlap);

std::pair<ORFVector, NodeStrandMap> call_ORFs(const AllPaths& all_paths,
                                             const GraphVector& graph_vector,
                                             const std::vector<std::string>& stop_codons_for,
                                             const std::vector<std::string>& start_codons_for,
                                             const int overlap,
                                             const size_t min_ORF_length,
                                             const bool is_ref,
                                             const fm_index_coll& fm_idx);

ORFVector sort_ORF_indexes(ORFNodeMap& ORF_node_map);

NodeStrandMap calculate_pos_strand(const ORFNodeMap& ORF_node_map);

#endif //CALL_ORFS_H
