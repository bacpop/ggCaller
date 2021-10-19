#ifndef BIFROST_API_GGCALLER_H
#define BIFROST_API_GGCALLER_H

#include "unitigDict.h"
#include "match_string.h"
#include "indexing.h"

// call_ORFs
void generate_ORFs(ORFNodeMap& ORF_node_map,
                   std::unordered_set<size_t>& hashes_to_remove,
                   const GraphVector& graph_vector,
                   const std::vector<std::string>& stop_codons,
                   const std::vector<std::string>& start_codons,
                   const std::vector<int>& unitig_path,
                   const int& overlap,
                   const size_t min_len,
                   const bool is_ref,
                   const fm_index_coll& fm_idx);

ORFCoords calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                           const std::vector<int>& nodelist,
                           const std::vector<std::vector<size_t>>& node_ranges,
                           const int& overlap);

ORFVector call_ORFs(const std::vector<PathVector>& all_paths,
                     const GraphVector& graph_vector,
                     const std::vector<std::string>& stop_codons_for,
                     const std::vector<std::string>& start_codons_for,
                     const int overlap,
                     const size_t min_ORF_length,
                     const bool is_ref,
                     const fm_index_coll& fm_idx);

ORFVector sort_ORF_indexes(ORFNodeMap& ORF_node_map,
                           const NodeStrandMap& pos_strand_map);

NodeStrandMap calculate_pos_strand(const ORFNodeMap& ORF_node_map);

void update_ORF_node_map (const GraphVector& graph_vector,
                          const size_t& ORF_hash,
                          ORFNodeVector& ORF_node_vector,
                          ORFNodeMap& ORF_node_map);

#endif //BIFROST_API_GGCALLER_H
