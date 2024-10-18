#ifndef GENE_OVERLAP_H
#define GENE_OVERLAP_H

#include "unitigDict.h"
#include "match_string.h"

inline std::pair<std::vector<int>, std::vector<indexPair>> combine_nodes(const std::vector<int>& ORF_nodes,
                                                                         const std::vector<indexPair>& ORF_coords,
                                                                         const std::vector<int>& TIS_nodes,
                                                                         const std::vector<indexPair>& TIS_coords);

void reverse_ORFNodeVector(const ColoredCDBG<MyUnitigMap>& ccdbg,
                           const std::vector<Kmer>& head_kmer_arr,
                           std::pair<std::vector<int>, std::vector<indexPair>>& ORF2_nodes,
                           int& ORF2_start_node,
                           int& ORF2_end_node,
                           std::pair<int, size_t>& ORF2_5p,
                           std::pair<int, size_t>& ORF2_3p);

std::tuple<bool, std::vector<size_t>, std::vector<size_t>> slice_ORFNodeVector(const std::pair<std::vector<int>, std::vector<indexPair>>& ORF1_nodes,
                                                                               const std::pair<std::vector<int>, std::vector<indexPair>>& ORF2_nodes,
                                                                               const int& ORF2_start_node,
                                                                               const int& ORF2_end_node,
                                                                               const bool is_ref,
                                                                               const fm_index_coll& fm_idx);

ORFOverlapMap calculate_overlaps(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                 const std::vector<Kmer>& head_kmer_arr,
                                 const ORFNodeMap& ORF_map,
                                 const int DBG_overlap,
                                 const size_t max_overlap,
                                 const bool is_ref,
                                 const fm_index_coll& fm_idx);


#endif //BIFROST_API_GGCALLER_H
