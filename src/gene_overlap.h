#ifndef GENE_OVERLAP_H
#define GENE_OVERLAP_H

#include "unitigDict.h"

void reverse_ORFNodeVector(const GraphVector& graph_vector,
                           std::pair<std::vector<int>, std::vector<indexPair>>& ORF2_nodes,
                           int& ORF2_start_node,
                           int& ORF2_end_node,
                           std::pair<int, size_t>& ORF2_5p,
                           std::pair<int, size_t>& ORF2_3p);

std::tuple<bool, std::vector<size_t>, std::vector<size_t>> slice_ORFNodeVector(const ORFNodeVector& ORF1_nodes,
                                                                               const std::pair<std::vector<int>, std::vector<indexPair>>& ORF2_nodes,
                                                                               const int& ORF2_start_node,
                                                                               const int& ORF2_end_node);

ORFOverlapMap calculate_overlaps(const GraphVector& graph_vector,
                                 const ORFVector& ORF_vector,
                                 const int DBG_overlap,
                                 const size_t max_overlap);


#endif //BIFROST_API_GGCALLER_H
