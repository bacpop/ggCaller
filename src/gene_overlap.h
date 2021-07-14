#ifndef GENE_OVERLAP_H
#define GENE_OVERLAP_H

#include "unitigDict.h"

// gene_overlap.cpp
ORFOverlapMap calculate_overlaps(const GraphVector& graph_vector,
                                 const std::pair<ORFVector, NodeStrandMap>& ORF_pair,
                                 const int DBG_overlap,
                                 const size_t max_overlap);

// checks how ORF ends are ordered if they share a common node
std::vector<std::pair<size_t,size_t>> order_node_ends(const GraphVector& graph_vector,
                                                      const std::unordered_set<size_t>& node_ORFs,
                                                      const int& id);

#endif //BIFROST_API_GGCALLER_H
