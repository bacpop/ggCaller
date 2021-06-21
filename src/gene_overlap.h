#ifndef GENE_OVERLAP_H
#define GENE_OVERLAP_H

#include "unitigDict.h"

// gene_overlap.cpp
ORFOverlapMap calculate_overlaps(const GraphVector& graph_vector,
                                 const std::pair<ORFVector, NodeStrandMap>& ORF_pair,
                                 const int DBG_overlap,
                                 const size_t max_overlap);

#endif //BIFROST_API_GGCALLER_H
