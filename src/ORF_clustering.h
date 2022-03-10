//
// Created by sth19 on 01/09/2021.
//

#ifndef GGCALLER_ORF_CLUSTERING_H
#define GGCALLER_ORF_CLUSTERING_H

#include "unitigDict.h"
#include "gene_overlap.h"
#include "translation.h"

ORFGroupTuple group_ORFs(const ColourORFVectorMap& colour_ORF_map,
                         const std::vector<Kmer>& head_kmer_arr);

ORFClusterMap produce_clusters(const ColourORFVectorMap& colour_ORF_map,
                               const ColoredCDBG<MyUnitigMap>& ccdbg,
                               const std::vector<Kmer>& head_kmer_arr,
                               const size_t& DBG_overlap,
                               const ORFMatrixVector& ORF_mat_vector,
                               const std::vector<std::unordered_set<size_t>>& ORF_group_vector,
                               std::vector<std::vector<std::pair<size_t, size_t>>>& centroid_vector,
                               const double& id_cutoff,
                               const double& len_diff_cutoff);

double align_seqs(const ORFNodeVector& ORF1_info,
                  const ORFNodeVector& ORF2_info,
                  const ColoredCDBG<MyUnitigMap>& ccdbg,
                  const std::vector<Kmer>& head_kmer_arr,
                  const size_t& DBG_overlap);

#endif //GGCALLER_ORF_CLUSTERING_H
