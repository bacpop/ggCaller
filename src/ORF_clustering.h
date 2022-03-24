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
                               const ORFMatrixMap& ORF_mat_map,
                               const std::vector<std::unordered_set<size_t>>& ORF_group_vector,
                               std::vector<int>& centroid_vector,
                               const robin_hood::unordered_map<size_t, size_t>& ID_hash_map,
                               const std::vector<std::pair<size_t, size_t>>& ORF_length_list,
                               const double& id_cutoff,
                               const double& len_diff_cutoff);

double align_seqs(const ORFNodeVector& ORF1_info,
                  const ORFNodeVector& ORF2_info,
                  const ColoredCDBG<MyUnitigMap>& ccdbg,
                  const std::vector<Kmer>& head_kmer_arr,
                  const size_t& DBG_overlap);

#endif //GGCALLER_ORF_CLUSTERING_H
