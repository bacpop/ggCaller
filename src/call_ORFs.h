#ifndef BIFROST_API_GGCALLER_H
#define BIFROST_API_GGCALLER_H

#include "ORF_scoring.h"
#include "unitigDict.h"
#include "match_string.h"
#include "indexing.h"


// call_ORFs
void generate_ORFs(const int& colour_ID,
                   ORFNodeMap& ORF_node_map,
                   const ColoredCDBG<MyUnitigMap>& ccdbg,
                   const std::vector<Kmer>& head_kmer_arr,
                   const std::vector<std::string>& stop_codons,
                   const std::vector<std::string>& start_codons,
                   const std::vector<int>& nodelist,
                   const int& overlap,
                   const size_t min_len,
                   const bool is_ref,
                   const fm_index_coll& fm_idx,
                   torch::jit::script::Module& ORF_model,
                   torch::jit::script::Module& TIS_model,
                   const double& minimum_ORF_score,
                   const bool no_filter,
                   tbb::concurrent_unordered_map<size_t, double>& all_ORF_scores,
                   tbb::concurrent_unordered_map<size_t, double>& all_TIS_scores);

ORFCoords calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                           const std::vector<int>& nodelist,
                           const std::vector<std::vector<size_t>>& node_ranges);

ORFVector call_ORFs(const int colour_ID,
                    std::vector<PathVector>& all_paths,
                    const ColoredCDBG<MyUnitigMap>& ccdbg,
                    const std::vector<Kmer>& head_kmer_arr,
                    const std::vector<std::string>& stop_codons_for,
                    const std::vector<std::string>& start_codons_for,
                    const int overlap,
                    const size_t min_ORF_length,
                    const bool is_ref,
                    const fm_index_coll& fm_idx,
                    torch::jit::script::Module& ORF_model,
                    torch::jit::script::Module& TIS_model,
                    const double& minimum_ORF_score,
                    const bool no_filter,
                    tbb::concurrent_unordered_map<size_t, double>& all_ORF_scores,
                    tbb::concurrent_unordered_map<size_t, double>& all_TIS_scores);

ORFVector sort_ORF_indexes(ORFNodeMap& ORF_node_map,
                           const NodeStrandMap& pos_strand_map,
                           const ColoredCDBG<MyUnitigMap>& ccdbg,
                           const std::vector<Kmer>& head_kmer_arr,
                           const bool is_ref);

NodeStrandMap calculate_pos_strand(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                   const std::vector<Kmer>& head_kmer_arr,
                                   const ORFNodeMap& ORF_node_map);

void update_ORF_node_map (const ColoredCDBG<MyUnitigMap>& ccdbg,
                          const std::vector<Kmer>& head_kmer_arr,
                          const size_t& ORF_hash,
                          ORFNodeVector& ORF_node_vector,
                          ORFNodeMap& ORF_node_map);

#endif //BIFROST_API_GGCALLER_H