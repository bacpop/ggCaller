#ifndef ORF_CONNECTION_H
#define ORF_CONNECTION_H

#include "unitigDict.h"

std::vector<robin_hood::unordered_set<size_t>> add_ORF_info (const std::vector<Kmer>& head_kmer_arr,
                                                             const std::set<size_t>& target_ORFs,
                                                             const ORFVector& ORF_vector);

//void remove_ORF_info (ColoredCDBG<MyUnitigMap>& ccdbg,
//                      const std::vector<Kmer>& head_kmer_arr,
//                      const size_t& colour_ID,
//                      const std::vector<size_t>& target_ORFs,
//                      const ORFVector& ORF_vector);

std::vector<size_t> order_ORFs_in_node(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                       const std::vector<Kmer>& head_kmer_arr,
                                       const robin_hood::unordered_set<size_t>& node_ORFs,
                                       const int node_id,
                                       const ORFVector& ORF_vector);

std::vector<std::pair<size_t, size_t>> pair_ORF_nodes (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                       const std::vector<Kmer>& head_kmer_arr,
                                                       const std::vector<robin_hood::unordered_set<size_t>> node_to_ORFs,
                                                       const size_t colour_ID,
                                                       const std::set<size_t>& target_ORFs,
                                                       const ORFVector& ORF_vector,
                                                       const size_t& max_ORF_path_length,
                                                       const int stream,
                                                       std::unordered_set<int>& prev_node_set,
                                                       const bool is_ref,
                                                       const int overlap);

std::vector<std::pair<size_t, size_t>> check_next_ORFs (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                        const std::vector<Kmer>& head_kmer_arr,
                                                        const std::vector<robin_hood::unordered_set<size_t>> node_to_ORFs,
                                                        const int& head_node,
                                                        const size_t stream_source,
                                                        const size_t current_colour,
                                                        const int stream,
                                                        const ORFVector& ORF_vector,
                                                        const size_t max_ORF_path_length,
                                                        std::unordered_set<int>& prev_node_set,
                                                        const bool is_ref,
                                                        const int overlap);

#endif //ORF_CONNECTION_H
