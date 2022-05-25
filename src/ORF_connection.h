#ifndef ORF_CONNECTION_H
#define ORF_CONNECTION_H

#include "unitigDict.h"
#include "match_string.h"

robin_hood::unordered_map<size_t, std::vector<size_t>> add_ORF_info (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                     const std::vector<Kmer>& head_kmer_arr,
                                                                     const robin_hood::unordered_set<size_t>& target_ORFs,
                                                                     const ORFNodeMap& gene_map);

std::set<std::pair<size_t, size_t>> pair_ORF_nodes (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                    const std::vector<Kmer>& head_kmer_arr,
                                                    robin_hood::unordered_map<size_t, std::vector<size_t>>& node_to_ORFs,
                                                    const size_t colour_ID,
                                                    const robin_hood::unordered_set<size_t>& target_ORFs,
                                                    const ORFNodeMap& gene_map,
                                                    const size_t& max_ORF_path_length,
                                                    const bool repeat,
                                                    const int stream,
                                                    robin_hood::unordered_set<size_t>& downstream_ORF_set,
                                                    robin_hood::unordered_set<size_t>& upstream_ORF_set,
                                                    const int overlap,
                                                    const bool is_ref,
                                                    const fm_index_coll& fm_idx);

std::set<std::pair<size_t, size_t>> check_next_ORFs (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                     const std::vector<Kmer>& head_kmer_arr,
                                                     const robin_hood::unordered_map<size_t, std::vector<size_t>>& node_to_ORFs,
                                                     const int& head_node,
                                                     const size_t stream_source,
                                                     const size_t current_colour,
                                                     const int stream,
                                                     const ORFNodeMap& gene_map,
                                                     const size_t max_ORF_path_length,
                                                     const bool repeat,
                                                     robin_hood::unordered_set<size_t>& downstream_ORF_set,
                                                     robin_hood::unordered_set<size_t>& upstream_ORF_set,
                                                     const int overlap,
                                                     const bool is_ref,
                                                     const fm_index_coll& fm_idx);

#endif //ORF_CONNECTION_H
