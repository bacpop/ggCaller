#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "unitigDict.h"

// traversal.cpp
PathVector recur_nodes_binary (const UnitigVector& graph_vector,
                               const std::vector<int>& head_kmer_list,
                               const uint8_t& codon_arr,
                               const size_t& colour_ID,
                               const std::unordered_set<int>& kmer_set,
                               const size_t& length,
                               const size_t& length_max,
                               const bool& repeat);

AllPaths traverse_graph(const UnitigVector& graph_vector,
                         const size_t& colour_ID,
                         const std::vector<size_t>& node_ids,
                         const bool repeat,
                         const size_t max_path_length);

// match_strings
fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx);

int seq_search(const seqan3::dna5_vector& query,
               const fm_index_coll& ref_idx);

bool check_colours(const std::string& path_sequence,
                   const fm_index_coll& fm_idx);

#endif //TRAVERSAL_H
