#ifndef GENE_REFINDING_H
#define GENE_REFINDING_H

#include "unitigDict.h"
#include "match_string.h"
#include "ORF_clustering.h"

void assign_seq(const GraphVector& graph_vector,
                const PathVector& unitig_complete_paths,
                const int kmer,
                const bool is_ref,
                const fm_index_coll& fm_idx,
                std::string& stream_seq,
                const size_t& ORF_end,
                const std::string& ORF_seq);

PathVector iter_nodes_length (const GraphVector& graph_vector,
                              const NodeTuple& head_node_tuple,
                              const size_t& current_colour,
                              const size_t& radius,
                              const bool& repeat);

std::string traverse_outward(const GraphVector& graph_vector,
                             const size_t& colour_ID,
                             const ORFNodeVector& ORF_info,
                             const size_t& radius,
                             bool is_ref,
                             const bool write_idx,
                             const int kmer,
                             const std::string& FM_fasta_file,
                             const bool repeat);

#endif //GENE_REFINDING_H