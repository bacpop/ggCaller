#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "unitigDict.h"
#include "indexing.h"
#include "gene_overlap.h"
#include "ORF_connection.h"

PathVector iter_nodes_binary (const ColoredCDBG<MyUnitigMap>& ccdbg,
                              const std::vector<Kmer>& head_kmer_arr,
                              const NodeTuple& head_node_tuple,
                              const size_t current_colour,
                              const size_t length_max,
                              const size_t overlap,
                              const bool repeat,
                              const bool is_ref,
                              const boost::dynamic_bitset<>& ref_set,
                              const fm_index_coll& fm_idx);

std::vector<PathVector> traverse_graph(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                       const std::vector<Kmer>& head_kmer_arr,
                                       const size_t colour_ID,
                                       const std::vector<size_t>& node_ids,
                                       const bool repeat,
                                       const size_t max_path_length,
                                       const size_t overlap,
                                       const bool is_ref,
                                       const boost::dynamic_bitset<>& ref_set,
                                       const fm_index_coll& fm_idx);

#endif //TRAVERSAL_H
