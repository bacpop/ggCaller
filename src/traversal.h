#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "unitigDict.h"
#include "indexing.h"
#include "gene_overlap.h"
#include "ORF_connection.h"

PathVector iter_nodes_binary (const GraphVector& graph_vector,
                              const NodeTuple& head_node_tuple,
                              const size_t& current_colour,
                              const size_t& length_max,
                              const bool& repeat);

std::vector<PathVector> traverse_graph(const GraphVector& graph_vector,
                                      const size_t& colour_ID,
                                      const std::vector<size_t>& node_ids,
                                      const bool repeat,
                                      const size_t max_path_length);

#endif //TRAVERSAL_H
