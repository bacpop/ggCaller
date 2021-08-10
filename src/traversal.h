#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "unitigDict.h"
#include "indexing.h"
#include "gene_overlap.h"

PathVector iter_nodes_binary (const GraphVector& graph_vector,
                              const NodeTuple& head_node_tuple,
                              const size_t& current_colour,
                              const size_t& length_max,
                              const bool& repeat);

AllPaths traverse_graph(const GraphVector& graph_vector,
                         const size_t& colour_ID,
                         const std::vector<size_t>& node_ids,
                         const bool repeat,
                         const size_t max_path_length);

std::vector<std::pair<size_t, size_t>> check_next_ORFs (const GraphVector& graph_vector,
                                                        const int& head_node,
                                                        const size_t& stream_source,
                                                        const size_t& current_colour,
                                                        const int& stream,
                                                        const ORFVector& ORF_vector,
                                                        std::unordered_set<int>& prev_node_set);

template <class T>
void clear_stack(std::stack<T>& to_clear);

#endif //TRAVERSAL_H
