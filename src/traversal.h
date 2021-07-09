#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "unitigDict.h"
#include "concurrent_hash_map.h"
#include "indexing.h"

// traversal.cpp
void iter_nodes_binary (const GraphVector& graph_vector,
                      NodeColourMap& node_colour_vector_traversed,
                      std::vector<AllPaths>& colour_graph_paths,
                      const NodeTuple& head_node_tuple,
                      const size_t& current_colour,
                      const size_t& length_max,
                      const bool& repeat);

//AllPaths traverse_graph(const GraphVector& graph_vector,
//                         const size_t& colour_ID,
//                         const std::vector<size_t>& node_ids,
//                         const bool repeat,
//                         const size_t max_path_length);

#endif //TRAVERSAL_H
