//
// Created by sth19 on 28/01/2022.
//

#ifndef GGCALLER_GENE_GRAPH_H
#define GGCALLER_GENE_GRAPH_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include "definitions.h"

using namespace boost;
using namespace std;

// balrog scores for overlap penalties
const float unidirectional_penalty_per_base = 3.895921717182765;  // 3' 5' overlap
const float convergent_penalty_per_base = 4.603432608883688;  // 3' 3' overlap
const float divergent_penalty_per_base = 3.3830814940689975;  // 5' 5' overlap

// typedefs for graph and edges for boost graph library
typedef property<edge_weight_t, float> Weight;
typedef property<vertex_index_t, size_t> VertexIndex;
typedef adjacency_list<vecS, vecS, bidirectionalS, no_property, Weight> GeneGraph;
typedef std::pair<int, int> Edge;
typedef graph_traits<GeneGraph>::vertex_descriptor VertexDescriptor;
typedef graph_traits<GeneGraph>::edge_descriptor EdgeDescriptor;
typedef graph_traits<GeneGraph>::out_edge_iterator OutEdgeIterator;
typedef graph_traits<GeneGraph>::in_edge_iterator InEdgeIterator;

// cycle detector struct
struct cycle_detector : public dfs_visitor<>
{
    cycle_detector( bool& has_cycle, EdgeDescriptor& e)
            : _has_cycle(has_cycle), _cycle_e(e) { }

    template <class Edge, class Graph>
    void back_edge(Edge e, Graph& g) {
        _has_cycle = true;
        _cycle_e = e;
    }
protected:
    bool& _has_cycle;
    EdgeDescriptor& _cycle_e;
};

std::vector<size_t> getPath(
        const GeneGraph& graph,
        const std::vector<VertexDescriptor>& pMap,
        const std::vector<float>& distances,
        const VertexDescriptor& source,
        const VertexDescriptor& destination,
        float& path_score,
        const std::vector<size_t>& vertex_mapping);

template <class T>
std::vector<size_t> traverse_components(const ORFVector& ORF_vector,
                                        const std::vector<size_t>& vertex_mapping,
                                        const std::unordered_set<size_t>& vertex_list,
                                        const GeneGraph& g,
                                        const float& minimum_path_score,
                                        const size_t numVertices,
                                        T weight_pmap);

std::vector<std::vector<size_t>> call_true_genes (const ORFVector& ORF_vector,
                                                  const ORFOverlapMap& overlap_map,
                                                  const float& minimum_path_score);

#endif //GGCALLER_GENE_GRAPH_H
