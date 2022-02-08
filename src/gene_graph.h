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
const double unidirectional_penalty_per_base = 3.895921717182765;  // 3' 5' overlap
const double convergent_penalty_per_base = 4.603432608883688;  // 3' 3' overlap
const double divergent_penalty_per_base = 3.3830814940689975;  // 5' 5' overlap

// typedefs for graph and edges for boost graph library
typedef property<edge_weight_t, double> Weight;
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

std::vector<VertexDescriptor> getPath(
        const GeneGraph& graph,
        const std::vector<VertexDescriptor>& pMap,
        const std::vector<double>& distances,
        const VertexDescriptor& source,
        const VertexDescriptor& destination,
        double& path_score,
        const std::vector<size_t>& vertex_mapping);

std::vector<VertexDescriptor> traverse_components(const std::unordered_map<size_t, double>& score_map,
                                                  const std::unordered_set<size_t>& vertex_list,
                                                  const GeneGraph& g,
                                                  const double& minimum_path_score,
                                                  const size_t numVertices);

std::vector<std::vector<size_t>> call_true_genes (const std::unordered_map<size_t, double>& score_map,
                                                  const ORFOverlapMap& overlap_map,
                                                  const double& minimum_path_score);

#endif //GGCALLER_GENE_GRAPH_H
