//
// Created by sth19 on 28/01/2022.
//

#include "gene_graph.h"

// recursively look for adjacent vertices
void find_adjacent(const VertexDescriptor& v,
                   const GeneGraph& g,
                   std::unordered_set<size_t>& curr_component)
{
    if (curr_component.find(v) == curr_component.end())
    {
        curr_component.insert(v);

        // iterate over out edges
        {
            OutEdgeIterator out_begin, out_end;
            boost::tie(out_begin, out_end) = boost::out_edges(v, g);
            for (; out_begin != out_end; out_begin++)
            {
                find_adjacent(target(*out_begin, g), g, curr_component);
            }
        }

        // iterate over in edges
        {
            InEdgeIterator in_begin, in_end;
            boost::tie(in_begin, in_end) = boost::in_edges(v, g);
            for (; in_begin != in_end; in_begin++)
            {
                find_adjacent(source(*in_begin, g), g, curr_component);
            }
        }
    }
}

// iterate over all vertices, finding all adjacent vertices recursively
std::vector<std::unordered_set<size_t>> get_components(const GeneGraph& g)
{
    std::unordered_set<size_t> visited;

    std::vector<std::unordered_set<size_t>> components;

    GeneGraph::vertex_iterator i, end;

    for (boost::tie(i, end) = boost::vertices(g); i != end; i++) {
        if (visited.find(*i) == visited.end())
        {
            visited.insert(*i);

            std::unordered_set<size_t> curr_component {*i};

            // iterate over out edges
            {
                OutEdgeIterator out_begin, out_end;
                boost::tie(out_begin, out_end) = boost::out_edges(*i, g);
                for (; out_begin != out_end; out_begin++)
                {
                    find_adjacent(target(*out_begin, g), g, curr_component);
                }
            }

            // iterate over in edges
            {
                InEdgeIterator in_begin, in_end;
                boost::tie(in_begin, in_end) = boost::in_edges(*i, g);
                for (; in_begin != in_end; in_begin++)
                {
                    find_adjacent(source(*in_begin, g), g, curr_component);
                }
            }

            visited.insert(curr_component.begin(), curr_component.end());

            components.push_back(std::move(curr_component));
        }
    }

    return components;
}

std::vector<size_t> getPath(
        const GeneGraph& graph,
        const std::vector<VertexDescriptor>& pMap,
        const std::vector<float>& distances,
        const VertexDescriptor& source,
        const VertexDescriptor& destination,
        float& path_score,
        const std::vector<size_t>& vertex_mapping,
        const std::vector<size_t>& component_vertex_mapping)
{
    std::vector<VertexDescriptor> path;

    // no path to destination
    if (pMap.at(destination) == destination)
    {
        path_score = 0;
        return path;
    }

    // if path to target, get path out
    VertexDescriptor current = destination;
    while (current != source)
    {
        path.push_back(vertex_mapping.at(component_vertex_mapping.at(current)));
        path_score += distances.at(current);
        current = pMap.at(current);
    }
    path.push_back(vertex_mapping.at(component_vertex_mapping.at(source)));

    // reverse as works back from beginning
    std::reverse(path.begin(), path.end());
    return path;
}

template <class T>
std::vector<size_t> traverse_components(const ORFNodeMap& ORF_map,
                                       const std::vector<size_t>& vertex_mapping,
                                       const std::vector<size_t>& component_vertex_mapping,
                                       const std::unordered_set<size_t>& vertex_list,
                                       const GeneGraph& g,
                                       const float& minimum_path_score,
                                       const size_t numVertices,
                                       T weight_pmap,
                                       const std::vector<Kmer>& head_kmer_arr)
{
    // set all vertices as source and sink
    std::vector<VertexDescriptor> start_vertices;
    std::vector<VertexDescriptor> end_vertices;

    // push all vertices onto vertex list to test all possible paths
    for (const auto& v : vertex_list)
    {
        start_vertices.push_back(v);
        end_vertices.push_back(v);
    }

    // catch issue if no vertices in start or end
    if (start_vertices.empty())
    {
        std::copy(vertex_list.begin(), vertex_list.end(), start_vertices.begin());
    }
    if (end_vertices.empty())
    {
        std::copy(vertex_list.begin(), vertex_list.end(), end_vertices.begin());
    }

    // set values of paths and scores
    std::vector<size_t> gene_path;
    float high_score = 0;

    for (const auto& start : start_vertices)
    {
        // get start score
        const float start_score = -std::get<4>(ORF_map.at(vertex_mapping.at(component_vertex_mapping.at(start))));

        // if starting node is highest scorer, set as path
        if (start_score < high_score)
        {
            high_score = start_score;
            gene_path = {vertex_mapping.at(component_vertex_mapping.at(start))};
        }

        std::vector<float> distances(numVertices);
        std::vector<VertexDescriptor> pMap(numVertices);

        // call to the algorithm, searching all paths from start
        bellman_ford_shortest_paths(g, numVertices,
                                    weight_map(weight_pmap).
                                            root_vertex(start).
                                            predecessor_map(make_iterator_property_map(pMap.begin(), get(vertex_index, g))).
                                            distance_map(make_iterator_property_map(distances.begin(), get(vertex_index, g))));

        // determine shortest path to each end
        for (const auto& end : end_vertices)
        {
            // pass if start and end are the same
            if (start == end)
            {
                continue;
            }

            float path_score = start_score;

            auto path = getPath(g, pMap, distances, start, end, path_score, vertex_mapping, component_vertex_mapping);

            if (path_score < high_score)
            {
                gene_path = std::move(path);
                high_score = path_score;
            }
            // if path_score is the same, take longest path
            else if (path_score == high_score)
            {
                if (path.size() > gene_path.size())
                {
                    gene_path = std::move(path);
                    high_score = path_score;
                }
                // if same size, then take lowest hash
                else if (path.size() == gene_path.size())
                {
                    size_t curr_hash;
                    {
                        std::string seq;
                        for (const auto& ORF_ID : gene_path)
                        {
                            const auto& ORF_info = ORF_map.at(ORF_ID);
                            for (const auto& node_traversed : std::get<0>(ORF_info))
                            {
                                const size_t node_ID = abs(node_traversed) - 1;
                                seq += head_kmer_arr.at(node_ID).toString();
                            }
                        }
                        curr_hash = hasher{}(seq);
                    }

                    size_t new_hash;
                    {
                        std::string seq;
                        for (const auto& ORF_ID : path)
                        {
                            const auto& ORF_info = ORF_map.at(ORF_ID);
                            for (const auto& node_traversed : std::get<0>(ORF_info))
                            {
                                const size_t node_ID = abs(node_traversed) - 1;
                                seq += head_kmer_arr.at(node_ID).toString();
                            }
                        }
                        new_hash = hasher{}(seq);
                    }

                    if (new_hash < curr_hash)
                    {
                        gene_path = std::move(path);
                        high_score = path_score;
                    }
                }
            }
        }
    }

    // check if higher than minimum_path_score
    if (high_score <= -(minimum_path_score))
    {
        return gene_path;
    } else
    {
        return {};
    }
}

std::vector<std::vector<size_t>> call_true_genes (const ORFNodeMap& ORF_map,
                                                  const ORFOverlapMap& overlap_map,
                                                  const float& minimum_path_score,
                                                  const std::vector<Kmer>& head_kmer_arr)
{
    std::vector<std::vector<size_t>> gene_paths;

    // scope for graph objects
    {
        // initialise graph
        GeneGraph g;
        std::vector<size_t> vertex_mapping;
        std::unordered_map<size_t, size_t> ORF_ID_mapping;
        std::vector<std::unordered_set<size_t>> components;
        size_t numVertices;

        // create IDs for vertex mappings
        int vertex_id = 0;
        // add all vertices
        for (const auto& entry : ORF_map)
        {
            vertex_mapping.push_back(entry.first);
            ORF_ID_mapping[entry.first] = vertex_id++;
            add_vertex(g);
        }

        // add all edges
        for (const auto& target : overlap_map)
        {
            for (const auto& source : target.second)
            {
                add_edge(ORF_ID_mapping.at(source.first), ORF_ID_mapping.at(target.first), g);
            }
        }

        // check for cycles
        bool has_cycle = false;
        EdgeDescriptor cycle_edge;
        cycle_detector vis(has_cycle, cycle_edge);
        depth_first_search(g, visitor(vis));

        // check if cycle present, if so continually remove until no more cycles
        while (has_cycle)
        {
            remove_edge(cycle_edge, g);
            has_cycle = false;
            depth_first_search(g, visitor(vis));
        }

        // get components from simple graph
        components = get_components(g);
        numVertices = num_vertices(g);

        // iterate over components using bellman ford algorithm
        for (const auto& component : components)
        {
            // make transative closure graph
            GeneGraph tc;
            // TODO this won't work with a component, need to generate a small graph as above most likely

            std::vector<size_t> component_vertex_mapping;
            std::unordered_map<size_t, size_t> component_ORF_ID_mapping;
            {
                GeneGraph component_g;
                int component_vertex_id = 0;
                for (const auto& entry : component)
                {
                    component_vertex_mapping.push_back(entry);
                    component_ORF_ID_mapping[entry] = component_vertex_id++;
                    add_vertex(component_g);
                }

                //iterate over component again and add edges if in original graph
                for (const auto& source : component)
                {
                    for (const auto& target : component)
                    {
                        auto edgePair = boost::edge(source, target, g);
                        if (edgePair.second)
                        {
                            add_edge(component_ORF_ID_mapping.at(source), component_ORF_ID_mapping.at(target), component_g);
                        }
                    }
                }

                transitive_closure(component_g, tc);
            }

            // iterate over edges and add scores
            std::vector<EdgeDescriptor> incompatible_edges;
            auto es = boost::edges(tc);
            for (auto eit = es.first; eit != es.second; ++eit)
            {
                // get sink and source and convert to ORF_IDs
                const auto& source = vertex_mapping.at(component_vertex_mapping.at(boost::source(*eit, tc)));
                const auto& target = vertex_mapping.at(component_vertex_mapping.at(boost::target(*eit, tc)));

                float edge_score = std::get<4>(ORF_map.at(target));

                if (overlap_map.find(target) == overlap_map.end())
                {
                    boost::put(boost::edge_weight_t(), tc, *eit, -edge_score);
                    continue;
                } else if (overlap_map.at(target).find(source) == overlap_map.at(target).end())
                {
                    boost::put(boost::edge_weight_t(), tc, *eit, -edge_score);
                    continue;
                }

                // parse the overlap type
                const auto& overlap_pair = overlap_map.at(target).at(source);
                const auto& overlap_type = overlap_pair.first;
                const auto& abs_overlap = overlap_pair.second;

                float penalty = 0;
                if (overlap_type == 'i')
                {
                    incompatible_edges.push_back(*eit);
                    continue;
                } else if (overlap_type == 'u')
                {
                    penalty = unidirectional_penalty_per_base * abs_overlap;
                } else if (overlap_type == 'c')
                {
                    penalty = convergent_penalty_per_base * abs_overlap;
                } else if (overlap_type == 'd')
                {
                    penalty = divergent_penalty_per_base * abs_overlap;
                }

                // negate the penalty from edge_score
                edge_score -= penalty;

                // add edge score, making negative for bellman-ford
                boost::put(boost::edge_weight_t(), tc, *eit, -edge_score);
            }

            // remove incompatible edges
            for (const auto& e : incompatible_edges)
            {
                remove_edge(e, tc);
            }

            // get weight map for graph
            const auto weight_pmap = get(boost::edge_weight_t(), tc);
            
            auto path = traverse_components(ORF_map, vertex_mapping, component_vertex_mapping, component, tc, minimum_path_score, numVertices, weight_pmap, head_kmer_arr);
            // ensure path is not empty
            if (!path.empty())
            {
                // need to convert to 
                gene_paths.push_back(std::move(path));
            }
        }
    }

    return gene_paths;
}
