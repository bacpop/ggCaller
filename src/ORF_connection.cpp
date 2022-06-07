#include "ORF_connection.h"
robin_hood::unordered_map<size_t, std::vector<size_t>> add_ORF_info (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                     const std::vector<Kmer>& head_kmer_arr,
                                                                     const robin_hood::unordered_set<size_t>& target_ORFs,
                                                                     ORFNodeMap& gene_map,
                                                                     const int& overlap)
{

    robin_hood::unordered_map<size_t, std::vector<std::pair<unsigned int, size_t>>> ORF_coords_in_nodes;

    for (const auto & source : target_ORFs)
    {
        // get graph information for source node
        auto & ORF_info = gene_map.at(source);

        // remove regions in unitig overlaps
        simplify_ORFNodeVector(ORF_info, overlap);

        // get start and stop node of ORF to add to node_to_ORFs
        std::unordered_set<size_t> front_back_index = {0, std::get<0>(ORF_info).size() - 1};

        for (const auto i : front_back_index)
        {
            const auto& node_id = std::get<0>(ORF_info).at(i);
            const auto& node_indices = std::get<1>(ORF_info).at(i);

            auto index_start = node_indices.first;
            auto index_end = node_indices.second;

            auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
            auto& um = um_pair.first;

            size_t node_end = um.size - 1;

            // if negative, get first reverse start ORF
            if (node_id < 0)
            {
                index_start = node_end - index_end;
            }

            // add the ORF ID to node_to_ORFs
            ORF_coords_in_nodes[abs(node_id) - 1].push_back({index_start, source});
        }
    }

    robin_hood::unordered_map<size_t, std::vector<size_t>> node_to_ORFs;

    // iterate over ORF_coords_in_nodes, sorting and then filling node_to_ORFs
    for (auto& entry : ORF_coords_in_nodes)
    {
        // sort based on position in node
        std::stable_sort(entry.second.begin(), entry.second.end());
        for (const auto& node_id_pair : entry.second)
        {
            node_to_ORFs[entry.first].push_back(node_id_pair.second);
        }
    }


    return node_to_ORFs;
}

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
                                                    const fm_index_coll& fm_idx)
{
    // initialise pair of vectors (first = upstream of start_ORF, second = downstream of start_ORF)
    std::set<std::pair<size_t, size_t>> ORF_edges;

    // iterate over each entry in end_ORFs
    for (const auto& source : target_ORFs)
    {
        // get ORF info
        const auto & source_info = gene_map.at(source);

        // get id for source_node and sink node for ORF
        auto source_node_id = std::get<0>(source_info).at(0);
        auto sink_node_id = std::get<0>(source_info).back();

        // decide whether to traverse from source (upstream) or sink (downstream)
        int start_node = (stream > 0) ? sink_node_id : source_node_id;

        // check if ORF has already been pair
        bool prev_paired = false;

        if (stream > 0)
        {
            if (downstream_ORF_set.find(source) != downstream_ORF_set.end())
            {
                prev_paired = true;
            }
        } else
        {
            if (upstream_ORF_set.find(source) != upstream_ORF_set.end())
            {
                prev_paired = true;
            }
        }

        // check if node has been traversed previously in opposite direction i.e. ORF has already been paired up/downstream
        if (!prev_paired)
        {
            // get source_node and sink_node info
            // get a reference to the unitig map object
            auto start_um_pair = get_um_data(ccdbg, head_kmer_arr, start_node);
            auto& start_um = start_um_pair.first;
            auto& start_um_data = start_um_pair.second;

            // check if there are any overlapping ORFs on current node and order...
            auto& start_node_ORFs = node_to_ORFs.at(abs(start_node) - 1);

            // if moving in reverse direction, then reverse start_node_ORFs
            if (start_node_ORFs.size() > 1)
            {
                if (start_node * stream < 0)
                {
                    std::reverse(start_node_ORFs.begin(), start_node_ORFs.end());
                }
                // find index of entry
                auto index = std::find(start_node_ORFs.begin(), start_node_ORFs.end(), source) - start_node_ORFs.begin();

                // if at end of path, then can't add anything
                if (index + 1 != start_node_ORFs.size())
                {
                    bool pair_assigned = false;
                    size_t correct_next_ORF;

                    if (is_ref)
                    {
                        // make vector to search with, one less than ORF_nodes
                        std::vector<int> search_vector;
                        const auto & source_info = gene_map.at(source);
                        const auto & source_nodes = std::get<0>(source_info);

                        if (stream == -1)
                        {
                            for (int i = source_nodes.size() - 1; i > 0; i--)
                            {
                                search_vector.push_back(source_nodes.at(i) * -1);
                            }
                        } else
                        {
                            search_vector.insert(search_vector.end(), source_nodes.begin(), source_nodes.end() - 1);
                        }

                        // check that next ORFs are correct
                        auto search_vector_temp = search_vector;
                        for (int j = index + 1; j < start_node_ORFs.size(); j++)
                        {
                            const auto& ORF = start_node_ORFs.at(j);
                            auto search_vector_temp_iter = search_vector_temp;
                            const auto & ORF_info = gene_map.at(ORF);
                            const auto & ORF_nodes = std::get<0>(ORF_info);

                            if (source_nodes.back() == ORF_nodes.at(0))
                            {
                                search_vector_temp_iter.insert(search_vector_temp_iter.end(), ORF_nodes.begin(), ORF_nodes.end());
                            } else
                            {
                                for (int i = ORF_nodes.size() - 1; i > -1; i--)
                                {
                                    search_vector_temp_iter.push_back(ORF_nodes.at(i) * -1);
                                }
                            }

                            const auto present = path_search(search_vector_temp_iter, fm_idx);
                            if (present.first)
                            {
                                correct_next_ORF = ORF;
                                pair_assigned = true;
                                break;
                            }
                        }
                    } else
                    {
                        correct_next_ORF = start_node_ORFs.at(index + 1);
                        pair_assigned = true;
                    }

                    if (pair_assigned)
                    {
                        if (source <= correct_next_ORF)
                        {
                            ORF_edges.insert({source, correct_next_ORF});
                        } else
                        {
                            ORF_edges.insert({correct_next_ORF, source});
                        }

                        // determine how next ORF is orientated
                        const auto & ORF_info = gene_map.at(correct_next_ORF);
                        const auto & ORF_nodes = std::get<0>(ORF_info);

                        // if match between strands, then must be upstream of current ORF
                        if (start_node * stream == ORF_nodes.at(0))
                        {
                            upstream_ORF_set.insert(correct_next_ORF);
                        } else
                        {
                            downstream_ORF_set.insert(correct_next_ORF);
                        }

                        // assign previously paired
                        if (stream > 0)
                        {
                            downstream_ORF_set.insert(source);
                        } else
                        {
                            upstream_ORF_set.insert(source);
                        }

                        // continue for loop
                        continue;
                    }
                }
            }

            // if going downstream (stream > 0), check that ORF is last, or upstream (stream < 0), check it is first
            auto next_ORFs = check_next_ORFs(ccdbg, head_kmer_arr, node_to_ORFs, start_node, source, colour_ID, stream, gene_map, max_ORF_path_length, repeat, downstream_ORF_set, upstream_ORF_set, overlap, is_ref, fm_idx);

            const auto & source_info = gene_map.at(source);

            ORF_edges.insert(make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
        }
    }
    return ORF_edges;
}

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
                                                     const fm_index_coll& fm_idx)
{
    // initialise return vector of upstream ORFs (vectors will break at branches in graph)
    std::set<std::pair<size_t, size_t>> connected_ORFs;

    // get the colour array of the head node
    auto um_pair = get_um_data(ccdbg, head_kmer_arr, head_node);
    auto& um = um_pair.first;
    auto& um_data = um_pair.second;

    auto colour_arr = um_data->full_colour();

    // generate path list, vector for path and the stack
    ORFStack ORF_stack;
    std::vector<int> node_vector;

    // create node set for identification of repeats
    std::unordered_set<int> node_set;
    node_set.insert(head_node * stream);

    // make vector to search with, one less than ORF_nodes
    std::vector<int> search_vector;
    if (is_ref)
    {
        const auto & source_info = gene_map.at(stream_source);
        const auto & ORF_nodes = std::get<0>(source_info);

        if (stream == -1)
        {
            for (int i = ORF_nodes.size() - 1; i > 0; i--)
            {
                search_vector.push_back(ORF_nodes.at(i) * -1);
            }
        } else
        {
            search_vector.insert(search_vector.end(), ORF_nodes.begin(), ORF_nodes.end() - 1);
        }
    }

    // create first item in stack, multiply by stream - upstream (stream = -1) or downstream (stream = 1) and add stream ORF
    ORF_stack.push(std::make_tuple(0, head_node * stream, colour_arr, 0));

    while(!ORF_stack.empty())
    {
        // pop node in stack
        auto path_tuple = ORF_stack.top();
        ORF_stack.pop();
        const size_t & pos_idx = std::get<0>(path_tuple);
        const auto& node_id = std::get<1>(path_tuple);
        colour_arr = std::get<2>(path_tuple);
        auto& path_length = std::get<3>(path_tuple);

        if (pos_idx != 0)
        {
            node_vector = std::vector<int> (node_vector.begin(), node_vector.begin() + pos_idx);
        }

        // add node to path
        node_vector.push_back(node_id);

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        // get unitig_dict entry in graph_vector
        um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        um = um_pair.first;
        um_data = um_pair.second;

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // reverse the strand of the unitig
        if (!strand)
        {
            um.strand = !um.strand;
        }

        // iterate over neighbours, recurring through incomplete paths
        for (auto& neighbour_um : um.getSuccessors())
        {
            // make copy of path_length
            auto temp_path_length = path_length;

            auto neighbour_da = neighbour_um.getData();
            auto neighbour_um_data = neighbour_da->getData(neighbour_um);
            const bool neighbour_strand = neighbour_um.strand;

            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const int neighbour_id = (neighbour_strand) ? neighbour_um_data->get_id() : neighbour_um_data->get_id() * -1;

            // check against fm-idx every other node, pass if not present
            if (is_ref && node_vector.size() % 2 == 1)
            {
                std::vector<int> check_vector = node_vector;
                check_vector.push_back(neighbour_id);
                std::pair<bool, bool> present = path_search(check_vector, fm_idx);
                if (!present.first)
                {
                    continue;
                }
            } else if (!is_ref && !repeat)
            {
                // if using reads, check if unitig has already been traversed, and pass if repeat not specified
                const bool is_in = std::find(node_vector.begin(), node_vector.end(), neighbour_id) != node_vector.end();
                if (is_in)
                {
                    continue;
                }
            }

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            auto neighbour_colour = neighbour_um_data->full_colour();
            updated_colours_arr &= neighbour_colour;

            // determine if neighbour is in same colour as iteration, if not pass
            if (!(bool)updated_colours_arr[current_colour])
            {
                continue;
            }

            // if not is_ref, check that unitig is shared in at least one other colour
            if (!is_ref)
            {
                if (neighbour_colour.count() < 2)
                {
                    continue;
                }
            }

            // check if node is traversed by end of an ORF
            const auto ORF_found = node_to_ORFs.find(abs(neighbour_id) - 1);
            if (ORF_found != node_to_ORFs.end())
            {
                //pull out next ORFs
                auto next_ORFs = ORF_found->second;

                // if moving in reverse direction, then reverse start_node_ORFs
                if (!neighbour_strand)
                {
                    std::reverse(next_ORFs.begin(), next_ORFs.end());
                }

                // hold correct_next_ORF
                size_t correct_next_ORF;
                bool pair_found = false;

                // iterate over ORFs in node, determine which are correct in terms of genome
                // check full gene paths against fm-index, if not present then pass
                if (is_ref)
                {
                    auto search_vector_temp = search_vector;
                    search_vector_temp.insert(search_vector_temp.end(), node_vector.begin(), node_vector.end());
                    for (const auto & ORF : next_ORFs)
                    {
                        auto search_vector_temp_iter = search_vector_temp;
                        const auto & ORF_info = gene_map.at(ORF);
                        const auto & ORF_nodes = std::get<0>(ORF_info);

                        if (neighbour_id == ORF_nodes.at(0))
                        {
                            search_vector_temp_iter.insert(search_vector_temp_iter.end(), ORF_nodes.begin(), ORF_nodes.end());
                        } else
                        {
                            for (int i = ORF_nodes.size() - 1; i > -1; i--)
                            {
                                search_vector_temp_iter.push_back(ORF_nodes.at(i) * -1);
                            }
                        }

                        const auto present = path_search(search_vector_temp_iter, fm_idx);
                        if (present.first)
                        {
                            correct_next_ORF = ORF;
                            pair_found = true;
                            break;
                        }
                    }
                }

                // if pair found, add to connected ORFs
                if (pair_found)
                {
                    // determine how next ORF is orientated
                    const auto & ORF_info = gene_map.at(correct_next_ORF);
                    const auto & ORF_nodes = std::get<0>(ORF_info);

                    // if match between strands, then must be upstream of current ORF
                    if (neighbour_id == ORF_nodes.at(0))
                    {
                        upstream_ORF_set.insert(correct_next_ORF);
                    } else
                    {
                        downstream_ORF_set.insert(correct_next_ORF);
                    }

                    // assign previously paired
                    if (stream > 0)
                    {
                        downstream_ORF_set.insert(stream_source);
                    } else
                    {
                        upstream_ORF_set.insert(stream_source);
                    }

                    // add first entry and stream_source
                    // add ORFs so that lowest ORF_ID is added first (to remove redundant edges)
                    std::pair<size_t, size_t> ORF_pair;
                    ORF_pair.first = (stream_source <= correct_next_ORF) ? stream_source : correct_next_ORF;
                    ORF_pair.second = (stream_source <= correct_next_ORF) ? correct_next_ORF : stream_source;
                    connected_ORFs.insert(std::move(ORF_pair));
                }
            }
            else
            {
                // add to stack if max_ORF_path_length not violated. If next node violates, still traverse to see if meet ORF
                if (temp_path_length <= max_ORF_path_length)
                {
                    temp_path_length += neighbour_um.size - overlap;
                    ORF_stack.push({new_pos_idx, neighbour_id, updated_colours_arr, temp_path_length});
                }
            }
        }
    }
    return connected_ORFs;
}

