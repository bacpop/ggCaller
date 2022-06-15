#include "gene_refinding.h"

std::pair<std::vector<std::vector<size_t>>, std::string> calculate_node_ranges(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                               const std::vector<Kmer>& head_kmer_arr,
                                                                               const int& overlap,
                                                                               const std::vector<int>& full_nodelist)
{
    // initialise node ranges to return
    std::vector<std::vector<size_t>> node_ranges;
    size_t node_start = 0;

    // generate path sequence by merging nodes in sequence
    std::string path_sequence;
    for (int i = 0; i < full_nodelist.size(); i++) {
        // parse out information from node integer value
        const auto& node = full_nodelist.at(i);

        // 0th entry is start index of node within the unitig, 1st entry is end index of node within unitig, 2nd is node end
        std::vector<size_t> node_range(3);

        // calculate length of unitig and get end coordinates
        // get a reference to the unitig map object
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, node);
        auto& um = um_pair.first;
        auto& um_data = um_pair.second;

        size_t node_end = um.size;

        std::string seq;

        if (node >= 0)
        {
            seq = um.referenceUnitigToString();
        } else
        {
            seq = reverse_complement(um.referenceUnitigToString());
        }

        if (i == 0) {
            // start index of node in path
            node_range[0] = node_start;
            // start index of next node in path (one past the end index)
            node_range[1] = node_end;
            // relative end of node i.e. node length, zero indexed
            node_range[2] = node_end - 1;
            node_start = node_end;
            path_sequence = seq;
        } else {
            // start index of node in path
            node_range[0] = node_start - overlap;
            // relative end of node i.e. node length, zero indexed
            node_range[2] = node_end - 1;
            node_end += node_start - overlap;
            // start index of next node in path (one past the end index)
            node_range[1] = node_end;
            node_start = node_end;

            // add to path sequencing
            path_sequence.append(seq.begin() + overlap, seq.end());
        }
        node_ranges.push_back(std::move(node_range));
    }


    return {node_ranges, path_sequence};
}

std::pair<std::vector<int>, bool> assign_seq(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                             const std::vector<Kmer>& head_kmer_arr,
                                             RefindPathVector& unitig_complete_paths,
                                             const int kmer,
                                             const bool is_ref,
                                             const fm_index_coll& fm_idx,
                                             const std::vector<int>& curr_nodelist,
                                             const bool downstream)
{
    // initialise path of nodes to return
    std::vector<int> nodelist;
    bool rev_comp = false;

    // sort the unitig_complete_paths by size
    std::stable_sort(unitig_complete_paths.rbegin(), unitig_complete_paths.rend());

    // hold size of previous path
    size_t longest_path = 0;

    // iterate over all the paths, determine which is the longest and real.
    for (const auto& unitig_path_pair : unitig_complete_paths)
    {
        // if shorter, break
        if (unitig_path_pair.first < longest_path)
        {
            break;
        }

        auto unitig_path = unitig_path_pair.second;

        // append the current nodelist to unitig_path, missing off first node
        auto temp_unitig_path = curr_nodelist;
        temp_unitig_path.insert(end(temp_unitig_path), begin(unitig_path) + 1, end(unitig_path));
        unitig_path = std::move(temp_unitig_path);

        // check if path is real, if not then pass
        std::pair<bool, bool> present;
        if (is_ref)
        {
            present = path_search(unitig_path, fm_idx);
            if (!present.first)
            {
                continue;
            }
            rev_comp = present.second;
        }

        // set longest path as current unitig
        longest_path = unitig_path_pair.first;

        nodelist = unitig_path;
    }

    // if not downstream, take reverse of rev_comp
    if (!downstream)
    {
        rev_comp = !rev_comp;
    }

    return {nodelist, rev_comp};
}

RefindPathVector iter_nodes_length (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                  const std::vector<Kmer>& head_kmer_arr,
                                  const NodeTuple& head_node_tuple,
                                  const size_t& current_colour,
                                  const size_t& radius,
                                  const bool& repeat,
                                  const bool& is_ref,
                                  const fm_index_coll& fm_idx,
                                  const int overlap)
{
    // generate path list, vector for path and the stack
    RefindPathVector path_list;
    std::vector<int> node_vector;
    NodeStack node_stack;

    // create first item in stack
    node_stack.push(head_node_tuple);

    while(!node_stack.empty())
    {
        // pop node in stack
        auto node_tuple = node_stack.top();
        node_stack.pop();

        // unpack tuple
        const size_t & pos_idx = std::get<0>(node_tuple);
        const int & node_id = std::get<1>(node_tuple);
        const boost::dynamic_bitset<> & colour_arr = std::get<3>(node_tuple);
        const size_t & path_length = std::get<4>(node_tuple);

        // add node to path
        node_vector = std::vector<int> (node_vector.begin(), node_vector.begin() + pos_idx);
        node_vector.push_back(node_id);

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        // get unitig_dict entry in graph_vector
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        auto& um = um_pair.first;
        auto& um_data = um_pair.second;

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
            auto neighbour_da = neighbour_um.getData();
            auto neighbour_um_data = neighbour_da->getData(neighbour_um);
            const bool neighbour_strand = neighbour_um.strand;

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            auto neighbour_colour = neighbour_um_data->full_colour();
            updated_colours_arr &= neighbour_colour;

            // determine if neighbour is in same colour as iteration, if not return and pass
            if (!(bool)updated_colours_arr[current_colour])
            {
                // update path_list to enable to return all nodes up to point where stop encountered
                path_list.push_back({path_length, node_vector});
                continue;
            } else
            {
                // if not is_ref, check that unitig is shared in at least one other colour
                if (!is_ref)
                {
                    if (neighbour_colour.count() < 2)
                    {
                        // update path_list to enable to return all nodes up to point where stop encountered
                        path_list.push_back({path_length, node_vector});
                        continue;
                    }
                }
            }

            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const int neighbour_id = (neighbour_strand) ? neighbour_um_data->get_id() : neighbour_um_data->get_id() * -1;

            // check against fm-idx every node, pass if not present
            if (is_ref)
            {
                std::vector<int> check_vector = node_vector;
                check_vector.push_back(neighbour_id);
                std::pair<bool, bool> present = path_search(check_vector, fm_idx);
                if (!present.first)
                {
                    // update path_list to enable to return all nodes up to point where stop encountered
                    path_list.push_back({path_length, node_vector});
                    continue;
                }
            } else if (!is_ref && !repeat)
            {
                // if using reads, check if unitig has already been traversed, and pass if repeat not specified
                const bool is_in = std::find(node_vector.begin(), node_vector.end(), neighbour_id) != node_vector.end();
                if (is_in)
                {
                    // update path_list to enable to return all nodes up to point where stop encountered
                    path_list.push_back({path_length, node_vector});
                    continue;
                }
            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + (neighbour_um.size - overlap);

            // if adding new node pushes path length over radius or neighbour is end contig, add neighbour and return
            if (updated_path_length >= radius || neighbour_um_data->end_contig(current_colour))
            {
                // create temporary path to account for reaching end of contig
                std::vector<int> return_path = node_vector;
                return_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::make_pair(updated_path_length, return_path));

                // if node is end_contig, continue to add to path
                if (neighbour_um_data->end_contig(current_colour) && !(updated_path_length >= radius))
                {
                    // if no previous conditions are satisfied, prepare tuple for stack
                    NodeTuple new_node_tuple(new_pos_idx, neighbour_id, 0, updated_colours_arr, updated_path_length);

                    // add to stack
                    node_stack.push(new_node_tuple);
                }
                // move onto next neighbour
                continue;
            }

            // if no previous conditions are satisfied, prepare tuple for stack
            NodeTuple new_node_tuple(new_pos_idx, neighbour_id, 0, updated_colours_arr, updated_path_length);

            // add to stack
            node_stack.push(new_node_tuple);
        }
    }
    return path_list;
}

RefindTuple traverse_outward(const ColoredCDBG<MyUnitigMap>& ccdbg,
                             const std::vector<Kmer>& head_kmer_arr,
                             const size_t& colour_ID,
                             const std::pair<std::vector<int>, std::vector<indexPair>>& ORF_info,
                             const size_t& radius,
                             const bool is_ref,
                             const int kmer,
                             const fm_index_coll& fm_idx,
                             const bool repeat)
{
    // initialise upstream and downstream strings
    std::string upstream_seq;
    std::string downstream_seq;

    // initialise path of nodes to use for coordinate generation and strand of path
    std::vector<int> full_nodelist;
    bool path_rev_comp = true;

    //bool for determining if upstream or downstream empty
    bool upstream_empty = true;
    bool downstream_empty = true;

    // traverse start node in reverse direction
    {
        // get node_id to traverse from, parse unitig_id from TIS if present.
        const int head_id = std::get<0>(ORF_info).at(0) * -1;

        // get reference to unitig_dict object
        // get unitig data
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, head_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // gather unitig information from graph_vector
        const std::bitset<3> codon_arr;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // calculate where in unitig ORF sits, to determine initial length of path
        size_t path_length;

        path_length = std::get<1>(ORF_info).at(0).first;

        // check if path_length already exceeds radius
        RefindPathVector unitig_complete_paths;
        if (path_length >= radius)
        {
            unitig_complete_paths.push_back({path_length, {head_id}});
        } else
        {
            // generate node tuple for iteration
            NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, path_length);

            // recur paths
            unitig_complete_paths = iter_nodes_length(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, radius, repeat, is_ref, fm_idx, kmer - 1);
        }

        if (!unitig_complete_paths.empty())
        {
            auto curr_path = std::get<0>(ORF_info);

            // Reverse curr_path node vector
            std::reverse(curr_path.begin(), curr_path.end());

            // reverse sign on each ID in ORF2_nodes
            for (auto & node_id : curr_path)
            {
                node_id = node_id * -1;
            }

            auto upstream_pair = std::move(assign_seq(ccdbg, head_kmer_arr, unitig_complete_paths, kmer,
                                                       is_ref, fm_idx, curr_path, false));

            auto& upstream_nodes = upstream_pair.first;

            // reverse upstream_nodelist
            if (!upstream_nodes.empty())
            {
                // set path strand as rev_comp output
                path_rev_comp = upstream_pair.second;

                // Reverse upstream_nodes node vector
                std::reverse(upstream_nodes.begin(), upstream_nodes.end());

                // reverse sign on each ID in ORF2_nodes
                for (auto & node_id : upstream_nodes)
                {
                    node_id = node_id * -1;
                }

                upstream_empty = false;

                // assign to full_nodelist
                full_nodelist = std::move(upstream_nodes);
            }
        }
    }


    // traverse end node in forward direction
    {
        // get the node id to traverse from
        const int& head_id = std::get<0>(ORF_info).back();

        // get reference to unitig_dict object
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, head_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // gather unitig information from graph_vector
        const std::bitset<3> codon_arr;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // calculate where in unitig ORF sits, to determine initial length of path
        size_t path_length = um.size - std::get<1>(ORF_info).back().second;

        // check if path_length already exceeds radius
        RefindPathVector unitig_complete_paths;
        if (path_length >= radius)
        {
            unitig_complete_paths.push_back({path_length, {head_id}});
        } else
        {
            // generate node tuple for iteration
            NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, path_length);

            // recur paths
            unitig_complete_paths = iter_nodes_length(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, radius, repeat, is_ref, fm_idx, kmer - 1);
        }

        if (!unitig_complete_paths.empty())
        {
            auto downstream_pair = std::move(assign_seq(ccdbg, head_kmer_arr, unitig_complete_paths, kmer,
                                                        is_ref, fm_idx, full_nodelist, true));

            auto downstream_nodes = downstream_pair.first;

            if (!downstream_nodes.empty())
            {
                // if downstream nodes has extended path, update the path_rev_comp
                path_rev_comp = downstream_pair.second;

                downstream_empty = false;

                full_nodelist = std::move(downstream_nodes);
            }
        }
    }

    // check if upstream or downstream seq present and append to full gene sequence
    if (!downstream_empty || !upstream_empty)
    {
        // calculate the positions of each node within the path
        const auto node_ranges_pair = calculate_node_ranges(ccdbg, head_kmer_arr, kmer - 1, full_nodelist);

        auto& full_node_ranges = node_ranges_pair.first;
        auto& path_seq = node_ranges_pair.second;

        return {path_seq, full_nodelist, full_node_ranges, path_rev_comp};
    } else
    {
        return {};
    }
}

RefindMap refind_in_nodes(const ColoredCDBG<MyUnitigMap>& ccdbg,
                          const std::vector<Kmer>& head_kmer_arr,
                          const size_t colour_ID,
                          const NodeSearchDict& node_search_dict,
                          const size_t radius,
                          const bool is_ref,
                          const int kmer,
                          const fm_index_coll& fm_idx,
                          const bool repeat)
{
    RefindMap refind_map;

    // iterate over nodes in node_search_dict
    for (const auto node_search : node_search_dict)
    {
        const auto& node = node_search.first;
        const auto& seq_search = node_search.second;

        // iterate over search sequences in node_search
        for (const auto& ORF_info : seq_search.second)
        {
            refind_map[node].push_back(std::move(traverse_outward(ccdbg, head_kmer_arr, colour_ID, ORF_info, radius, is_ref,
                                                         kmer, fm_idx, repeat)));
        }
    }
    return refind_map;
}