#include "gene_refinding.h"

std::vector<std::vector<size_t>> calculate_node_ranges(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                       const std::vector<Kmer>& head_kmer_arr,
                                                       const int& overlap,
                                                       const std::vector<int>& full_nodelist)
{
    // initialise node ranges to return
    std::vector<std::vector<size_t>> node_ranges;
    size_t node_start = 0;

    // generate path sequence by merging nodes in sequence
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

        if (i == 0) {
            // start index of node in path
            node_range[0] = node_start;
            // start index of next node in path (one past the end index)
            node_range[1] = node_end;
            // relative end of node i.e. node length, zero indexed
            node_range[2] = node_end - 1;
            node_start = node_end;
        } else {
            // start index of node in path
            node_range[0] = node_start - overlap;
            // relative end of node i.e. node length, zero indexed
            node_range[2] = node_end - 1;
            node_end += node_start - overlap;
            // start index of next node in path (one past the end index)
            node_range[1] = node_end;
            node_start = node_end;
        }
        node_ranges.push_back(std::move(node_range));
    }

    return node_ranges;
}

std::pair<std::vector<int>, bool> assign_seq(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                             const std::vector<Kmer>& head_kmer_arr,
                                             const PathVector& unitig_complete_paths,
                                             const int kmer,
                                             const bool is_ref,
                                             const fm_index_coll& fm_idx,
                                             std::string& stream_seq,
                                             const size_t& ORF_end,
                                             const std::string& ORF_seq)
{
    // initialise path of nodes to return
    std::vector<int> nodelist;
    bool rev_comp = false;

    // iterate over all the paths, determine which is the longest and real.
    // If multiple, choose from that with the lowest hash for first kmer
    for (const auto& unitig_path : unitig_complete_paths)
    {
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

        // initilise path_sequence
        std::string path_sequence = "";

        // generate the path sequence
        for (int i = 0; i < unitig_path.size(); i++)
        {
            const auto& node = unitig_path.at(i);

            // get sequence of node
            // get a reference to the unitig map object
            auto um_pair = get_um_data(ccdbg, head_kmer_arr, node);
            auto& um = um_pair.first;
            auto& um_data = um_pair.second;

            // reverse sequence if strand is negative
            std::string seq;
            if (node >= 0)
            {
                seq = um.referenceUnitigToString();
            } else
            {
                seq = reverse_complement(um.referenceUnitigToString());
            }

            if (i == 0)
            {
                path_sequence.append(seq.begin() + ORF_end + 1, seq.end());
            } else
            {
                path_sequence.append(seq.begin() + kmer - 1, seq.end());
            }
        }

        // add ORF seq to path_seq
        path_sequence = ORF_seq + path_sequence;

        // compare the path_sequence to current stream_seq
        if (path_sequence.size() > stream_seq.size())
        {
            stream_seq = path_sequence;
            nodelist = unitig_path;
        } else if (path_sequence.size() == stream_seq.size() && path_sequence != stream_seq)
        {
            // if equal size, get the hash of the last kmer in each and assign to highest
            const size_t path_hash = hasher{}(path_sequence.substr(path_sequence.size() - kmer, kmer));
            const size_t stream_hash = hasher{}(stream_seq.substr(stream_seq.size() - kmer, kmer));

            if (path_hash > stream_hash)
            {
                stream_seq = path_sequence;
                nodelist = unitig_path;
            }
        }
    }
    return {nodelist, rev_comp};
}

PathVector iter_nodes_length (const ColoredCDBG<MyUnitigMap>& ccdbg,
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
    PathVector path_list;
    std::vector<int> node_vector;
    NodeStack node_stack;

    // create node set for identification of repeats
    std::unordered_set<int> node_set;
    node_set.insert(std::get<1>(head_node_tuple));

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

        // slice path, unless at first node
        if (pos_idx != 0)
        {
            node_vector = std::vector<int> (node_vector.begin(), node_vector.begin() + pos_idx);
            // check if path is real, if not then pass
            std::pair<bool, bool> present;
            if (is_ref)
            {
                present = path_search(node_vector, fm_idx);
                if (!present.first)
                {
                    continue;
                }
            }
        }

        // add node to path
        node_vector.push_back(node_id);

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        // get unitig_dict entry in graph_vector
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        auto& um = um_pair.first;
        auto& um_data = um_pair.second;

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // initialise bool to determine if valid neigbour found
        bool neighbour_found = false;

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

            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const int neighbour_id = (neighbour_strand) ? neighbour_um_data->get_id() : neighbour_um_data->get_id() * -1;

            // check if unitig has already been traversed, and pass if repeat not specified
            const bool is_in = node_set.find(neighbour_id) != node_set.end();
            if (!repeat && is_in)
            {
                continue;
            }

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            updated_colours_arr &= neighbour_um_data->full_colour();

            // determine if neighbour is in same colour as iteration, if not return and pass
            if (!updated_colours_arr[current_colour])
            {
                continue;
            } else
            {
                neighbour_found = true;
            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + (neighbour_um.size - overlap);

            // if adding new node pushes path length over radius or neighbour is end contig, add neighbour and return
            if (updated_path_length >= radius || neighbour_um_data->end_contig())
            {
                // create temporary path to account for reaching end of contig
                std::vector<int> return_path = node_vector;
                return_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(return_path));
                continue;
            }

            // if no previous conditions are satisfied, prepare tuple for stack
            NodeTuple new_node_tuple(new_pos_idx, neighbour_id, 0, updated_colours_arr, updated_path_length);

            // add to stack
            node_stack.push(new_node_tuple);

            // add node to node_set
            node_set.insert(neighbour_id);
        }

        // if neighbour not found, push back the node_vector as is
        if (!neighbour_found)
        {
            // update path_list to enable to return all nodes up to point where stop encountered
            path_list.push_back(node_vector);
        }
    }
    return path_list;
}

RefindTuple traverse_outward(const ColoredCDBG<MyUnitigMap>& ccdbg,
                             const std::vector<Kmer>& head_kmer_arr,
                             const size_t& colour_ID,
                             const ORFNodeVector& ORF_info,
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

    // generate a string of the ORF to check against FM-index
    std::string ORF_seq = generate_sequence_nm(std::get<3>(ORF_info), std::get<4>(ORF_info), kmer - 1, ccdbg, head_kmer_arr);
    ORF_seq += generate_sequence_nm(std::get<0>(ORF_info), std::get<1>(ORF_info), kmer - 1, ccdbg, head_kmer_arr);

    const auto original_ORF = ORF_seq;

    // traverse start node in reverse direction
    {
        // get node_id to traverse from, parse unitig_id from TIS if present.
        const bool TIS_present = (!std::get<3>(ORF_info).empty()) ? true : false;
        const int head_id = ((TIS_present) ? std::get<3>(ORF_info).at(0) : std::get<0>(ORF_info).at(0)) * -1;

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
        size_t ORF_end;

        if (TIS_present)
        {
            path_length = std::get<4>(ORF_info).at(0).first + 16;

            // reverse the ORF end
            // get absolute last node index (same as unitig length minus 1 as zero indexed)
            auto end_um_pair = get_um_data(ccdbg, head_kmer_arr, std::get<3>(ORF_info).at(0));
            const auto& end_um = end_um_pair.first;
            const size_t node_end = end_um.size - 1;

            // get difference from original start to absolute last node index
            ORF_end = node_end - std::get<4>(ORF_info).at(0).first;
        } else
        {
            path_length = std::get<1>(ORF_info).at(0).first;

            // reverse the ORF end
            // get absolute last node index (same as unitig length minus 1 as zero indexed)
            auto end_um_pair = get_um_data(ccdbg, head_kmer_arr, std::get<0>(ORF_info).at(0));
            const auto& end_um = end_um_pair.first;
            const size_t node_end = end_um.size - 1;

            // get difference from original start to absolute last node index
            ORF_end = node_end - std::get<1>(ORF_info).at(0).first;
        }

        // check if path_length already exceeds radius
        PathVector unitig_complete_paths;
        if (path_length >= radius)
        {
            unitig_complete_paths.push_back({head_id});
        } else
        {
            // generate node tuple for iteration
            NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, path_length);

            // recur paths
            unitig_complete_paths = iter_nodes_length(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, radius, repeat, is_ref, fm_idx, kmer - 1);
        }

        if (!unitig_complete_paths.empty())
        {
            auto upstream_pair = std::move(assign_seq(ccdbg, head_kmer_arr, unitig_complete_paths, kmer,
                                                       is_ref, fm_idx, upstream_seq, ORF_end, reverse_complement(ORF_seq)));

            auto& upstream_nodes = upstream_pair.first;

            // reverse upstream_nodelist
            if (!upstream_seq.empty())
            {
                // set path strand as rev_comp output
                path_rev_comp = upstream_pair.second;

                // Reverse ORF2 node vector
                std::reverse(upstream_nodes.begin(), upstream_nodes.end());

                // reverse sign on each ID in ORF2_nodes
                for (auto & node_id : upstream_nodes)
                {
                    node_id = node_id * -1;
                }
            }

            // assign to full_nodelist
            full_nodelist = std::move(upstream_nodes);
        }
    }

    // assign the upstream_seq to ORF_seq to use for upstream checks
    if (!upstream_seq.empty())
    {
        ORF_seq = reverse_complement(upstream_seq);
    }

    // scope for TIS + ORF node_list
    {
        // merge TIS and ORF node lists
        auto temp_nodelist = std::get<3>(ORF_info);
        if (std::get<0>(ORF_info).size() > 1)
        {
            int i = 0;
            for (; i < std::get<0>(ORF_info).size(); i++)
            {
                if (std::find(temp_nodelist.begin(), temp_nodelist.end(), std::get<0>(ORF_info).at(i)) == temp_nodelist.end())
                {
                    break;
                }
            }
            temp_nodelist.insert(temp_nodelist.end(), std::get<0>(ORF_info).begin() + i, std::get<0>(ORF_info).end());
        }
        if (full_nodelist.empty())
        {
            full_nodelist = temp_nodelist;
        }
        else if (temp_nodelist.size() > 1)
        {
            full_nodelist.insert(full_nodelist.end(), temp_nodelist.begin() + 1, temp_nodelist.end());
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
        size_t ORF_end = std::get<1>(ORF_info).back().second;

        // check if path_length already exceeds radius
        PathVector unitig_complete_paths;
        if (path_length >= radius)
        {
            unitig_complete_paths.push_back({head_id});
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
                                                        is_ref, fm_idx, downstream_seq, ORF_end, ORF_seq));

            auto downstream_nodes = downstream_pair.first;

            if (downstream_nodes.size() > 1)
            {
                // if downstream nodes has extended path, update the path_rev_comp
                path_rev_comp = downstream_pair.second;

                full_nodelist.insert(full_nodelist.end(), downstream_nodes.begin() + 1, downstream_nodes.end());
            }
        }
    }

    // check if upstream or downstream seq present and append to full gene sequence
    if (!upstream_seq.empty() || !downstream_seq.empty())
    {
        // append upstream seq to ORF_seq
        ORF_seq = downstream_seq;

        // calculate the positions of each node within the path
        const auto full_node_ranges = calculate_node_ranges(ccdbg, head_kmer_arr, kmer - 1, full_nodelist);

        return {ORF_seq, full_nodelist, full_node_ranges, path_rev_comp};
    } else
    {
        return {};
    }
}

RefindMap refind_in_nodes(const ColoredCDBG<MyUnitigMap>& ccdbg,
                          const std::vector<Kmer>& head_kmer_arr,
                          const size_t colour_ID,
                          const std::unordered_map<int, std::unordered_map<std::string, ORFNodeVector>>& node_search_dict,
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

        // iterate over search sequences in node_search
        for (const auto& seq_search : node_search.second)
        {
            const auto& search = seq_search.first;
            const auto& ORF_info = seq_search.second;

            refind_map[node][search] = traverse_outward(ccdbg, head_kmer_arr, colour_ID, ORF_info, radius, is_ref,
                                                         kmer, fm_idx, repeat);
        }
    }
    return refind_map;
}