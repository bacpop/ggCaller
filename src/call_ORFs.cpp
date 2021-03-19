// ggCaller header
#include "ggCaller_classes.h"

// generate ORFs from paths
ORFNodeMap generate_ORFs(const unitigMap& graph_map,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const std::vector<std::pair<size_t, bool>>& unitig_path,
                         const int overlap,
                         const size_t min_len,
                         const std::vector<bool>& path_colours)
{
    // initialise path sequence and ORF list
    std::string path_sequence;
    std::vector<std::string> ORF_list;

    // here generate vector/map which contains the start/end of each node in basepairs. This will then be used to determine which nodes an ORF traverses
    std::vector<size_t> nodelist;
    std::vector<std::vector<size_t>> node_ranges;
    size_t node_start = 0;

    // generate path sequence by merging nodes in sequence
    for (const auto& node : unitig_path)
    {
        // add node to node list for path coordinates
        nodelist.push_back(node.first);
        // 0th entry is start index of node within the unitig, 1st entry is end index of node within unitig, 2nd entry is node length. 3rd is strandedness (1 for forward, 0 for reverse)
        std::vector<size_t> node_range(4);

        // add unitig strand to node_range
        node_range[3] = node.second;

        // if strand is negative, calculate reverse complement
        std::string unitig_seq;
        if (node.second)
        {
            unitig_seq = graph_map.at(node.first).unitig_seq;
        } else {
            unitig_seq = reverse_complement(graph_map.at(node.first).unitig_seq);
        }

        // calculate length of unitig and get end coordinates
        size_t node_end = unitig_seq.size();

        if (path_sequence.empty())
        {
            path_sequence = unitig_seq;
            // start index of node in path
            node_range[0] = node_start;
            // start index of next node (one past the end index)
            node_range[1] = node_end;
            // absolute end index of node
            node_range[2] = node_end - 1;
            node_start = node_end;
        } else {
            path_sequence.append(unitig_seq.begin()+overlap, unitig_seq.end());
            // start index of node in path
            node_range[0] = node_start - overlap;
            // absolute end index of node
            node_range[2] = node_end - 1;
            node_end += node_start - overlap;
            // start index of next node in path (one past the end index)
            node_range[1] = node_end;
            node_start = node_end;

        }
        node_ranges.push_back(std::move(node_range));
    }

    // generate codon indexes using FindIndex function for start and stop codons
    std::vector<size_t> start_codon_indices;
    std::vector<size_t> stop_codon_indices;
    std::vector<std::size_t> found_indices;

    for (const auto& codon : stop_codons)
    {
        found_indices = findIndex(path_sequence, codon, 0, 0, false);
        stop_codon_indices.insert(stop_codon_indices.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();
    }
    for (const auto& codon : start_codons)
    {
        found_indices = findIndex(path_sequence, codon, 0, 0, false);
        start_codon_indices.insert(start_codon_indices.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();
    }

    // create pair for each paired start and stop codon
    std::vector<std::pair<std::size_t, std::size_t>> ORF_index_pairs;

    // sort codon index arrays into ascending order
    std::sort(start_codon_indices.begin(), start_codon_indices.end());
    std::sort(stop_codon_indices.begin(), stop_codon_indices.end());

    // generate dictionaries for start and stop codon indices for each frame
    std::unordered_map<size_t, std::vector<size_t>> start_codon_dict;
    std::unordered_map<size_t, std::vector<size_t>> stop_codon_dict;

    // initialise vectors to determine indexes in each frame of path
    std::vector<size_t> frame1;
    std::vector<size_t> frame2;
    std::vector<size_t> frame3;

    // parse indexes of start codons
    for (auto& index : start_codon_indices)
    {
        if (index % 3 == 0 || index == 0)
        {
            frame1.push_back(index);
        } else if (index % 3 == 1 || index == 1)
        {
            frame2.push_back(index);
        } else if (index % 3 == 2 || index == 2)
        {
            frame3.push_back(index);
        }
    }
    // update start codon indexes
    start_codon_dict[0] = std::move(frame1);
    start_codon_dict[1] = std::move(frame2);
    start_codon_dict[2] = std::move(frame3);
    // clear all index vectors
    frame1.clear();
    frame2.clear();
    frame3.clear();
    start_codon_indices.clear();

    // parse indexes of stop codons
    for (auto& index : stop_codon_indices)
    {
        if (index % 3 == 0 || index == 0)
        {
            frame1.push_back(index);
        } else if (index % 3 == 1 || index == 1)
        {
            frame2.push_back(index);
        } else if (index % 3 == 2 || index == 2)
        {
            frame3.push_back(index);
        }
    }
    // update stop codon indexes
    stop_codon_dict[0] = std::move(frame1);
    stop_codon_dict[1] = std::move(frame2);
    stop_codon_dict[2] = std::move(frame3);
    // clear all index vectors
    frame1.clear();
    frame2.clear();
    frame3.clear();
    stop_codon_indices.clear();

    // iterate through frames, pair sequential start+stop codons after first stop codon
    for (int modulus = 0; modulus < 3; modulus++)
    {
        // initialise stop and start index iterators
        auto start_index = start_codon_dict[modulus].begin();
        auto stop_index = stop_codon_dict[modulus].begin();

        // iterate over start and stop indices until you reach the end
        while (stop_index != stop_codon_dict[modulus].end() && start_index != start_codon_dict[modulus].end())
        {
            if (*start_index < *stop_index)
            {
                // stop index is indexed at first base, therefore end of ORF is two bases after
                std::pair<size_t, size_t> codon_pair(*start_index, *stop_index + 2);
                ORF_index_pairs.push_back(std::move(codon_pair));

                // iterate start index
                start_index++;
            } else {
                // start is greater than stop, so iterate stop_iterator
                stop_index++;
            }
        }
    }

    // ORF dictionary for locating ORF position in graph. Nodes traversed by ORF is ordered map to ensure ordering is kept for overlap comparisons.
    ORFNodeMap ORF_map;

    // generate sequences for ORFs from codon pairs
    for (const auto& codon_pair : ORF_index_pairs)
    {
        // add one as codon_pair is zero-indexed
        size_t ORF_len = (codon_pair.second - codon_pair.first) + 1;

        if (ORF_len >= min_len)
        {
            // initialise items for a tuple containing a vector of each node name, corresponding vector of positions traversed in the node and node strand
            std::vector<size_t> ORF_node_id;
            std::vector<indexTriplet> ORF_node_coords;
            std::vector<bool> ORF_node_strand;

            // pull 16bp upstream of start codon for TIS model if possible
            if (codon_pair.first >= 16)
            {
                for (size_t i = 0; i < nodelist.size(); i++)
                {
                    size_t traversed_node_start;
                    size_t traversed_node_end;
                    bool start_assigned = false;
                    bool end_assigned = false;

                    // if start of ORF is below node range, then check ORF traverses node
                    if (codon_pair.first < node_ranges[i][0]){
                        traversed_node_start = 0;
                        start_assigned = true;
                    } else if (codon_pair.first >= node_ranges[i][0] && codon_pair.first < node_ranges[i][1]){
                        traversed_node_start = codon_pair.first - node_ranges[i][0];
                        // check that the start difference to end is greater than the overlap. If not, then sequence is covered in next node traversal
                        if ((node_ranges[i][2] - traversed_node_start) >= overlap)
                        {
                            start_assigned = true;
                        }
                    }

                    // if end of ORF is above node range, then check if ORF traversed
                    if (codon_pair.second >= node_ranges[i][1]) {
                        traversed_node_end = node_ranges[i][2];
                        end_assigned = true;
                    } else if (codon_pair.second >= node_ranges[i][0] && codon_pair.second < node_ranges[i][1]){
                        traversed_node_end = codon_pair.second - node_ranges[i][0];
                        // check that the end is greater than the overlap. If not, then sequence is already covered in prior node traversal
                        if (traversed_node_end >= overlap)
                        {
                            end_assigned = true;
                        }
                    }

                    // if the ORF traverses node, update coordinates
                    if (start_assigned && end_assigned){
                        size_t node_end = node_ranges[i][2];
                        indexTriplet node_coords = std::make_tuple(traversed_node_start, traversed_node_end, node_end);
                        ORF_node_id.push_back(nodelist[i]);
                        ORF_node_coords.push_back(std::move(node_coords));
                        ORF_node_strand.push_back(node_ranges[i][3]);
                    }

                }
                // create ORF_node_vector
                ORFNodeVector ORF_node_vector = std::make_tuple(ORF_node_id, ORF_node_coords, ORF_node_strand);

                // add create colours/ORF_node_vector return object
                std::pair<std::vector<bool>, ORFNodeVector> colour_path_pair(path_colours, std::move(ORF_node_vector));

                // generate ORF string
                std::string ORF = path_sequence.substr((codon_pair.first - 16), (ORF_len + 16));

                //std::string ORF = path_sequence.substr((codon_pair.first), (ORF_len));
                ORF_map.emplace(std::move(ORF), std::move(colour_path_pair));
            }
        }
    }
    return ORF_map;
}

std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>> filter_artificial_ORFS(ORFNodeMap& ORF_node_paths,
                                                                                    const std::vector<std::string>& fasta_files,
                                                                                    const bool write_index)
{
    // run call strings for ORF_node_paths in place
    call_strings(ORF_node_paths, fasta_files, write_index);

    // generate a colours dictionary for gene overlap analysis
    auto return_tuple = sort_ORF_colours(ORF_node_paths);
    return return_tuple;
}

std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>> sort_ORF_colours(ORFNodeMap& ORF_node_paths)
{
    ORFColoursMap ORF_colours_map;
    ORFIDMap ORF_ID_map;
    std::unordered_set<size_t> ORF_colours_set;

    // generate string for colours and IDs for each ORF
    size_t ORF_ID = 0;
    for (const auto& ORF : ORF_node_paths)
    {
        std::vector<size_t> colours_vector;
        const auto& ORF_colours = ORF.second.first;
        for (size_t i = 0; i < ORF_colours.size(); i++)
        {
            if (ORF_colours[i])
            {
                colours_vector.push_back(i);
            }
        }
        // generate new ORFIDMap entry, taking ORF sequence and ORF PathVector
        ORF_ID_map[ORF_ID] = std::make_pair(ORF.first, ORF.second.second);

        for (const auto& colours : colours_vector)
        {
            ORF_colours_map[colours].push_back(ORF_ID);
            ORF_colours_set.insert(colours);
        }

        // iterate ORF_ID
        ORF_ID++;
    }

    // convert set to vector to enable parallelisation
    std::vector<size_t> ORF_colours_vector(ORF_colours_set.begin(), ORF_colours_set.end());

    // clear ORF_node_paths
    ORF_node_paths.clear();

    // generate return tuple
    const auto return_tuple = std::make_tuple(ORF_colours_map, ORF_ID_map, ORF_colours_vector);
    return return_tuple;
}

// calculate the relative strand of each node traversed in an ORF, per colour
ColourNodeStrandMap calculate_pos_strand(const ORFNodeMap& ORF_node_paths)
{
    // initialise map to store orientation of nodes (colour is an ID, not a string
    robin_hood::unordered_map<size_t, std::vector<NodeStrandMap>> pos_strand_vector_map;

    // move node strand calculation here, only look at nodes found in ORFs and merge the colours
    for (const auto& ORF_nodes : ORF_node_paths)
    {
        // unpack tuple, colours is ORF_nodes.second.first, path is ORF_nodes.second.second
        const auto& nodes = std::get<0>(ORF_nodes.second.second);
        const auto& node_strands = std::get<2>(ORF_nodes.second.second);
        const auto& ORF_colours = ORF_nodes.second.first;

        // generate string for colours
        std::vector<std::size_t> colours_vector;
        for (int i = 0; i < ORF_colours.size(); i++)
        {
            if (ORF_colours[i])
            {
                colours_vector.push_back(i);
            }
        }

        // iterate over the single colours parsed from the ORF_node colours
        for (const auto& colours : colours_vector)
        {
            // create entry if colour not present in pos_strand_map_private
            if (pos_strand_vector_map.find(colours) == pos_strand_vector_map.end())
            {
                std::vector<NodeStrandMap> empty_vector;
                pos_strand_vector_map[colours] = std::move(empty_vector);
            }

            // make reference to pos_strand_map_private[colours]
            std::vector<NodeStrandMap>& pos_strand_vector = pos_strand_vector_map[colours];

            // create new map to store node info
            NodeStrandMap new_map;
            for (int node_index = 0; node_index < nodes.size(); node_index++)
            {
                new_map[nodes.at(node_index)] = node_strands.at(node_index);
            }

            // create vector to keep track of which maps need to be added to, and whether this is in positive (true) or negative (false) orientation
            std::vector<std::pair<size_t, bool>> maps_to_add;

            // if pos_strand_vector is empty, add first entry of nodes in orientation as is
            if (pos_strand_vector.empty())
            {
                pos_strand_vector.push_back(std::move(new_map));
            }
            // if not, search for nodes in existing maps
            else
            {
                for (size_t i = 0; i < pos_strand_vector.size(); i++)
                {
                    for (const auto& node : new_map)
                    {
                        // if found, note which maps to add the nodes to
                        if (pos_strand_vector[i].find(node.first) != pos_strand_vector[i].end())
                        {
                            // if the strand is the same across the map and the node, add as is
                            if (pos_strand_vector[i][node.first] == node.second)
                            {
                                std::pair in_map(i, true);
                                maps_to_add.push_back(std::move(in_map));
                            } else{
                                std::pair in_map(i, false);
                                maps_to_add.push_back(std::move(in_map));
                            }
                            break;
                        }
                    }
                }
                // if nodes are not found, then add new map
                if (maps_to_add.empty())
                {
                    pos_strand_vector.push_back(new_map);
                }
                // else, go through and add the new map entries to the existing maps in pos_strand_vector_private
                else {
                    for (const auto& map_index : maps_to_add)
                    {
                        if (map_index.second == true)
                        {
                            // if not present, and node is in strand, add as is
                            pos_strand_vector[map_index.first].insert(new_map.begin(), new_map.end());
                        } else
                        {
                            for (const auto& node : new_map)
                            {
                                // if not present, and node is in opposite strand, add reversed strand
                                if (pos_strand_vector[map_index.first].find(node.first) == pos_strand_vector[map_index.first].end())
                                {
                                    pos_strand_vector[map_index.first][node.first] = !node.second;
                                }
                            }
                        }
                    }
                    // if maps to add is greater than one, then means that the detected maps can be merged
                    if (maps_to_add.size() > 1)
                    {
                        std::vector<size_t> to_remove;
                        for (size_t i = 1; i < maps_to_add.size(); i++)
                        {
                            // if maps to add orientation is the same between the two maps, then can be merged as is. Merge all into first map in maps_to_add
                            if (maps_to_add[0].second == maps_to_add[i].second)
                            {
                                pos_strand_vector[maps_to_add[0].first].insert(pos_strand_vector[maps_to_add[i].first].begin(), pos_strand_vector[maps_to_add[i].first].end());

                            } else
                            {
                                // if not in same orientation, need to reverse the strands of all nodes in the merging vector
                                for (const auto& node : pos_strand_vector[maps_to_add[i].first])
                                {
                                    pos_strand_vector[maps_to_add[0].first][node.first] = !node.second;
                                }
                            }
                            // set up to remove the merged maps from the pos_strand_vector_private vector
                            to_remove.push_back(maps_to_add[i].first);
                        }
                        // remove the merged maps
                        // reverse to_remove so that indexes aren't affected
                        std::reverse(to_remove.begin(), to_remove.end());
                        for (const auto& index : to_remove)
                        {
                            pos_strand_vector.erase(pos_strand_vector.begin() + index);
                        }
                    }
                }
            }
        }
    }

    // iterate over pos_strand_vector, merging any maps that have matching ORFs. If no matching ORFs, merge as is. Continue until only single item present in pos_strand_vector
    for (auto& colours : pos_strand_vector_map)
    {
        auto& pos_strand_vector = colours.second;
        while (pos_strand_vector.size() > 1)
        {
            // merge all vectors with the first map as is as no nodes should be shared
            for (size_t index = 1; index < pos_strand_vector.size(); index++)
            {
                pos_strand_vector[0].insert(pos_strand_vector[index].begin(), pos_strand_vector[index].end());
            }

            // remove all but first entry
            pos_strand_vector.erase(pos_strand_vector.begin() + 1, pos_strand_vector.begin() + pos_strand_vector.size());
        }
    }

    ColourNodeStrandMap pos_strand_map;

    // create flat map for return
    for (const auto& colours : pos_strand_vector_map)
    {
        pos_strand_map[colours.first] = colours.second[0];
    }

    return pos_strand_map;
}



std::pair<ORFNodeMap, ColourNodeStrandMap> call_ORFs(const PathPair& path_pair,
                                                     const unitigMap& graph_map,
                                                     const std::vector<std::string>& stop_codons_for,
                                                     const std::vector<std::string>& start_codons_for,
                                                     const int overlap,
                                                     const size_t min_ORF_length)
{
    //initialise ORF_nodes_paths to return
    ORFNodeMap ORF_node_paths;

    #pragma omp parallel
    {
        ORFNodeMap ORF_node_paths_private;

        // iterate over head_kmer_strings
        #pragma omp for nowait
        for (auto it = path_pair.second.begin(); it < path_pair.second.end(); it++)
        {
            const auto& unitig_paths = (path_pair.first).at(*it);
            // iterate over paths following head_kmer
            for (const auto& path : unitig_paths)
            {
                // CALL ORFS
                const auto& path_colours = path.second;
                // iterate over each start codon, generate all ORFs within the path
                for (const auto& start_codon : start_codons_for)
                {
                    std::vector<std::string> start_codon_vector{start_codon};
                    const auto ORF_map = std::move(generate_ORFs(graph_map, stop_codons_for, start_codon_vector, path.first, overlap, min_ORF_length, path_colours));

                    // check if item in ORF_node_maps already. If not, add colours array. If yes, update the colours array.
                    for (const auto& ORF : ORF_map)
                    {
                        // check if ORF is already in ORF_node_paths_private
                        if (ORF_node_paths_private.find(ORF.first) == ORF_node_paths_private.end())
                        {
                            ORF_node_paths_private[ORF.first] = ORF.second;
                        } else {
                            // if not, check if colours need updating
                            if (ORF_node_paths_private[ORF.first].first != path_colours)
                            {
                                std::vector<bool> updated_colours = std::move(add_colours_array(ORF_node_paths_private[ORF.first].first, path_colours));
                                ORF_node_paths_private[ORF.first].first = std::move(updated_colours);
                            }
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            // Update ORF_node_paths

            // go through all private ORFs, update colours as before in ORF_node_paths
            for (const auto& ORF : ORF_node_paths_private)
            {
                if (ORF_node_paths.find(ORF.first) == ORF_node_paths.end())
                {
                    ORF_node_paths[ORF.first] = ORF.second;
                } else {
                    if (ORF_node_paths[ORF.first].first != ORF.second.first)
                    {
                        std::vector<bool> updated_colours = std::move(add_colours_array(ORF_node_paths[ORF.first].first, ORF.second.first));
                        ORF_node_paths[ORF.first].first = std::move(updated_colours);
                    }
                }
            }
        }
    }

    // generate pos_strand_map to determine relative strands of each node for each colour
    auto pos_strand_map = std::move(calculate_pos_strand(ORF_node_paths));

    const auto ORF_pair = std::make_pair(ORF_node_paths, pos_strand_map);
    return ORF_pair;
}
