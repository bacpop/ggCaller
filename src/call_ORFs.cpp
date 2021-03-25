// ggCaller header
#include "ggCaller_classes.h"

// generate ORFs from paths
ORFNodeMap generate_ORFs(const unitigMap& graph_map,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const std::vector<int>& unitig_path,
                         const int& overlap,
                         const size_t min_len,
                         const bool is_ref,
                         const std::vector<fm_index_coll>& seq_idx,
                         const std::vector<bool>& path_colours,
                         const size_t& nb_colours)
{
    // ORF dictionary for locating ORF position in graph. Nodes traversed by ORF is ordered map to ensure ordering is kept for overlap comparisons.
    ORFNodeMap ORF_map;

    // initialise path sequence
    std::string path_sequence;

    // here generate vector/map which contains the start/end of each node in basepairs. This will then be used to determine which nodes an ORF traverses
    std::vector<int> nodelist;
    std::vector<std::vector<size_t>> node_ranges;
    size_t node_start = 0;

    // testing
//    std::set<std::pair<std::string, bool>> node_set;
//    std::vector<std::pair<std::string, bool>> node_vector;

    // generate path sequence by merging nodes in sequence
    for (const auto& node : unitig_path)
    {
        // parse out information from node integer value
        bool strand = (node >= 0) ? true : false;

        // add node to node list for path coordinates
        nodelist.push_back(node);
        // 0th entry is start index of node within the unitig, 1st entry is end index of node within unitig, 2nd is node end, 3rd entry is node length. 4th is strandedness (1 for forward, 0 for reverse), not used!
        std::vector<size_t> node_range(3);

        // add unitig strand to node_range
        //node_range[3] = strand;

        //testing
//        std::string node_str = graph_map.at(node.first).head_kmer;
//        std::pair<std::string, bool> node_pair(node_str, node.second);
//        node_set.insert(node_pair);
//        node_vector.push_back(node_pair);

        // if strand is negative, calculate reverse complement
        std::string unitig_seq;
        if (strand)
        {
            unitig_seq = graph_map.at(abs(node)).unitig_seq;
        } else {
            unitig_seq = reverse_complement(graph_map.at(abs(node)).unitig_seq);
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

    // testing
//    int test = 0;
//    //std::pair<std::string, bool> head("GGCATGTTGAAAGCGAAGCCATGAACCAGGC", true);
//    std::pair<std::string, bool> query1("GCTTGATATCCATTTCATCCCGCCCGGTCAT", false);
//    std::pair<std::string, bool> query2("TCATCTTATGCCTTTTTTATCTTCATTTGGG", false);
//    if (node_set.find(query1) != node_set.end() && node_set.find(query2) != node_set.end() && path_colours[45])
//    {
//        test = 1;
//    }


    // generate codon indexes using FindIndex function for start and stop codons
    std::vector<size_t> start_codon_indices;
    std::vector<size_t> stop_codon_indices;
    std::vector<std::size_t> found_indices;

    // iterate over each stop codon
    for (const auto& codon : stop_codons)
    {
        found_indices = findIndex(path_sequence, codon, 0, 0, false);
        stop_codon_indices.insert(stop_codon_indices.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();
    }

    // iterate over each start codon
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

    // generate sequences for ORFs from codon pairs
    for (const auto& codon_pair : ORF_index_pairs)
    {

        // add one as codon_pair is zero-indexed
        size_t ORF_len = (codon_pair.second - codon_pair.first) + 1;

        // check ORF is longer than minimum length
        if (ORF_len >= min_len)
        {
            // initialise TIS_present and ORF_colours variables
            auto ORF_colours = path_colours;
            bool TIS_present = true;

            // segmentated to avoid storing ORF seqs longer than required
            {
                // generate ORF sequence. If not present, then ignore
                std::string ORF_seq;

                // get ORF seqeunce and pull 16bp upstream of start codon for TIS model if possible. If not, do not and set TIS_present as false
                if (codon_pair.first >= 16)
                {
                    ORF_seq = path_sequence.substr((codon_pair.first - 16), (ORF_len + 16));
                } else
                {
                    ORF_seq = path_sequence.substr((codon_pair.first), (ORF_len));
                    TIS_present = false;
                }

                // check that ORF present in stated colours
                // check path sequence is real if is_ref
                if (is_ref)
                {
                    // check the colours of the ORF
                    check_colours(ORF_colours, ORF_seq, seq_idx, nb_colours);

                    // accumulate path_colours
                    int sum_colours = accumulate(ORF_colours.begin(), ORF_colours.end(), 0);

                    // check if sum_colours is 0, if so, pass on the ORF
                    if (sum_colours == 0)
                    {
                        continue;
                    }
                }
            }

            // work out coordinates for ORF in node space
            auto ORF_coords = std::move(calculate_coords(codon_pair, nodelist, node_ranges, overlap));

            // work coordinates for TIS in node space
            std::tuple<std::string, std::vector<int>, std::vector<indexTriplet>> TIS_coords;
            if (TIS_present)
            {
                std::pair<std::size_t, std::size_t> TIS_pair(codon_pair.first - 16, codon_pair.first - 1);
                TIS_coords = std::move(calculate_coords(TIS_pair, nodelist, node_ranges, overlap));
            }

            // unpack tuples
            auto& ORF_path_ID = std::get<0>(ORF_coords);
            auto& ORF_node_id = std::get<1>(ORF_coords);
            auto& ORF_node_coords = std::get<2>(ORF_coords);
            //auto& ORF_node_strand = std::get<3>(ORF_coords);

            auto& TIS_path_id = std::get<0>(TIS_coords);
            auto& TIS_node_id = std::get<1>(TIS_coords);
            auto& TIS_node_coords = std::get<2>(TIS_coords);
            //auto& TIS_node_strand = std::get<3>(TIS_coords);

            // add TIS_path_id to ORF_path_ID to ensure ORFs with same sequence but different TIS are not merged
            ORF_path_ID += TIS_path_id;

            // create ORF_node_vector, populate with results from node traversal
            ORFNodeVector ORF_node_vector = std::make_tuple(ORF_node_id, ORF_node_coords, ORF_len, TIS_node_id, TIS_node_coords);

            // add create colours/ORF_node_vector return object
            std::pair<std::vector<bool>, ORFNodeVector> colour_path_pair(ORF_colours, std::move(ORF_node_vector));

            ORF_map.emplace(std::move(ORF_path_ID), std::move(colour_path_pair));
        }
    }
    return ORF_map;
}

// calculate node coordinates for an ORF
std::tuple<std::string, std::vector<int>, std::vector<indexTriplet>> calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                                                                                        const std::vector<int>& nodelist,
                                                                                        const std::vector<std::vector<size_t>>& node_ranges,
                                                                                        const int& overlap)
{
    // initialise items for a tuple containing a vector of each node name, corresponding vector of positions traversed in the node and node strand
    std::vector<int> ORF_node_id;
    std::vector<indexTriplet> ORF_node_coords;
    //std::vector<bool> ORF_node_strand;

    // generate a unique string for the ORF for storage using nodes traversed (will match across matching ORFs)
    std::string ORF_path_ID;

    for (size_t i = 0; i < nodelist.size(); i++)
    {
        size_t traversed_node_start;
        size_t traversed_node_end;
        bool start_assigned = false;
        bool end_assigned = false;

        // if start of ORF is below node range, then check ORF traverses node
        if (codon_pair.first < node_ranges[i][0])
        {
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
        if (codon_pair.second >= node_ranges[i][1])
        {
            traversed_node_end = node_ranges[i][2];
            end_assigned = true;
        } else if (codon_pair.second >= node_ranges[i][0] && codon_pair.second < node_ranges[i][1])
        {
            traversed_node_end = codon_pair.second - node_ranges[i][0];
            // check that the end is greater than the overlap. If not, then sequence is already covered in prior node traversal
            if (traversed_node_end >= overlap)
            {
                end_assigned = true;
            }
        }

        // if the ORF traverses node, update coordinates
        if (start_assigned && end_assigned)
        {
            size_t node_end = node_ranges[i][2];
            indexTriplet node_coords = std::make_tuple(traversed_node_start, traversed_node_end, node_end);
            ORF_node_id.push_back(nodelist[i]);
            ORF_node_coords.push_back(std::move(node_coords));
            //ORF_node_strand.push_back(node_ranges[i][3]);

            // add to ORF_path_ID with node ID and strand. If the ORF_path_ID is empty, add the start position to make the ORF unique
            if (ORF_path_ID.empty())
            {
                ORF_path_ID += std::to_string(std::get<0>(node_coords));
            }
            ORF_path_ID += std::to_string(nodelist[i]);
        }
    }

    const auto return_tuple = std::make_tuple(ORF_path_ID, ORF_node_id, ORF_node_coords);
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
        // generate new ORFIDMap entry, taking ORF PathVector only
        ORF_ID_map[ORF_ID] = std::move(ORF.second.second);

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
            for (const auto& n : nodes)
            {
                // assigned new_map entry. If node is positive, strand is true, if negative, strand is false
                new_map[abs(n)] = (n >= 0) ? true : false;
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



std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>, ColourNodeStrandMap> call_ORFs(const PathPair& path_pair,
                                                                                            const unitigMap& graph_map,
                                                                                             const std::vector<std::string>& stop_codons_for,
                                                                                             const std::vector<std::string>& start_codons_for,
                                                                                             const int overlap,
                                                                                             const size_t min_ORF_length,
                                                                                             const bool is_ref,
                                                                                             const std::vector<fm_index_coll>& seq_idx,
                                                                                             const size_t& nb_colours)
{
    //initialise ORF_nodes_paths to add ORF sequences to
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

                // generate copy of path colours to check for artificial sequences
                auto path_colours = path.second;

                // generate all ORFs within the path for start and stop codon pairs
                const auto ORF_map = std::move(generate_ORFs(graph_map, stop_codons_for, start_codons_for, path.first, overlap, min_ORF_length, is_ref, seq_idx, path_colours, nb_colours));

                // check if item in ORF_node_maps already. If not, add colours array. If yes, update the colours array.
                for (const auto& ORF : ORF_map)
                {
                    // check if ORF is already in ORF_node_paths_private
                    if (ORF_node_paths_private.find(ORF.first) == ORF_node_paths_private.end())
                    {
                        ORF_node_paths_private[ORF.first] = ORF.second;
                    } else {
                        // if not, check if colours need updating
                        const auto& ORF_colours = ORF.second.first;
                        if (ORF_node_paths_private[ORF.first].first != ORF_colours)
                        {
                            std::vector<bool> updated_colours = std::move(add_colours_array(ORF_node_paths_private[ORF.first].first, ORF_colours));
                            ORF_node_paths_private[ORF.first].first = std::move(updated_colours);
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
                    const auto& ORF_colours = ORF.second.first;
                    if (ORF_node_paths[ORF.first].first != ORF_colours)
                    {
                        std::vector<bool> updated_colours = std::move(add_colours_array(ORF_node_paths[ORF.first].first, ORF_colours));
                        ORF_node_paths[ORF.first].first = std::move(updated_colours);
                    }
                }
            }
        }
    }

    // generate pos_strand_map to determine relative strands of each node for each colour
    auto pos_strand_map = std::move(calculate_pos_strand(ORF_node_paths));

    // group colours of ORFs together
    std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>> ORF_colours_tuple = std::move(sort_ORF_colours(ORF_node_paths));

    const auto ORF_tuple = std::make_tuple(std::get<0>(ORF_colours_tuple), std::get<1>(ORF_colours_tuple), std::get<2>(ORF_colours_tuple), pos_strand_map);
    return ORF_tuple;
}
