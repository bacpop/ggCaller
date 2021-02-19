// ggCaller header
#include "ggCaller_classes.h"

// generate ORFs from paths
ORFNodeMap generate_ORFs(const ColoredCDBG<>& ccdbg,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const std::vector<std::pair<std::string, bool>>& unitig_path,
                         const int& overlap,
                         const size_t& min_len)
{
    // initialise path sequence and ORF list
    std::string path_sequence;
    std::vector<std::string> ORF_list;

    // initialise if nodes belong to positive or negative strand in current orientation
    bool positive = true;

    // here generate vector/map which contains the start/end of each node in basepairs. This will then be used to determine which nodes an ORF traverses
    std::vector<std::string> nodelist;
    std::vector<std::vector<size_t>> node_ranges;
    size_t node_start = 0;

    // generate path sequence by merging nodes in sequence
    for (const auto& node : unitig_path)
    {
        const Kmer kmer = Kmer(node.first.c_str());
        auto unitig = ccdbg.find(kmer, true);
        unitig.strand = node.second;

        // add node to node list for path coordinates
        nodelist.push_back(node.first);
        // 0th entry is start index of node within the unitig, 1st entry is end index of node within unitig, 2nd entry is node length. 3rd is strandedness (1 for forward, 0 for reverse)
        std::vector<size_t> node_range(4);

        // add unitig strand to node_range
        node_range[3] = unitig.strand;

        // if strand is negative, calculate reverse complement
        std::string unitig_seq;
        if (unitig.strand)
        {
            unitig_seq = unitig.referenceUnitigToString();
        } else {
            unitig_seq = reverse_complement(unitig.referenceUnitigToString());
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
        node_ranges.push_back(node_range);
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

    std::vector<size_t> frame1;
    std::vector<size_t> frame2;
    std::vector<size_t> frame3;

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
        if (!stop_codon_dict[modulus].empty())
        {
            // initialise first stop codon
            size_t last_index = stop_codon_dict[modulus][0];

            // cycle through start codons
            for (const auto& start_index : start_codon_dict[modulus])
            {
                // check that start index is in correct frame and further in sequence than last stop codon index
                if (start_index > last_index)
                {
                    // cycle through stop codons
                    for (const auto& stop_index : stop_codon_dict[modulus])
                    {
                        // check that stop index is in correct frame and further in sequence than last start codon index
                        if (stop_index > start_index)
                        {
                            // stop index is indexed at first base, therefore end of ORF is two bases after
                            ORF_index_pairs.push_back(std::make_pair(start_index, stop_index + 2));
                            last_index = start_index;
                            break;
                        }
                    }
                }
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
            std::vector<std::string> ORF_node_id;
            std::vector<indexTriplet> ORF_node_coords;
            std::vector<bool> ORF_node_strand;

            // pull 16bp upstream of start codon for TIS model if possible
            if (codon_pair.first >= 16) {
                for (size_t i = 0; i < nodelist.size(); i++)
                {
                    size_t traversed_node_start;
                    size_t traversed_node_end;
                    bool start_assigned = false;
                    bool end_assigned = false;

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

                // generate ORF string
                std::string ORF = path_sequence.substr((codon_pair.first - 16), (ORF_len + 16));

                //std::string ORF = path_sequence.substr((codon_pair.first), (ORF_len));
                ORF_map.emplace(std::move(ORF), std::move(ORF_node_vector));
            }
        }
    }
    return ORF_map;
}

std::pair<ORFColoursMap, std::vector<std::string>> filter_artificial_ORFS(SeqORFMap& all_ORFs,
                                                                          ORFNodeMap& ORF_node_paths,
                                                                          const std::vector<std::string>& fasta_files,
                                                                          const bool write_index)
{
    // call strings edits all_ORFs in place
    // run call strings for all ORFs
    call_strings(all_ORFs, ORF_node_paths, fasta_files, write_index);

    // generate a colours dictionary for gene overlap analysis
    auto return_tuple = sort_ORF_colours(all_ORFs);
    return return_tuple;
}

std::pair<ORFColoursMap, std::vector<std::string>> sort_ORF_colours(const SeqORFMap& all_ORFs)
{
    ORFColoursMap ORF_colours_map;
    std::unordered_set<std::string> ORF_colours_set;

    // generate string for colours
    for (const auto& ORF_seq : all_ORFs)
    {
        std::vector<std::string> colours_vector;
        for (int i = 0; i < ORF_seq.second.size(); i++)
        {
            std::string colours(ORF_seq.second.size(), '0');
            if (ORF_seq.second[i])
            {
                colours[i] = '1';
                colours_vector.push_back(colours);
            }
        }
        for (const auto& colours : colours_vector)
        {
            ORF_colours_map[colours].push_back(ORF_seq.first);
            ORF_colours_set.insert(colours);
        }
    }

    // convert set to vector to enable parallelisation
    std::vector<std::string> ORF_colours_vector(ORF_colours_set.begin(), ORF_colours_set.end());

    // generate return tuple
    std::pair<ORFColoursMap, std::vector<std::string>> return_pair = std::make_pair(ORF_colours_map, ORF_colours_vector);
    return return_pair;
}

std::tuple<SeqORFMap, ORFNodeMap, std::unordered_map<std::string, NodeStrandMap>> call_ORFs(const ColoredCDBG<>& ccdbg,
                                                                                            const PathTuple& path_tuple,
                                                                                            const std::vector<std::string>& stop_codons_for,
                                                                                            const std::vector<std::string>& start_codons_for,
                                                                                            const int& overlap,
                                                                                            const size_t& min_ORF_length)
{
    //initialise all_ORFs to return
    SeqORFMap all_ORFs;
    ORFNodeMap ORF_node_paths;

    #pragma omp parallel
    {
        SeqORFMap all_ORFs_private;
        ORFNodeMap ORF_node_paths_private;
        //std::unordered_map<std::string, std::vector<NodeStrandMap>> pos_strand_map_private;

        // iterate over head_kmer_strings
        #pragma omp for nowait
        for (auto it = std::get<1>(path_tuple).begin(); it < std::get<1>(path_tuple).end(); it++)
        {
            const auto unitig_paths = (std::get<0>(path_tuple)).at(*it);
            // iterate over paths following head_kmer
            for (const auto& path : unitig_paths)
            {
                // CALL ORFS
                const std::vector<bool> path_colour = path.second;
                // iterate over each start codon, generate all ORFs within the path
                for (const auto& start_codon : start_codons_for)
                {
                    std::vector<std::string> start_codon_vector{start_codon};
                    auto ORF_map = generate_ORFs(ccdbg, stop_codons_for, start_codon_vector, path.first, overlap, min_ORF_length);

                    // check if item in all_ORFs already. If not, add colours array. If yes, update the colours array.
                    for (const auto& ORF : ORF_map)
                    {
                        ORF_node_paths_private[ORF.first] = ORF.second;

                        if (all_ORFs_private.find(ORF.first) == all_ORFs_private.end())
                        {
                            all_ORFs_private[ORF.first] = path_colour;
                        } else {
                            std::vector<bool> updated_colours = add_colours_array(all_ORFs_private[ORF.first], path_colour);
                            all_ORFs_private[ORF.first] = updated_colours;
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            // Update ORF_node_paths and pos_strand_vector
            ORF_node_paths.insert(ORF_node_paths_private.begin(), ORF_node_paths_private.end());

            // go through all private ORFs, update colours as before in all_ORFs
            for (const auto& ORF_seq : all_ORFs_private)
            {
                if (all_ORFs.find(ORF_seq.first) == all_ORFs.end())
                {
                    all_ORFs[ORF_seq.first] = ORF_seq.second;
                } else {
                    std::vector<bool> updated_colours = add_colours_array(all_ORFs[ORF_seq.first], ORF_seq.second);
                    all_ORFs[ORF_seq.first] = updated_colours;
                }
            }
        }
    }

    // initialise map to store orientation of nodes
    std::unordered_map<std::string, std::vector<NodeStrandMap>> pos_strand_vector_map;

    // move node strand calculation here, only look at nodes found in ORFs and merge the colours
    for (const auto& ORF_nodes : ORF_node_paths)
    {
        // unpack tuple
        const auto& nodes = std::get<0>(ORF_nodes.second);
        const auto& node_strands = std::get<2>(ORF_nodes.second);

        // generate string for colours
        std::vector<std::string> colours_vector;
        for (int i = 0; i < all_ORFs.at(ORF_nodes.first).size(); i++)
        {
            std::string colours(all_ORFs.at(ORF_nodes.first).size(), '0');
            if (all_ORFs.at(ORF_nodes.first)[i])
            {
                colours[i] = '1';
                colours_vector.push_back(colours);
            }
        }

        // iterate over the single colours parsed from the ORF_node colours
        for (const auto& colours : colours_vector)
        {
            // create entry if colour not present in pos_strand_map_private
            if (pos_strand_vector_map.find(colours) == pos_strand_vector_map.end())
            {
                std::vector<NodeStrandMap> empty_vector;
                pos_strand_vector_map[colours] = empty_vector;
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
                pos_strand_vector.push_back(new_map);
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

    std::unordered_map<std::string, NodeStrandMap> pos_strand_map;

    // create flat map for return
    for (const auto& colours : pos_strand_vector_map)
    {
        pos_strand_map[colours.first] = colours.second[0];
    }

    const auto ORF_tuple = std::make_tuple(all_ORFs, ORF_node_paths, pos_strand_map);
    return ORF_tuple;
}
