// ggCaller header

#include "call_ORFs.h"

// generate ORFs from paths
void generate_ORFs(const int& colour_ID,
                   ORFNodeMap& ORF_node_map,
                   std::unordered_set<size_t>& hashes_to_remove,
                   const ColoredCDBG<MyUnitigMap>& ccdbg,
                   const std::vector<Kmer>& head_kmer_arr,
                   const std::vector<std::string>& stop_codons,
                   const std::vector<std::string>& start_codons,
                   const std::vector<int>& nodelist,
                   const int& overlap,
                   const size_t min_len,
                   const bool is_ref,
                   const fm_index_coll& fm_idx,
                   torch::jit::script::Module& TIS_model,
                   const float& minimum_ORF_score,
                   const bool no_filter,
                   const size_t nb_colours,
                   tbb::concurrent_unordered_map<size_t, float>& all_TIS_scores)
{
    // Set as present and not-reverse complement if is_ref. If false positive slipped through, remove
    std::pair<bool, bool> present(true, false);
    if (is_ref)
    {
        present = path_search(nodelist, fm_idx);
        if (!present.first)
        {
            return;
        }
    }

    // initialise path sequence
    std::string path_sequence;

    // here generate vector/map which contains the start/end of each node in basepairs. This will then be used to determine which nodes an ORF traverses
    std::vector<std::vector<size_t>> node_ranges;
    size_t node_start = 0;

    // identify the frames of the stop codons in path
    const int source_node = nodelist.at(0);
    const int sink_node = nodelist.back();

    // create vector with each frame
    std::vector<bool> stop_frames(3, 0);

    // get a reference to the unitig map object
    auto start_um_pair = get_um_data(ccdbg, head_kmer_arr, source_node);
    auto& start_um = start_um_pair.first;
    auto& start_um_data = start_um_pair.second;

    auto end_um_pair = get_um_data(ccdbg, head_kmer_arr, sink_node);
    auto& end_um = end_um_pair.first;
    auto& end_um_data = end_um_pair.second;

    // determine if start or end of the path has a end-unitig. If so, set all frames to true, else calculate correct frames
    if (start_um_data->end_contig(colour_ID) || end_um_data->end_contig(colour_ID)) {
        stop_frames[0] = 1;
        stop_frames[1] = 1;
        stop_frames[2] = 1;
    } else {
        // get codon_arr from start_node
        const std::bitset<3> codon_arr = (source_node >= 0) ? start_um_data->get_codon_arr(true, true, 0)
                                                     : start_um_data->get_codon_arr(true, false, 0);

        // check if each bit is set, starting from 0 and add to stop_frames
        for (size_t i = 0; i < 3; i++) {
            if (codon_arr[i]) {
                stop_frames[i] = 1;
            }
        }
    }

    // generate path sequence by merging nodes in sequence
    for (const auto &node : nodelist) {
        // add node to node list for path coordinates
        // 0th entry is start index of node within the unitig, 1st entry is end index of node within unitig, 2nd is node end, 3rd entry is node length.
        std::vector<size_t> node_range(3);

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

        // calculate length of unitig and get end coordinates
        size_t node_end = seq.size();

        if (path_sequence.empty()) {
            path_sequence = seq;
            // start index of node in path
            node_range[0] = node_start;
            // start index of next node (one past the end index)
            node_range[1] = node_end;
            // absolute end index of node
            node_range[2] = node_end - 1;
            node_start = node_end;
        } else {
            path_sequence.append(seq.begin() + overlap, seq.end());
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

    // generate dictionaries for start and stop codon indices for each frame
    std::unordered_map<size_t, std::vector<size_t>> start_codon_dict;
    std::unordered_map<size_t, std::vector<size_t>> stop_codon_dict;


    // scope for start_codon_indices, stop_codon_indices and frame vectors
    {
        // generate codon indexes using FindIndex function for start and stop codons
        std::vector<size_t> start_codon_indices;
        std::vector<size_t> stop_codon_indices;
        std::vector<std::size_t> found_indices;

        // iterate over each stop codon
        for (const auto &codon : stop_codons) {
            found_indices = findIndex(path_sequence, codon, 0, 0, false);
            stop_codon_indices.insert(stop_codon_indices.end(), std::make_move_iterator(found_indices.begin()),
                                      std::make_move_iterator(found_indices.end()));
            found_indices.clear();
        }

        // iterate over each start codon
        for (const auto &codon : start_codons) {
            found_indices = findIndex(path_sequence, codon, 0, 0, false);
            start_codon_indices.insert(start_codon_indices.end(), std::make_move_iterator(found_indices.begin()),
                                       std::make_move_iterator(found_indices.end()));
            found_indices.clear();
        }


        // sort codon index arrays into ascending order
        std::sort(start_codon_indices.begin(), start_codon_indices.end());
        std::sort(stop_codon_indices.begin(), stop_codon_indices.end());
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
    }

    // create pair for each paired start and stop codon
    std::vector<std::unordered_map<size_t, std::vector<std::pair<size_t, std::pair<std::size_t, std::size_t>>>>> ORF_index_pairs(3);

    // iterate through frames, pair sequential start+stop codons after first stop codon
    // for each stop codon frame, determine first start codon and last stop codon and pair
    for (int frame = 0; frame < stop_frames.size(); frame++)
    {
        if (stop_frames[frame])
        {
            // initialise stop and start index iterators
            auto start_index = start_codon_dict[frame].begin();
            auto stop_index = stop_codon_dict[frame].begin();

            // iterate over start and stop indices until you reach the end
            while (stop_index != stop_codon_dict[frame].end() && start_index != start_codon_dict[frame].end())
            {
                if (*start_index < *stop_index)
                {
                    // stop index is indexed at first base, therefore end of ORF is two bases after
                    std::pair<size_t, size_t> codon_pair(*start_index, *stop_index + 2);
                    size_t ORF_len = codon_pair.second - codon_pair.first + 1;

                    // check if length is large enough, otherwise ignore
                    if (ORF_len >= min_len)
                    {
                        std::pair<size_t, std::pair<std::size_t, std::size_t>> ORF_entry = std::make_pair(ORF_len, codon_pair);

                        // add to correct frame in ORF_index_pairs
                        ORF_index_pairs[frame][codon_pair.second].push_back(std::move(ORF_entry));
                    }

                    // iterate start index
                    start_index++;
                } else {
                    // start is greater than stop, so iterate stop_iterator
                    stop_index++;
                }
            }
        }
    }

    // generate sequences for ORFs from codon pairs
    for (const auto& frame : ORF_index_pairs)
    {
        // check if any ORFs present, if not pass.
        if (!frame.empty())
        {
            // iterate over the entries in each frame, starting with the largest.
            for (auto stop_codon : frame)
            {
                // sort the ORFs with the current stop codon by size
                sort(stop_codon.second.rbegin(), stop_codon.second.rend());

                if (!no_filter)
                {
                    // Get TIS score for all entries
                    std::vector<float> TIS_scores(stop_codon.second.size());
                    std::vector<std::tuple<std::string, std::string, size_t>> TIS_list(stop_codon.second.size());

                    // create vector of all entries
                    for (int i = 0; i < stop_codon.second.size(); i++)
                    {
                        const auto& ORF_entry = stop_codon.second.at(i);
                        const auto& ORF_len = ORF_entry.first;
                        const auto& codon_pair = ORF_entry.second;
                        std::string ORF_seq = path_sequence.substr((codon_pair.first), (ORF_len));
                        std::string TIS_seq = "";
                        if (codon_pair.first >= 16)
                        {
                            TIS_seq = path_sequence.substr((codon_pair.first - 16), 16);
                        }
                        std::string Downstream_seq = ORF_seq.substr (0,19);

                        TIS_list[i] = std::make_tuple(TIS_seq, Downstream_seq, ORF_len);
                    }

                    auto score_vector = score_TIS(TIS_list, TIS_model, minimum_ORF_score, all_TIS_scores);

                    // set variables to determine highest scoring ORF with same stop codon
                    std::pair<size_t, size_t> best_codon = {0,0};
                    size_t best_hash;
                    bool best_TIS_present;
                    size_t best_ORF_len = 0;
                    ORFCoords best_ORF_coords;
                    float best_overall_score = 0;
                    float best_TIS_score = 0;

                    // set maximum ORF length for stop codon frame
                    const size_t max_ORF_length = stop_codon.second.at(0).first;

                    // unpack all ORFs with same stop codon
                    for (int i = 0; i < stop_codon.second.size(); i++)
                    {
                        const auto& ORF_entry = stop_codon.second.at(i);
                        const auto& ORF_len = ORF_entry.first;
                        const auto& codon_pair = ORF_entry.second;
                        const auto& TIS_score = score_vector.at(i);

                        // get ORF coords for current iteration
                        ORFCoords ORF_coords = std::move(calculate_coords(codon_pair, nodelist, node_ranges));

                        // determine first 'complete unitig' i.e. with start not in overlap region
                        int index = 0;
                        for (; index < ORF_coords.second.size(); index++)
                        {
                            if (ORF_coords.second.at(index).second - ORF_coords.second.at(index).first >= overlap)
                            {
                                break;
                            }
                        }

                        // determine coverage for start_site
                        auto um_pair = get_um_data(ccdbg, head_kmer_arr, ORF_coords.first[index]);
                        auto& um_data = um_pair.second;
                        float start_coverage = (float)um_data->full_colour().count() / (float)nb_colours;

                        // calculate length as proportion of previous ORF in list
                        float prev_len_prop = (float)ORF_len / (float)max_ORF_length;

                        // generate score based on start coverage multiplied by dataset size, TIS score and size proportion to previous ORF
                        const float overall_score = start_coverage * TIS_score * prev_len_prop;

                        // determine if score is better and start site is better supported
                        if (overall_score > best_overall_score)
                        {
                            best_codon = codon_pair;
                            best_TIS_present = (codon_pair.first >= 16);
                            best_ORF_len = ORF_len;
                            best_overall_score = overall_score;
                            best_TIS_score = TIS_score;
                            best_ORF_coords = std::move(ORF_coords);

                            // generate hash
                            std::string ORF_seq = path_sequence.substr((codon_pair.first), (ORF_len));
                            std::string TIS_seq = "";
                            if (best_TIS_present)
                            {
                                TIS_seq = path_sequence.substr((codon_pair.first - 16), 16);
                            }
                            ORF_seq = TIS_seq + ORF_seq;
                            best_hash = hasher{}(ORF_seq);
                        }
                    }

                    // check if codon found
                    if (best_codon.second)
                    {
                        // if TIS is present, add the non-TIS hash to to_remove
                        if (best_TIS_present)
                        {
                            size_t hash_to_remove = hasher{}(path_sequence.substr((best_codon.first), (best_ORF_len)));
                            hashes_to_remove.insert(hash_to_remove);
                        }

                        // if ORF looks unlikely to be real, do not add
                        if (best_overall_score)
                        {
                            // create ORF_node_vector, populate with results from node traversal (add true on end for relative strand if !is_ref).
                            ORFNodeVector ORF_node_vector = std::make_tuple(best_ORF_coords.first, best_ORF_coords.second, best_ORF_len, !present.second, best_TIS_score);

                            // think about if there is no TIS, then can ignore ORF?
                            update_ORF_node_map(ccdbg, head_kmer_arr, best_hash, ORF_node_vector, ORF_node_map);
                        }
                    }
                } else
                {
                    // unpack all ORFs with same stop codon
                    for (const auto& ORF_entry : stop_codon.second)
                    {
                        const auto& ORF_len = ORF_entry.first;
                        const auto& codon_pair = ORF_entry.second;

                        // get ORF seqeunce and pull 16bp upstream of start codon for TIS model if possible. If not, do not and set TIS_present as false
                        std::string TIS_seq;
                        if (codon_pair.first >= 16)
                        {
                            TIS_seq = path_sequence.substr((codon_pair.first - 16), 16);
                            size_t hash_to_remove = hasher{}(path_sequence.substr((codon_pair.first), (ORF_len)));
                            hashes_to_remove.insert(hash_to_remove);
                        }

                        // generate hash for ORF sequence and ORF sequence
                        std::string ORF_seq = path_sequence.substr((codon_pair.first), (ORF_len));

                        ORF_seq = TIS_seq + ORF_seq;

                        size_t ORF_hash = hasher{}(ORF_seq);

                        // If ORF is real, continue and work out coordinates for ORF in node space
                        ORFCoords ORF_coords = std::move(calculate_coords(codon_pair, nodelist, node_ranges));

                        // create ORF_node_vector, populate with results from node traversal (add true on end for relative strand if !is_ref).
                        ORFNodeVector ORF_node_vector = std::make_tuple(ORF_coords.first, ORF_coords.second, ORF_len, !present.second, 0);

                        // think about if there is no TIS, then can ignore ORF?
                        update_ORF_node_map(ccdbg, head_kmer_arr, ORF_hash, ORF_node_vector, ORF_node_map);
                    }
                }
            }
        }
    }
}

// calculate node coordinates for an ORF
ORFCoords calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                            const std::vector<int>& nodelist,
                            const std::vector<std::vector<size_t>>& node_ranges)
{
    // initialise items for a tuple containing a vector of each node name, corresponding vector of positions traversed in the node and node strand
    std::vector<int> ORF_node_id;
    std::vector<indexPair> ORF_node_coords;

    for (size_t i = 0; i < nodelist.size(); i++)
    {
        unsigned int traversed_node_start;
        unsigned int traversed_node_end;
        bool start_assigned = false;
        bool end_assigned = false;

        // if start of ORF is below node range, then check ORF traverses node
        if (codon_pair.first < node_ranges[i][0])
        {
            traversed_node_start = 0;
            start_assigned = true;
        } else if (codon_pair.first >= node_ranges[i][0] && codon_pair.first < node_ranges[i][1]){
            traversed_node_start = codon_pair.first - node_ranges[i][0];
            start_assigned = true;
        }

        // if end of ORF is above node range, then check if ORF traversed
        if (codon_pair.second >= node_ranges[i][1])
        {
            traversed_node_end = node_ranges[i][2];
            end_assigned = true;
        } else if (codon_pair.second >= node_ranges[i][0] && codon_pair.second < node_ranges[i][1])
        {
            traversed_node_end = codon_pair.second - node_ranges[i][0];
            end_assigned = true;
        }

        // if the ORF traverses node, update coordinates
        if (start_assigned && end_assigned)
        {
            indexPair node_coords = std::make_pair(traversed_node_start, traversed_node_end);
            ORF_node_id.push_back(nodelist[i]);
            ORF_node_coords.push_back(std::move(node_coords));
        }
        // gone past last node covering ORF so assign 3p and end index to previous node
        else if (start_assigned && !end_assigned)
        {
            break;
        }
    }

    const ORFCoords return_tuple = std::make_pair(ORF_node_id, ORF_node_coords);
    return return_tuple;
}

// updates ORF_node_map if identical ORFs with different node sets detected
void update_ORF_node_map (const ColoredCDBG<MyUnitigMap>& ccdbg,
                          const std::vector<Kmer>& head_kmer_arr,
                          const size_t& ORF_hash,
                          ORFNodeVector& ORF_node_vector,
                          ORFNodeMap& ORF_node_map)
{
    // if not present, add as is
    if (ORF_node_map.find(ORF_hash) == ORF_node_map.end())
    {
        ORF_node_map.emplace(ORF_hash, std::move(ORF_node_vector));
    } else
    {
        //reference entry in ORF_node_map
        auto& current_entry = ORF_node_map[ORF_hash];

        // if present, check which of the ORFs has the longest path in terms ORF,
        // as used for overlap information, then update path information
        size_t current_ORF_size = std::get<0>(current_entry).size();
        size_t new_ORF_size = std::get<0>(ORF_node_vector).size();

        // if more information stored in new ORF entry, replace and update with new information
        if (new_ORF_size > current_ORF_size)
        {
            current_entry = std::move(ORF_node_vector);
        } else if (new_ORF_size == current_ORF_size)
        {
            // if ORFs are also the same length, check if first/last entries are identical across the two entries
            // to deal with cases where a different node is traversed at beginning/end
            // for first, check if TIS present, if so take that as first entry, if not take first ORF node
            const auto& current_first = std::get<0>(current_entry).at(0);
            const auto& new_first = std::get<0>(ORF_node_vector).at(0);
            const auto& current_last = std::get<0>(current_entry).back();
            const auto& new_last = std::get<0>(ORF_node_vector).back();

            // check if the nodes are identical, if so pass
            if (current_first == new_first && current_last == new_last)
            {
                return;
            }

            // get a reference to the unitig map object
            auto curr_um_pair = get_um_data(ccdbg, head_kmer_arr, current_first);
            auto& curr_um = curr_um_pair.first;

            auto new_um_pair = get_um_data(ccdbg, head_kmer_arr, new_first);
            auto& new_um = new_um_pair.first;

            // get hashes of first nodes to ensure stable assignment across runs
            const size_t current_first_hash = hasher{}(curr_um.getUnitigHead().toString());
            const size_t new_first_hash = hasher{}(new_um.getUnitigHead().toString());

            // compare first/last entries. Assign the new ORF entry to that with the highest ID
            if (new_first_hash > current_first_hash)
            {
                current_entry = std::move(ORF_node_vector);
            } else if (new_first_hash == current_first_hash)
            {
                // get a reference to the unitig map object
                curr_um_pair = get_um_data(ccdbg, head_kmer_arr, current_last);
                curr_um = curr_um_pair.first;

                new_um_pair = get_um_data(ccdbg, head_kmer_arr, new_last);
                new_um = new_um_pair.first;

                // if no clear winner, check last nodes
                const size_t current_last_hash = hasher{}(curr_um.getUnitigHead().toString());
                const size_t new_last_hash = hasher{}(new_um.getUnitigHead().toString());

                if (new_last_hash > current_last_hash)
                {
                    current_entry = std::move(ORF_node_vector);
                }
            }
        }
    }
}

// converts ORF entries into a vector and assigns relative strand
ORFNodeRobMap sort_ORF_indexes(ORFNodeMap& ORF_node_map,
                           const NodeStrandMap& pos_strand_map,
                           const ColoredCDBG<MyUnitigMap>& ccdbg,
                           const std::vector<Kmer>& head_kmer_arr,
                           const bool is_ref)
{
    ORFNodeRobMap ORF_map;

    // generate ID for each ORF and assign strand
    size_t ORF_ID = 0;
    for (auto& ORF : ORF_node_map)
    {
        // if strand not known from contigs, assign based on ORF overlap
        if (!is_ref)
        {
            // assign strand to ORF by iterating over all nodes, and assigning strand to most supported by nodes
            int num_pos = 0;
            int num_neg = 0;
            // iterate over ORF nodes...
            for (const auto& node_id : std::get<0>(ORF.second))
            {
                const bool strand = (node_id > 0) ? true : false;

                // get a reference to the unitig map object
                auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
                auto& um = um_pair.first;

                const size_t node_hash = hasher{}(um.getUnitigHead().toString());
                if (strand != pos_strand_map.at(node_hash))
                {
                    num_neg++;
                } else
                {
                    num_pos++;
                }
            }

            // if the number of negative strands are greater than positive, assign overall strand as negative
            if (num_neg > num_pos)
            {
                std::get<3>(ORF.second) = false;
            }
        }

        // move entry from map to vector
        ORF_map[ORF_ID] = std::move(ORF.second);

        // iterate ORF_ID
        ORF_ID++;
    }

    // clear ORF_node_map
    ORF_node_map.clear();

    return ORF_map;
}

// calculate the relative strand of each node traversed in an ORF
NodeStrandMap calculate_pos_strand(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                   const std::vector<Kmer>& head_kmer_arr,
                                   const ORFNodeMap& ORF_node_map)
{
    // initialise map to store orientation of nodes (colour is an ID, not a string
    std::vector<NodeStrandMap> pos_strand_vector;

    for (const auto& ORF : ORF_node_map)
    {
        // create new map to store node info
        NodeStrandMap new_map;

        // repeat for ORF nodes
        for (const auto& node_id : std::get<0>(ORF.second))
        {
            // get a reference to the unitig map object
            auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
            auto& um = um_pair.first;

            const size_t node_hash = hasher{}(um.getUnitigHead().toString());
            // assigned new_map entry. If node is positive, strand is true, if negative, strand is false
            new_map[node_hash] = (node_id > 0) ? true : false;
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
            // else, go through and add the new map entries to the existing maps in pos_strand_vector
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

    // iterate over pos_strand_vector, merging any maps that have matching ORFs. If no matching ORFs, merge as is. Continue until only single item present in pos_strand_vector
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

    // if no ORFs present, need to just return an empty NodeStrandMap
    if (pos_strand_vector.empty())
    {
        NodeStrandMap empty_map;
        return empty_map;
    } else
    {
        // return only first entry of pos_strand_vector, as this structure is now flat
        return pos_strand_vector[0];
    }
}