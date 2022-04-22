// ggCaller header
#include "traversal.h"
PathVector iter_nodes_binary (const ColoredCDBG<MyUnitigMap>& ccdbg,
                              const std::vector<Kmer>& head_kmer_arr,
                              const NodeTuple& head_node_tuple,
                              const size_t current_colour,
                              const size_t length_max,
                              const size_t overlap,
                              const bool repeat,
                              const bool is_ref,
                              const boost::dynamic_bitset<>& ref_set,
                              const fm_index_coll& fm_idx)
{
    // generate path list, vector for path and the stack
    PathVector path_list;
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
        const auto & codon_arr = std::get<2>(node_tuple);
        const boost::dynamic_bitset<>& colour_arr = std::get<3>(node_tuple);
        const size_t & path_length = std::get<4>(node_tuple);

        // slice path, unless at first node
        if (pos_idx != 0)
        {
            node_vector = std::vector<int> (node_vector.begin(), node_vector.begin() + pos_idx);
        }

        // add node to path
        node_vector.push_back(node_id);

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        // get unitig data
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
//            if (!is_ref)
//            {
//                if (neighbour_colour.count() < 2)
//                {
//                    continue;
//                }
//            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + (neighbour_um.size - overlap);

            if (updated_path_length > length_max)
            {
                // if path is only 2 nodes in length return as will contain very large sequences
                if (new_pos_idx == 1)
                {
                    // check if end node is reference sequence
                    if (!is_ref)
                    {
                        auto ref_check = neighbour_colour;
                        ref_check &= ref_set;
                        if (ref_check.none())
                        {
                            continue;
                        }
                    }

                    // create temporary path to account for reaching end of contig
                    std::vector<int> return_path = node_vector;
                    return_path.push_back(neighbour_id);

                    // update path_list to enable this to be returned
                    path_list.push_back(std::move(return_path));
                }
                continue;
            }

            // calculate modulus for path and updated cached path reading frame
            int modulus = path_length % 3;
            std::bitset<3> updated_codon_arr = (codon_arr & ~(neighbour_um_data->get_codon_arr(false, neighbour_strand, modulus)));

            // return path if end of contig found or if stop indexes paired
            if (neighbour_um_data->end_contig(current_colour) || updated_codon_arr.none())
            {
                // check if end node is reference sequence
                if (!is_ref)
                {
                    auto ref_check = neighbour_colour;
                    ref_check &= ref_set;
                    if (ref_check.none())
                    {
                        continue;
                    }
                }

                // create temporary path to account for reaching end of contig
                std::vector<int> return_path = node_vector;
                return_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(return_path));

                // if node is end_contig, continue to add to path
                if (neighbour_um_data->end_contig(current_colour) && !updated_codon_arr.none())
                {
                    // if no previous conditions are satisfied, prepare tuple for stack
                    NodeTuple new_node_tuple(new_pos_idx, neighbour_id, updated_codon_arr, updated_colours_arr, updated_path_length);

                    // add to stack
                    node_stack.push(new_node_tuple);
                }

                // move onto next neighbour
                continue;
            }

            // if no previous conditions are satisfied, prepare tuple for stack
            NodeTuple new_node_tuple(new_pos_idx, neighbour_id, updated_codon_arr, updated_colours_arr, updated_path_length);

            // add to stack
            node_stack.push(new_node_tuple);
        }
    }
    return path_list;
}

ORFNodeRobMap traverse_graph(const ColoredCDBG<MyUnitigMap>& ccdbg,
                         const std::vector<Kmer>& head_kmer_arr,
                         const size_t colour_ID,
                         const std::vector<size_t>& node_ids,
                         const bool repeat,
                         const size_t max_path_length,
                         const size_t overlap,
                         const bool is_ref,
                         const boost::dynamic_bitset<>& ref_set,
                         const fm_index_coll& fm_idx,
                         const std::vector<std::string>& stop_codons_for,
                         const std::vector<std::string>& start_codons_for,
                         const size_t min_ORF_length,
                         torch::jit::script::Module& ORF_model,
                         const double& minimum_ORF_score,
                         const bool no_filter,
                         tbb::concurrent_unordered_map<size_t, float>& all_ORF_scores)
{
    //initialise ORF_nodes_paths to add ORF sequences to
    ORFNodeMap ORF_node_map;
    std::unordered_set<size_t> hashes_to_remove;

    // set for any end contigs
    std::bitset<3> full_binary;
    full_binary.set();

    // traverse nodes in forward direction
    for (const auto& node_id : node_ids)
    {
        // get unitig data
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // check if stop codons present. If not, pass
        if (!um_data->forward_stop() && !um_data->end_contig(colour_ID))
        {
            continue;
        }

        // generate integer version of unitig_id for recursion
        const int head_id = (int) node_id;

        // gather unitig information from graph_vector
        std::bitset<3> codon_arr = um_data->get_codon_arr(true, true, 0);
        const size_t unitig_len = um.size;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // if end contig, set to full array
        if (um_data->end_contig(colour_ID))
        {
            codon_arr = full_binary;
        }

        // check if end node is reference sequence
        if (!is_ref)
        {
            auto ref_check = colour_arr;
            ref_check &= ref_set;
            if (ref_check.none())
            {
                continue;
            }
        }

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, overlap, repeat, is_ref, ref_set, fm_idx);

        // iterate over paths, calling ORFs
        if (!unitig_complete_paths.empty())
        {
            // iterate over all_paths
            for (int i = 0; i < unitig_complete_paths.size(); i++)
            {
                // generate all ORFs within the path for start and stop codon pairs
                generate_ORFs(colour_ID, ORF_node_map, hashes_to_remove, ccdbg, head_kmer_arr, stop_codons_for, start_codons_for, unitig_complete_paths[i], overlap, min_ORF_length, is_ref, fm_idx, ORF_model, minimum_ORF_score, no_filter, all_ORF_scores);
            }
        }
    }

    // traverse nodes in reverse direction
    for (const auto& node_id : node_ids)
    {
        // get unitig data
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // check if stop codons present. If not, pass
        if (!um_data->reverse_stop() && !um_data->end_contig(colour_ID))
        {
            continue;
        }

        // generate integer version of unitig_id for recursion, multiplied by -1 to indicate reversal
        const int head_id = ((int) node_id) * -1;

        // gather unitig information from graph_vector
        std::bitset<3> codon_arr = um_data->get_codon_arr(true, false, 0);
        const size_t unitig_len = um.size;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // if end contig, set to full array
        if (um_data->end_contig(colour_ID))
        {
            codon_arr = full_binary;
        }

        // check if end node is reference sequence
        if (!is_ref)
        {
            auto ref_check = colour_arr;
            ref_check &= ref_set;
            if (ref_check.none())
            {
                continue;
            }
        }

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, overlap, repeat, is_ref, ref_set, fm_idx);

        // iterate over paths, calling ORFs
        if (!unitig_complete_paths.empty())
        {
            // iterate over all_paths
            for (int i = 0; i < unitig_complete_paths.size(); i++)
            {
                // generate all ORFs within the path for start and stop codon pairs
                generate_ORFs(colour_ID, ORF_node_map, hashes_to_remove, ccdbg, head_kmer_arr, stop_codons_for, start_codons_for, unitig_complete_paths[i], overlap, min_ORF_length, is_ref, fm_idx, ORF_model, minimum_ORF_score, no_filter, all_ORF_scores);
            }
        }
    }

    // remove hashes to remove from ORF_hash
    for (const auto& hash : hashes_to_remove)
    {
        auto it = ORF_node_map.find(hash);
        if (it != ORF_node_map.end())
        {
            ORF_node_map.erase(it);
        }
    }

    // generate pos_strand_map to determine relative strands of each node for each colour if not given by alignment of contigs
    NodeStrandMap pos_strand_map;
    if (!is_ref)
    {
        pos_strand_map = std::move(calculate_pos_strand(ccdbg, head_kmer_arr, ORF_node_map));
    }

    // group colours of ORFs together
    ORFNodeRobMap ORF_map = std::move(sort_ORF_indexes(ORF_node_map, pos_strand_map, ccdbg, head_kmer_arr, is_ref));

    return ORF_map;
}