// ggCaller header
#include "traversal.h"
// define mutex for safe addition to robinhood_maps
//std::mutex mtx2;

PathVector iter_nodes_binary (const ColoredCDBG<MyUnitigMap>& ccdbg,
                              const std::vector<Kmer>& head_kmer_arr,
                              const NodeTuple& head_node_tuple,
                              const size_t& current_colour,
                              const size_t& length_max,
                              const size_t& overlap,
                              const bool& repeat,
                              const bool& is_ref)
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
        const uint8_t & codon_arr = std::get<2>(node_tuple);
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
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // iterate over neighbours, recurring through incomplete paths
        for (const auto& neighbour : um_data->get_neighbours(strand))
        {
            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const auto& neighbour_id = std::get<0>(neighbour);
            const auto& frame = std::get<1>(neighbour);
            const auto& colour_set = std::get<2>(neighbour);

            // if is_ref, determine if edge is correct
            if (is_ref)
            {
                if (colour_set.find(current_colour) == colour_set.end())
                {
                    continue;
                }
            }

            // check if unitig has already been traversed, and pass if repeat not specified
            const bool is_in = node_set.find(neighbour_id) != node_set.end();
            if (!repeat && is_in)
            {
                continue;
            }

            // get reference to unitig_dict object for neighbour
            const auto neighbour_um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
            const auto& neighbour_um = neighbour_um_pair.first;
            const auto& neighbour_um_data = neighbour_um_pair.second;

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            updated_colours_arr &= neighbour_um_data->full_colour();

            // determine if neighbour is in same colour as iteration, if not pass
            if (!updated_colours_arr[current_colour])
            {
                continue;
            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + (neighbour_um.size - overlap);

            if (updated_path_length > length_max)
            {
                // if path is only 2 nodes in length return as will contain very large sequences
                if (new_pos_idx == 1)
                {
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
            uint8_t updated_codon_arr = (codon_arr & ~(frame.at(modulus)));

            // return path if end of contig found or if stop indexes paired
            if (neighbour_um_data->end_contig() || updated_codon_arr == 0)
            {
                // create temporary path to account for reaching end of contig
                std::vector<int> return_path = node_vector;
                return_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(return_path));

                // move onto next neighbour
                continue;
            }

            // if no previous conditions are satisfied, prepare tuple for stack
            NodeTuple new_node_tuple(new_pos_idx, neighbour_id, updated_codon_arr, updated_colours_arr, updated_path_length);

            // add to stack
            node_stack.push(new_node_tuple);

            // add node to node_set
            node_set.insert(neighbour_id);
        }
    }
    return path_list;
}

std::vector<PathVector> traverse_graph(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                       const std::vector<Kmer>& head_kmer_arr,
                                       const size_t& colour_ID,
                                       const std::vector<size_t>& node_ids,
                                       const bool repeat,
                                       const size_t max_path_length,
                                       const size_t& overlap,
                                       const bool is_ref)
{
    // initialise all_paths
    std::vector<PathVector> all_paths;

    // traverse nodes in forward direction
    for (const auto& node_id : node_ids)
    {
        // get unitig data
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // check if stop codons present. If not, pass
        if (!um_data->forward_stop())
        {
            continue;
        }

        // generate integer version of unitig_id for recursion
        const int head_id = (int) node_id;

        // gather unitig information from graph_vector
        const uint8_t codon_arr = um_data->get_codon_arr(true, true, 0);
        const size_t unitig_len = um.size;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, overlap, repeat, is_ref);

        if (!unitig_complete_paths.empty())
        {
            all_paths.push_back(std::move(unitig_complete_paths));
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
        if (!um_data->reverse_stop())
        {
            continue;
        }

        // generate integer version of unitig_id for recursion, multiplied by -1 to indicate reversal
        const int head_id = ((int) node_id) * -1;

        // gather unitig information from graph_vector
        const uint8_t codon_arr = um_data->get_codon_arr(true, false, 0);
        const size_t unitig_len = um.size;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, overlap, repeat, is_ref);

        if (!unitig_complete_paths.empty())
        {
            all_paths.push_back(std::move(unitig_complete_paths));
        }
    }

    return all_paths;
}