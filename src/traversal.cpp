// ggCaller header
#include "traversal.h"
// define mutex for safe addition to robinhood_maps
//std::mutex mtx2;

PathVector iter_nodes_binary (const GraphVector& graph_vector,
                              const NodeTuple& head_node_tuple,
                               const size_t& current_colour,
                               const size_t& length_max,
                               const bool& repeat)
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
        const sdsl::bit_vector & colour_arr = std::get<3>(node_tuple);
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

        // get unitig_dict entry in graph_vector
        const auto& node_dict = graph_vector.at(abs(node_id) - 1);

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // iterate over neighbours, recurring through incomplete paths
        for (const auto& neighbour : node_dict.get_neighbours(strand))
        {
            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const auto& neighbour_id = neighbour.first;
            const auto& frame = neighbour.second;

            // check if unitig has already been traversed, and pass if repeat not specified
            const bool is_in = node_set.find(neighbour_id) != node_set.end();
            if (!repeat && is_in)
            {
                continue;
            }

            // get reference to unitig_dict object for neighbour
            const auto& neighbour_dict = graph_vector.at(abs(neighbour_id) - 1);

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            updated_colours_arr &= neighbour_dict.full_colour();

            // determine if neighbour is in same colour as iteration, if not pass
            if (!updated_colours_arr[current_colour])
            {
                continue;
            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + neighbour_dict.size().second;

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
            if (neighbour_dict.end_contig() || updated_codon_arr == 0)
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

AllPaths traverse_graph(const GraphVector& graph_vector,
                         const size_t& colour_ID,
                         const std::vector<size_t>& node_ids,
                         const bool repeat,
                         const size_t max_path_length)
{
    // initialise all_paths
    AllPaths all_paths;

    // traverse nodes in forward direction
    for (const auto& node_id : node_ids)
    {
        // parse unitig_id. Zero based, so take 1
        const auto unitig_id = node_id - 1;

        // get reference to unitig_dict object
        const auto& unitig_dict = graph_vector.at(unitig_id);

        // check if stop codons present. If not, pass
        if (!unitig_dict.forward_stop())
        {
            continue;
        }

        // generate integer version of unitig_id for recursion
        const int head_id = (int) node_id;

        // gather unitig information from graph_vector
        const uint8_t codon_arr = unitig_dict.get_codon_arr(true, true, 0);
        const size_t unitig_len = unitig_dict.size().first;
        const sdsl::bit_vector colour_arr = unitig_dict.full_colour();

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(graph_vector, head_node_tuple, colour_ID, max_path_length, repeat);

        if (!unitig_complete_paths.empty())
        {
            all_paths.push_back(std::move(unitig_complete_paths));
        }
    }

    // traverse nodes in reverse direction
    for (const auto& node_id : node_ids)
    {
        // parse unitig_id. Zero based, so take 1
        const auto unitig_id = node_id - 1;

        // get reference to unitig_dict object
        const auto& unitig_dict = graph_vector.at(unitig_id);

        // check if stop codons present. If not, pass
        if (!unitig_dict.reverse_stop())
        {
            continue;
        }

        // generate integer version of unitig_id for recursion, multiplied by -1 to indicate reversal
        const int head_id = ((int) node_id) * -1;

        // gather unitig information from graph_vector
        const uint8_t codon_arr = unitig_dict.get_codon_arr(true, false, 0);
        const size_t unitig_len = unitig_dict.size().first;
        const sdsl::bit_vector colour_arr = unitig_dict.full_colour();

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(graph_vector, head_node_tuple, colour_ID, max_path_length, repeat);

        if (!unitig_complete_paths.empty())
        {
            all_paths.push_back(std::move(unitig_complete_paths));
        }
    }

    return all_paths;
}

std::vector<std::pair<size_t, size_t>> check_next_ORFs (const GraphVector& graph_vector,
                                                        const int& head_node,
                                                        const size_t& stream_source,
                                                        const size_t& current_colour,
                                                        const int& stream,
                                                        const ORFVector& ORF_vector,
                                                        const std::unordered_set<size_t>& uninode_ORFs,
                                                        std::unordered_set<int>& prev_node_set)
{
    // initialise return vector of upstream ORFs (vectors will break at branches in graph)
    std::vector<std::pair<size_t, size_t>> connected_ORFs;

    // set fully_traversed to false
    //bool fully_traversed = false;

    // get the colour array of the head node
    auto colour_arr = graph_vector.at(abs(head_node) - 1).full_colour();

    // check if start ORF is uninode, and therefore should be added to fully_traversed
//    if (uninode_ORFs.find(stream_source) != uninode_ORFs.end())
//    {
//        fully_traversed = true;
//    }

    // generate path list, vector for path and the stack
    ORFStack ORF_stack;

    // create node set for identification of repeats
    std::unordered_set<int> node_set;
    node_set.insert(head_node * stream);

    // create first item in stack, multiply by stream - upstream (stream = -1) or downstream (stream = 1) and add stream ORF
    ORF_stack.push(std::make_tuple(head_node * stream, node_set, colour_arr));

    while(!ORF_stack.empty())
    {
        // pop node in stack
        auto path_tuple = ORF_stack.top();
        const auto& node_id = std::get<0>(path_tuple);
        node_set = std::get<1>(path_tuple);
        colour_arr = std::get<2>(path_tuple);
        ORF_stack.pop();

        // get unitig_dict entry in graph_vector
        const auto& node_dict = graph_vector.at(abs(node_id) - 1);

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // iterate over neighbours
        for (const auto& neighbour : node_dict.get_neighbours(strand))
        {
            // make copy of node_set
            auto temp_node_set = node_set;

            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const auto& neighbour_id = neighbour.first;

            // check if unitig has already been traversed, and pass if repeat not specified
            const bool is_in = temp_node_set.find(neighbour_id) != temp_node_set.end();
            if (is_in)
            {
                continue;
            }

            // get reference to unitig_dict object for neighbour
            const auto& neighbour_dict = graph_vector.at(abs(neighbour_id) - 1);

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            updated_colours_arr &= neighbour_dict.full_colour();

//            // test for colour vector
//            std::vector<bool> updated_colours_arr_real(updated_colours_arr.size());
//            std::vector<bool> colours_arr_real(updated_colours_arr.size());
//            std::vector<bool> neighbour_dict_colours_arr_real(updated_colours_arr.size());
//            for (int i = 0; i < updated_colours_arr.size(); i++)
//            {
//                updated_colours_arr_real[i] = (updated_colours_arr[i]) ? true : false;
//                colours_arr_real[i] = (colour_arr[i]) ? true : false;
//                neighbour_dict_colours_arr_real[i] = ((neighbour_dict.full_colour())[i]) ? true : false;
//            }

            // determine if neighbour is in same colour as iteration, if not pass
            if (!updated_colours_arr[current_colour])
            {
                continue;
            }

            // add node to temp_node_set
            temp_node_set.insert(neighbour_id);

            // check if node is traversed by end of an ORF
            if (!neighbour_dict.ORFs_empty(current_colour))
            {
                // pull out next ORFs and order them
                const auto& next_ORFs = neighbour_dict.get_ORFs(current_colour);
                const auto ordered_ORFs = order_ORFs(graph_vector, next_ORFs, neighbour_id, ORF_vector);

                // bool for determining if paired ORFs are the same, and if they are already connected as end_ORFs
                bool skip_pair = false;

                // make sure that 5' and 3' of ORF is traversed if it is not uninode
                if (std::find(ordered_ORFs.begin(), ordered_ORFs.end(), stream_source) != ordered_ORFs.end())
                {
                    //fully_traversed = true;
                    // check if stream source is at the back of the vector, if so, then need to continue traversing
                    if (ordered_ORFs.back() == stream_source)
                    {
                        skip_pair = true;
                    }
                }

                // if the node is not fully traversed, do not add the edge
                // connecting this and the next, but do add for remaining in ordered_ORFs.
                if (!skip_pair)
                {
                    // add stream_source and first entry
                    connected_ORFs.push_back({stream_source, ordered_ORFs.at(0)});

                    // update previous nodes encountered
                    prev_node_set.insert(std::make_move_iterator(temp_node_set.begin()), std::make_move_iterator(temp_node_set.end()));
                } else
                {
                    ORF_stack.push({neighbour_id, temp_node_set, updated_colours_arr});
                }

                // add remaining ordered ORFs to connected_nodes vector
                for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
                {
                    // if ORFs are overlapping, don't want to add connections between an ORF and end of stream source
                    if (ordered_ORFs.at(i + 1) != stream_source)
                    {
                        connected_ORFs.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
                    }
                }
            } else
            {
                // add to stack
                ORF_stack.push({neighbour_id, temp_node_set, updated_colours_arr});
            }
        }
    }
    return connected_ORFs;
}

template <class T>
void clear_stack(std::stack<T>& to_clear)
{
    while (!to_clear.empty())
    {
        to_clear.pop();
    }
}