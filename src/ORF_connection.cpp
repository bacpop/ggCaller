#include "ORF_connection.h"

void add_ORF_info (GraphVector& graph_vector,
                  const size_t& colour_ID,
                  const std::vector<size_t>& target_ORFs,
                  const ORFVector& ORF_vector)
{
    for (const auto & source : target_ORFs)
    {
        {
            // get graph information for source node
            const auto & ORF_info = ORF_vector.at(source);

            // get start node IDs (-1 as graph is zero-based). Only add 5 prime node.
            const auto source_node_id = abs(std::get<0>(ORF_info).at(0)) - 1;
            const auto sink_node_id = abs(std::get<0>(ORF_info).back()) - 1;

            // add ORF information to graph
            graph_vector.at(source_node_id).set_ORFs(colour_ID, source);
            // check if ORF is present only on single node
            if (source_node_id != sink_node_id)
            {
                graph_vector.at(sink_node_id).set_ORFs(colour_ID, source);
            }
        }
    }
}

std::vector<size_t> order_ORFs_in_node(const GraphVector& graph_vector,
                                       const std::unordered_set<size_t>& node_ORFs,
                                       const int& node_id,
                                       const ORFVector& ORF_vector)
{
    std::vector<size_t> overlapping_ORF_IDs;
    std::vector<std::pair<bool, indexPair>> overlapping_ORF_coords;
    // iterate over ORFs in node_ORFs
    for (const auto& ORF_ID : node_ORFs)
    {
        // add to return vector overlapping ORF vector
        overlapping_ORF_IDs.push_back(ORF_ID);

        // get reference to ORF_vector entry
        const auto & ORF_info = ORF_vector.at(ORF_ID);

        // get the index of the node in ORFNodeVector for that ORF
        auto it = find(std::get<0>(ORF_info).begin(), std::get<0>(ORF_info).end(), node_id);

        // if not present, search for reversed node
        if (it == std::get<0>(ORF_info).end())
        {
            it = find(std::get<0>(ORF_info).begin(), std::get<0>(ORF_info).end(), (node_id * -1));
        }

        // get strand from sign of node id (true if positive, false if negative)
        bool strand = (*it > 0) ? true : false;

        // get index of node in ORF coords
//        size_t index = std::distance(std::get<0>(ORF_info).begin(), it);
        size_t index = it - std::get<0>(ORF_info).begin();

        // add coords for node traversal to overlapping_ORF_coords
        overlapping_ORF_coords.push_back(std::pair<bool, indexPair>(strand, std::get<1>(ORF_info).at(index)));
    }

    // initialse vector to work out order of nodes
    std::vector<std::pair<size_t,size_t>> ordered_ORFs;

    // may not be upstream ORFs called, so check if overlapping_ORF_coords is empty
    if (!overlapping_ORF_coords.empty())
    {
        // ensure all coordinates are in the same strand, set the strand of the target node
        bool overall_strand = (node_id > 0) ? true : false;

        // get length of node if reversal is needed
        size_t node_end = graph_vector.at(abs(node_id) - 1).size().first - 1;

        // iterate over entries and flip coords if needed (ignore first as this is the reference)
        for (int i = 0; i < overlapping_ORF_coords.size(); i++)
        {
            if (overlapping_ORF_coords.at(i).first != overall_strand)
            {
                // get difference from original end to absolute last node index
                size_t reversed_end = node_end - overlapping_ORF_coords.at(i).second.first;
                // get difference from original end to absolute last node index
                size_t reversed_start = node_end - overlapping_ORF_coords.at(i).second.second;
                // reassigned the entry in-place in ORF2_nodes.second
                overlapping_ORF_coords[i].second = std::make_pair(reversed_start, reversed_end);
            }

            // add to ordered_ORFs
            ordered_ORFs.push_back(std::pair<size_t,size_t>(overlapping_ORF_coords.at(i).second.first, overlapping_ORF_IDs.at(i)));
        }

        // sort based on first entry in ordered_ORFs
        sort(ordered_ORFs.begin(), ordered_ORFs.end());

        // return ordered_ORFs with ORF_IDs only
        overlapping_ORF_IDs.clear();
        std::transform(ordered_ORFs.begin(), ordered_ORFs.end(), std::back_inserter(overlapping_ORF_IDs),
                       [] (auto const& pair) {return pair.second; });
    }

    return overlapping_ORF_IDs;
}

std::vector<std::pair<size_t, size_t>> pair_ORF_nodes (const GraphVector& graph_vector,
                                                      const size_t& colour_ID,
                                                      const std::vector<size_t>& target_ORFs,
                                                      const ORFVector& ORF_vector,
                                                      const size_t& max_ORF_path_length,
                                                      const int& stream,
                                                      std::unordered_set<int>& prev_node_set,
                                                      const bool is_ref)
{
    // initialise pair of vectors (first = upstream of start_ORF, second = downstream of start_ORF)
    std::vector<std::pair<size_t, size_t>> ORF_edges;

    // iterate over each entry in end_ORFs
    for (const auto& source : target_ORFs)
    {
        // get ORF info
        const auto & source_info = ORF_vector.at(source);

        // get id for source_node and sink node for ORF
        auto source_node_id = std::get<0>(source_info).at(0);
        auto sink_node_id = std::get<0>(source_info).back();

        // decide whether to traverse from source (upstream) or sink (downstream)
        int start_node = (stream > 0) ? sink_node_id : source_node_id;

        // check if node has been traversed previously in opposite direction i.e. ORF has already been paired up/downstream
        if (prev_node_set.find(start_node * stream * -1) == prev_node_set.end())
        {
            // get source_node and sink_node info
            const auto& start_node_info = graph_vector.at(abs(start_node) - 1);
//            const auto& sink_node_info = graph_vector.at(abs(sink_node_id) - 1);

            // check if there are any overlapping ORFs on current node and order...
            const auto& start_node_ORFs = start_node_info.get_ORFs(colour_ID);
//            const auto& sink_node_ORFs = sink_node_info.get_ORFs(colour_ID);


            // check if there are overlapping ORFs on the source and sink node, as these will not be picked up
            std::vector<size_t> ordered_ORFs = order_ORFs_in_node(graph_vector, start_node_ORFs, start_node, ORF_vector);
            if (ordered_ORFs.size() >= 2)
            {
                // if 2 or more ORFs overlapping, order and then add each connecting edge
                for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
                {
                    ORF_edges.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
                }
            }

            // now traverse graph using DFS, finding next ORFs upstream and downstream of source and sink node respectively
            // check that sink hasn't been traversed in reverse (downstream) or source hasn't been traverse forward (downstream

            // if going downstream (stream > 0), check that ORF is last, or upstream (stream < 0), check it is first
            if ((stream > 0 && source == ordered_ORFs.back()) || stream < 0 && source == ordered_ORFs.at(0))
            {
                auto next_ORFs = check_next_ORFs(graph_vector, start_node, source, colour_ID, stream, ORF_vector, max_ORF_path_length, prev_node_set, is_ref);
                ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
            }

            prev_node_set.insert(start_node);
        }
    }
    return ORF_edges;
}

std::vector<std::pair<size_t, size_t>> check_next_ORFs (const GraphVector& graph_vector,
                                                        const int& head_node,
                                                        const size_t& stream_source,
                                                        const size_t& current_colour,
                                                        const int& stream,
                                                        const ORFVector& ORF_vector,
                                                        const size_t& max_ORF_path_length,
                                                        std::unordered_set<int>& prev_node_set,
                                                        const bool is_ref)
{
    // initialise return vector of upstream ORFs (vectors will break at branches in graph)
    std::vector<std::pair<size_t, size_t>> connected_ORFs;

    // get the colour array of the head node
    auto colour_arr = graph_vector.at(abs(head_node) - 1).full_colour();

    // generate path list, vector for path and the stack
    ORFStack ORF_stack;

    // create node set for identification of repeats
    std::unordered_set<int> node_set;
    node_set.insert(head_node * stream);

    // create first item in stack, multiply by stream - upstream (stream = -1) or downstream (stream = 1) and add stream ORF
    ORF_stack.push(std::make_tuple(head_node * stream, colour_arr, 0));

    while(!ORF_stack.empty())
    {
        // pop node in stack
        auto path_tuple = ORF_stack.top();
        const auto& node_id = std::get<0>(path_tuple);
        colour_arr = std::get<1>(path_tuple);
        auto& path_length = std::get<2>(path_tuple);
        ORF_stack.pop();

        // get unitig_dict entry in graph_vector
        const auto& node_dict = graph_vector.at(abs(node_id) - 1);

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // iterate over neighbours
        for (const auto& neighbour : node_dict.get_neighbours(strand))
        {
            // make copy of path_length
            auto temp_path_length = path_length;

            // parse neighbour information.
            const auto& neighbour_id = std::get<0>(neighbour);
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
            if (is_in)
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

            // add node to temp_node_set
            node_set.insert(neighbour_id);

            // check if node is traversed by end of an ORF
            if (!neighbour_dict.ORFs_empty(current_colour))
            {
                // pull out next ORFs and order them
                const auto& next_ORFs = neighbour_dict.get_ORFs(current_colour);
                const auto ordered_ORFs = order_ORFs_in_node(graph_vector, next_ORFs, neighbour_id, ORF_vector);

                // add first entry and stream_source
                // add ORFs so that lowest ORF_ID is added first (to remove redundant edges)
                std::pair<size_t, size_t> ORF_pair;
                ORF_pair.first = (stream_source < ordered_ORFs.at(0)) ? stream_source : ordered_ORFs.at(0);
                ORF_pair.second = (stream_source < ordered_ORFs.at(0)) ? ordered_ORFs.at(0) : stream_source;
                connected_ORFs.push_back(std::move(ORF_pair));

                // pair all entries
                if (ordered_ORFs.size() > 1)
                {
                    for (int i = 0; i < ordered_ORFs.size() - 1; i++)
                    {
                        ORF_pair.first = (ordered_ORFs.at(i) < ordered_ORFs.at(i + 1)) ? ordered_ORFs.at(i) : ordered_ORFs.at(i + 1);
                        ORF_pair.second = (ordered_ORFs.at(i) < ordered_ORFs.at(i + 1)) ? ordered_ORFs.at(i + 1) : ordered_ORFs.at(i);
                        connected_ORFs.push_back(std::move(ORF_pair));;
                    }
                }
            }
            else
            {
                // add to stack if max_ORF_path_length not violated. If next node violates, still traverse to see if meet ORF
                if (temp_path_length <= max_ORF_path_length)
                {
                    temp_path_length += neighbour_dict.size().second;
                    ORF_stack.push({neighbour_id, updated_colours_arr, temp_path_length});
                }
            }
        }
    }
    // update previous nodes encountered
    prev_node_set.insert(std::make_move_iterator(node_set.begin()), std::make_move_iterator(node_set.end()));
    return connected_ORFs;
}

