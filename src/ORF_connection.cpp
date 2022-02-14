#include "ORF_connection.h"
robin_hood::unordered_map<size_t, robin_hood::unordered_set<size_t>> add_ORF_info (const std::vector<Kmer>& head_kmer_arr,
                                                                                     const robin_hood::unordered_set<size_t>& target_ORFs,
                                                                                     const ORFNodeMap& gene_map)
{
    robin_hood::unordered_map<size_t, robin_hood::unordered_set<size_t>> node_to_ORFs;

    for (const auto & source : target_ORFs)
    {
        {
            // get graph information for source node
            const auto & ORF_info = gene_map.at(source);

            const auto source_node_id = std::get<0>(ORF_info).at(0);
            const auto sink_node_id = std::get<0>(ORF_info).back();

            // add the ORF ID to node_to_ORFs
            node_to_ORFs[abs(source_node_id) - 1].insert(source);
            node_to_ORFs[abs(sink_node_id) - 1].insert(source);
        }
    }

    return node_to_ORFs;
}

std::vector<size_t> order_ORFs_in_node(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                       const std::vector<Kmer>& head_kmer_arr,
                                       const robin_hood::unordered_set<size_t>& node_ORFs,
                                       const int node_id,
                                       const ORFNodeMap& gene_map)
{
    std::vector<size_t> overlapping_ORF_IDs;
    std::vector<std::pair<bool, indexPair>> overlapping_ORF_coords;
    // iterate over ORFs in node_ORFs
    for (const auto& ORF_ID : node_ORFs)
    {
        // add to return vector overlapping ORF vector
        overlapping_ORF_IDs.push_back(ORF_ID);

        // get reference to gene_map entry
        const auto & ORF_info = gene_map.at(ORF_ID);

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
        // get a reference to the unitig map object
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        auto& um = um_pair.first;

        size_t node_end = um.size - 1;

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

std::set<std::pair<size_t, size_t>> pair_ORF_nodes (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                    const std::vector<Kmer>& head_kmer_arr,
                                                    const robin_hood::unordered_map<size_t, robin_hood::unordered_set<size_t>>& node_to_ORFs,
                                                    const size_t colour_ID,
                                                    const robin_hood::unordered_set<size_t>& target_ORFs,
                                                    const ORFNodeMap& gene_map,
                                                    const size_t& max_ORF_path_length,
                                                    const int stream,
                                                    std::unordered_set<int>& prev_node_set,
                                                    const int overlap,
                                                    const bool is_ref,
                                                    const fm_index_coll& fm_idx)
{
    // initialise pair of vectors (first = upstream of start_ORF, second = downstream of start_ORF)
    std::set<std::pair<size_t, size_t>> ORF_edges;

    // iterate over each entry in end_ORFs
    for (const auto& source : target_ORFs)
    {
        // get ORF info
        const auto & source_info = gene_map.at(source);

        // get id for source_node and sink node for ORF
        auto source_node_id = std::get<0>(source_info).at(0);
        auto sink_node_id = std::get<0>(source_info).back();

        // decide whether to traverse from source (upstream) or sink (downstream)
        int start_node = (stream > 0) ? sink_node_id : source_node_id;

        // check if node has been traversed previously in opposite direction i.e. ORF has already been paired up/downstream
        if (prev_node_set.find(start_node * stream * -1) == prev_node_set.end())
        {
            // get source_node and sink_node info
            // get a reference to the unitig map object
            auto start_um_pair = get_um_data(ccdbg, head_kmer_arr, start_node);
            auto& start_um = start_um_pair.first;
            auto& start_um_data = start_um_pair.second;

            // check if there are any overlapping ORFs on current node and order...
            const auto& start_node_ORFs = node_to_ORFs.at(abs(start_node) - 1);

            // check if there are overlapping ORFs on the source and sink node, as these will not be picked up
            std::vector<size_t> ordered_ORFs = order_ORFs_in_node(ccdbg, head_kmer_arr, start_node_ORFs, start_node, gene_map);
            if (ordered_ORFs.size() >= 2)
            {
                // if 2 or more ORFs overlapping, order and then add each connecting edge
                for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
                {
                    ORF_edges.insert({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
                }
            }

            // now traverse graph using DFS, finding next ORFs upstream and downstream of source and sink node respectively
            // check that sink hasn't been traversed in reverse (downstream) or source hasn't been traverse forward (downstream

            // if going downstream (stream > 0), check that ORF is last, or upstream (stream < 0), check it is first
            if ((stream > 0 && source == ordered_ORFs.back()) || stream < 0 && source == ordered_ORFs.at(0))
            {
                auto next_ORFs = check_next_ORFs(ccdbg, head_kmer_arr, node_to_ORFs, start_node, source, colour_ID, stream, gene_map, max_ORF_path_length, prev_node_set, overlap, is_ref, fm_idx);
                ORF_edges.insert(make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
            }

            prev_node_set.insert(start_node);
        }
    }
    return ORF_edges;
}

std::set<std::pair<size_t, size_t>> check_next_ORFs (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                     const std::vector<Kmer>& head_kmer_arr,
                                                     const robin_hood::unordered_map<size_t, robin_hood::unordered_set<size_t>>& node_to_ORFs,
                                                     const int& head_node,
                                                     const size_t stream_source,
                                                     const size_t current_colour,
                                                     const int stream,
                                                     const ORFNodeMap& gene_map,
                                                     const size_t max_ORF_path_length,
                                                     std::unordered_set<int>& prev_node_set,
                                                     const int overlap,
                                                     const bool is_ref,
                                                     const fm_index_coll& fm_idx)
{
    // initialise return vector of upstream ORFs (vectors will break at branches in graph)
    std::set<std::pair<size_t, size_t>> connected_ORFs;

    // get the colour array of the head node
    auto um_pair = get_um_data(ccdbg, head_kmer_arr, head_node);
    auto& um = um_pair.first;
    auto& um_data = um_pair.second;

    auto colour_arr = um_data->full_colour();

    // generate path list, vector for path and the stack
    ORFStack ORF_stack;
    std::vector<int> node_vector;

    // create node set for identification of repeats
    std::unordered_set<int> node_set;
    node_set.insert(head_node * stream);

    // make vector to search with, one less than ORF_nodes
    std::vector<int> search_vector;
    if (is_ref)
    {
        const auto & source_info = gene_map.at(stream_source);
        const auto & ORF_nodes = std::get<0>(source_info);

        if (stream == -1)
        {
            for (int i = ORF_nodes.size() - 1; i > 0; i--)
            {
                search_vector.push_back(ORF_nodes.at(i) * -1);
            }
        } else
        {
            search_vector.insert(search_vector.end(), ORF_nodes.begin(), ORF_nodes.end() - 1);
        }
    }


    // create first item in stack, multiply by stream - upstream (stream = -1) or downstream (stream = 1) and add stream ORF
    ORF_stack.push(std::make_tuple(0, head_node * stream, colour_arr, 0));

    while(!ORF_stack.empty())
    {
        // pop node in stack
        auto path_tuple = ORF_stack.top();
        ORF_stack.pop();
        const size_t & pos_idx = std::get<0>(path_tuple);
        const auto& node_id = std::get<1>(path_tuple);
        colour_arr = std::get<2>(path_tuple);
        auto& path_length = std::get<3>(path_tuple);

        if (pos_idx != 0)
        {
            node_vector = std::vector<int> (node_vector.begin(), node_vector.begin() + pos_idx);
        }

        // add node to path
        node_vector.push_back(node_id);

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        // get unitig_dict entry in graph_vector
        um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        um = um_pair.first;
        um_data = um_pair.second;

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
            // make copy of path_length
            auto temp_path_length = path_length;

            auto neighbour_da = neighbour_um.getData();
            auto neighbour_um_data = neighbour_da->getData(neighbour_um);
            const bool neighbour_strand = neighbour_um.strand;


            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const int neighbour_id = (neighbour_strand) ? neighbour_um_data->get_id() : neighbour_um_data->get_id() * -1;

            // check if unitig has already been traversed, and pass if repeat not specified
            const bool is_in = node_set.find(neighbour_id) != node_set.end();
            if (is_in)
            {
                continue;
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
            if (!is_ref)
            {
                if (neighbour_colour.count() < 2)
                {
                    continue;
                }
            }

            // add node to temp_node_set
            node_set.insert(neighbour_id);

            // check if node is traversed by end of an ORF
            const auto ORF_found = node_to_ORFs.find(abs(neighbour_id) - 1);
            if (ORF_found != node_to_ORFs.end())
            {
                //pull out next ORFs
                auto next_ORFs = ORF_found->second;
                // check full gene paths against fm-index, if not present then pass
                if (is_ref)
                {
                    std::vector<size_t> to_remove;
                    auto search_vector_temp = search_vector;
                    search_vector_temp.insert(search_vector_temp.end(), node_vector.begin(), node_vector.end());
                    for (const auto & ORF : next_ORFs)
                    {
                        auto search_vector_temp_iter = search_vector_temp;
                        const auto & ORF_info = gene_map.at(ORF);
                        const auto & ORF_nodes = std::get<0>(ORF_info);

                        if (neighbour_id == ORF_nodes.at(0))
                        {
                            search_vector_temp_iter.insert(search_vector_temp_iter.end(), ORF_nodes.begin(), ORF_nodes.end());
                        } else
                        {
                            for (int i = ORF_nodes.size() - 1; i > -1; i--)
                            {
                                search_vector_temp_iter.push_back(ORF_nodes.at(i) * -1);
                            }
                        }

                        const auto present = path_search(search_vector_temp_iter, fm_idx);
                        if (!present.first)
                        {
                            to_remove.push_back(ORF);
                        }
                    }

                    if (next_ORFs.size() == to_remove.size())
                    {
                        continue;
                    } else
                    {
                        for (const auto& ORF : to_remove)
                        {
                            next_ORFs.erase(ORF);
                        }
                    }
                }

                // pull out next ORFs and order them
                const auto ordered_ORFs = order_ORFs_in_node(ccdbg, head_kmer_arr, next_ORFs, neighbour_id, gene_map);

                // add first entry and stream_source
                // add ORFs so that lowest ORF_ID is added first (to remove redundant edges)
                std::pair<size_t, size_t> ORF_pair;
                ORF_pair.first = (stream_source < ordered_ORFs.at(0)) ? stream_source : ordered_ORFs.at(0);
                ORF_pair.second = (stream_source < ordered_ORFs.at(0)) ? ordered_ORFs.at(0) : stream_source;
                connected_ORFs.insert(std::move(ORF_pair));

                // pair all entries
                if (ordered_ORFs.size() > 1)
                {
                    for (int i = 0; i < ordered_ORFs.size() - 1; i++)
                    {
                        ORF_pair.first = (ordered_ORFs.at(i) < ordered_ORFs.at(i + 1)) ? ordered_ORFs.at(i) : ordered_ORFs.at(i + 1);
                        ORF_pair.second = (ordered_ORFs.at(i) < ordered_ORFs.at(i + 1)) ? ordered_ORFs.at(i + 1) : ordered_ORFs.at(i);
                        connected_ORFs.insert(std::move(ORF_pair));;
                    }
                }
            }
            else
            {
                // add to stack if max_ORF_path_length not violated. If next node violates, still traverse to see if meet ORF
                if (temp_path_length <= max_ORF_path_length)
                {
                    temp_path_length += neighbour_um.size - overlap;
                    ORF_stack.push({new_pos_idx, neighbour_id, updated_colours_arr, temp_path_length});
                }
            }
        }
    }
    // update previous nodes encountered
    prev_node_set.insert(std::make_move_iterator(node_set.begin()), std::make_move_iterator(node_set.end()));
    return connected_ORFs;
}

