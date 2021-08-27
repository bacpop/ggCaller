#include "ORF_connection.h"

//PathOverlapMap overlapping_paths(const GraphVector& graph_vector,
//                                 const PathVector& all_paths)
//{
//    // initialise overlap map for path. This describes how overlapping path sits relative to current path
//    // e.g. in case a, the overlapping path sits downstream to current path in forward orientation
//    // a == downstream, forward
//    // b == upstream, forward
//    // c == downstream, reverse
//    // d == upstream, reverse
//    PathOverlapMap path_overlap_map;
//
//    // intialise Eigen Triplet
//    std::vector<ET> tripletList;
//
//    // iterate over each ORF sequence with specific colours combination
//    for (size_t path_index = 0; path_index < all_paths.size(); path_index++)
//    {
//        // iterate over nodes traversed by ORF
//        const auto& path_nodes = all_paths.at(path_index);
//        for (const auto& node_traversed : path_nodes)
//        {
//            // add to triplet list, with temp_ORF_ID (row), node id (column) and set value as 1 (true)
//            // convert node_traversed to size_t, minus 1 as unitigs are one-based, needs to be zero based
//            const size_t abs_node_id = abs(node_traversed) - 1;
//            tripletList.push_back(ET(path_index, abs_node_id, 1));
//        }
//    }
//
//    // initialise sparse matrix
//    Eigen::SparseMatrix<double> mat(all_paths.size(), graph_vector.size());
//    mat.setFromTriplets(tripletList.begin(), tripletList.end());
//
//    // conduct transposition + matrix multiplication to calculate ORFs sharing nodes
//    auto path_overlap_mat = ((mat * mat.transpose()).pruned()).eval();
//
//    // iterate over non-zero entries in matrix and calculate overlaps
//    for (int outit = 0; outit < path_overlap_mat.outerSize(); ++outit)
//    {
//        for (Eigen::SparseMatrix<double>::InnerIterator init(path_overlap_mat, outit); init; ++init)
//        {
//            // iterate through the bottom half of the symmetric matrix and ignore line of symmetry, pass if row <= col
//            if(init.row() <= init.col())
//            {
//                continue;
//            }
//
//            // assign path1 as the longer of the two paths
//            const bool col_longer = ((all_paths.at(init.col()).size() >= all_paths.at(init.row()).size()) ? true : false);
//            const size_t& path1_ID = (col_longer ? init.col() : init.row());
//            const size_t& path2_ID = (col_longer ? init.row() : init.col());
//
//            // emplace in path_overlap_map
//            path_overlap_map.try_emplace(path1_ID, 2);
//            path_overlap_map.try_emplace(path2_ID, 2);
//
//            // work out how paths are orientated
//            const auto& path1 = all_paths.at(path1_ID);
//            auto path2 = all_paths.at(path2_ID);
//
//            auto& path2_start = path2.at(0);
//            auto& path2_end = path2.back();
//
//            // check if start of path2 is present in forward traversal
//            auto it_start = std::find(path1.begin(), path1.end(), path2_start);
//            auto it_end = std::find(path1.begin(), path1.end(), path2_end);
//            bool start_overlap = true;
//            bool end_overlap = true;
//            bool reversed = false;
//            if (it_start == path1.end())
//            {
//                // check if end of path1 is present in reverse traversal
//                it_start = std::find(path1.begin(), path1.end(), path2_start * -1);
//                if (it_start == path1.end())
//                {
//                    start_overlap = false;
//                } else
//                {
//                    reversed = true;
//                }
//            }
//            if (it_end == path1.end())
//            {
//                // check if end of path2 is present in reverse traversal
//                it_end = std::find(path1.begin(), path1.end(), path2_end * -1);
//                // if still not found, paths do not overlap
//                if (it_end == path1.end())
//                {
//                    end_overlap = false;
//                } else
//                {
//                    reversed = true;
//                }
//            }
//
//            // check if neither start or end found, if so then pass
//            if (!start_overlap && !end_overlap)
//            {
//                continue;
//            }
//
//            // if reversed is true, iterate through path2 coordinates and reverse
//            // check if both strands are the same. If not, reverse the nodes of ORF2 and their within-node coordinates
//            if (reversed)
//            {
//                // Reverse ORF2 node vector
//                std::reverse(path2.begin(), path2.end());
//
//                // reverse sign on each ID in ORF2_nodes
//                for (auto & node_id : path2)
//                {
//                    node_id = node_id * -1;
//                }
//
//                // reverse start and end overlap to correct for reversal
//                bool temp_start = start_overlap;
//                start_overlap = end_overlap;
//                end_overlap = temp_start;
//            }
//
//            // boolean to determine if check for path2 within path1 is correct
//            bool overlap_complete = false;
//
//            // check if path2 sits within path1
//            if (start_overlap && end_overlap)
//            {
//                // initialise start and end index lists for ORF2 within ORF1
//                std::vector<size_t> start_index_list;
//                auto iter = find(path1.begin(), path1.end(), path2_start);
//                while(iter != path1.end())
//                {
//                    size_t index = iter - path1.begin();
//                    start_index_list.push_back(std::move(index));
//                    iter = find(iter + 1, path1.end(), path2_start);
//                }
//                std::vector<size_t> end_index_list;
//                iter = find(path1.begin(), path1.end(), path2_end);
//                while(iter != path1.end())
//                {
//                    size_t index = iter - path1.begin();
//                    end_index_list.push_back(std::move(index));
//                    iter = find(iter + 1, path1.end(), path2_end);
//                }
//                for (const auto& start_index : start_index_list)
//                {
//                    for (const auto& end_index : end_index_list)
//                    {
//                        // check moving in forward direction
//                        std::vector<int> path1_nodes_sliced;
//                        if (start_index <= end_index)
//                        {
//                            // slice ORF1 node vector from first entry to last
//                            path1_nodes_sliced = std::vector<int> (path1.begin() + start_index, path1.begin() + end_index + 1);
//                        }
//                        // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
//                        if (path1_nodes_sliced == path2)
//                        {
//                            // think about reversal here! We need to know how both paths are orientated relative to eachother and store this info!
//                            if (!reversed)
//                            {
//                                // path1 sits upstream of path2 in forward orientation, overlapping from index 0
//                                path_overlap_map[path2_ID][0].push_back({(path1_ID + 1) * -1, 0, path2.size() - 1});
//                                // path2 sits downstream of path1 in forward orientation, overlapping from start_index
//                                path_overlap_map[path1_ID][1].push_back({(path2_ID + 1), start_index, end_index});
//                                int test = 0;
//                            } else
//                            {
//                                // path1 sits upstream of path2 in reverse orientation, overlapping from start to end
//                                path_overlap_map[path2_ID][0].push_back({(path1_ID + 1), 0, path2.size() - 1});
//                                //path2 sits downstream of path1 in reverse orientation, overlapping from start_index to end_index
//                                path_overlap_map[path1_ID][1].push_back({(path2_ID + 1) * -1, start_index, end_index});
//                                int test = 0;
//                            }
//
//                            overlap_complete = true;
//                            break;
//                        }
//                    }
//                    // if chain is complete, don't iterate over any more start codons
//                    if (overlap_complete)
//                    {
//                        break;
//                    }
//                }
//            }
//
//            // may be case that path2 not within path1 even if start and end are found, so go and check start and end of path2
//            if (!overlap_complete)
//            {
//                // check start of path2 overlaps with end of path1
//                if (start_overlap)
//                {
//                    std::vector<size_t> start_index_list;
//                    auto iter = find(path1.begin(), path1.end(), path2_start);
//                    while(iter != path1.end())
//                    {
//                        size_t index = iter - path1.begin();
//                        start_index_list.push_back(std::move(index));
//                        iter = find(iter + 1, path1.end(), path2_start);
//                    }
//                    for (const auto& start_index : start_index_list)
//                    {
//                        // check moving in forward direction
//                        std::vector<int> path1_nodes_sliced;
//                        std::vector<int> path2_nodes_sliced;
//
//                        // slice ORF1 node vector from first entry to last
//                        path1_nodes_sliced = std::vector<int> (path1.begin() + start_index, path1.end());
//                        // check if ORF1_slice is too large to slice ORF1 (meaning likely ORF2 is reversed
//                        if (path1_nodes_sliced.size() <= path2.size())
//                        {
//                            // slice ORF1 node vector from first entry to the equivalent length of path2_nodes_sliced
//                            path2_nodes_sliced = std::vector<int>(path2.begin(), path2.begin() + path1_nodes_sliced.size());
//                        }
//
//                        // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
//                        if (path1_nodes_sliced == path2_nodes_sliced)
//                        {
//                            if (!reversed)
//                            {
//                                // path1 sits upstream of path2 in forward orientation, overlapping from start
//                                path_overlap_map[path2_ID][0].push_back({(path1_ID + 1) * -1, 0, path2_nodes_sliced.size() - 1});
//                                //path2 sits downstream of path1 in forward orientation, overlapping from start_index to end
//                                path_overlap_map[path1_ID][1].push_back({(path2_ID + 1), start_index, path1.size() - 1});
//                                int test = 0;
//                            } else
//                            {
//                                // path1 sits downstream of path2 in reverse orientation, overlapping from start of overlap to end
//                                path_overlap_map[path2_ID][1].push_back({(path1_ID + 1) * -1, path2.size() - path2_nodes_sliced.size(), path2.size() - 1});
//                                // path2 sits downstream of path1 in reverse orientation, overlapping from start_index to end
//                                path_overlap_map[path1_ID][1].push_back({(path2_ID + 1) * -1, start_index, path1.size() - 1});
//                                int test = 0;
//                            }
//                            break;
//                        }
//                    }
//                }
//                // check end of path2 overlaps with start of path1
//                else
//                {
//                    std::vector<size_t> end_index_list;
//                    auto iter = find(path1.begin(), path1.end(), path2_end);
//                    while(iter != path1.end())
//                    {
//                        size_t index = iter - path1.begin();
//                        end_index_list.push_back(std::move(index));
//                        iter = find(iter + 1, path1.end(), path2_end);
//                    }
//                    for (const auto& end_index : end_index_list)
//                    {
//                        // check moving in forward direction
//                        std::vector<int> path1_nodes_sliced;
//                        std::vector<int> path2_nodes_sliced;
//
//                        // slice ORF1 node vector from first entry to last
//                        path1_nodes_sliced = std::vector<int> (path1.begin(), path1.begin() + end_index + 1);
//                        // check that slice is small enough to correctly slice ORF2_nodes
//                        if (path1_nodes_sliced.size() <= path2.size())
//                        {
//                            // slice ORF2 node vector from last entry - end_index to the equivalent length of path2_nodes_sliced
//                            path2_nodes_sliced = std::vector<int>(path2.end() - end_index - 1, path2.end());
//                        }
//
//                        // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
//                        if (path1_nodes_sliced == path2_nodes_sliced)
//                        {
//                            if (!reversed)
//                            {
//                                // path1 sits downstream of path2 in forward orientation, overlapping from start of overlap to end
//                                path_overlap_map[path2_ID][1].push_back({(path1_ID + 1), path2.size() - path2_nodes_sliced.size(), path2.size() - 1});
//                                // path2 sits upstream of path1 in forward orientation, overlapping from start to end index
//                                path_overlap_map[path1_ID][0].push_back({(path2_ID + 1) * -1, 0, end_index});
//                                int test = 0;
//                            } else
//                            {
//                                // path1 sits upstream of path2 in reverse orientation, overlapping from start to the end of overlap
//                                path_overlap_map[path2_ID][0].push_back({(path1_ID + 1), 0, path2_nodes_sliced.size() - 1});
//                                //path2 sits upstream of path1 in reverse orientation, overlapping from start to end_index
//                                path_overlap_map[path1_ID][0].push_back({(path2_ID + 1), 0, end_index});
//                                int test = 0;
//                            }
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    return path_overlap_map;
//}
//
//// traverse path_overlap_map to find next ORF upstream/downstream
//std::vector<std::pair<size_t, size_t>> pair_ORF_paths (const PathOverlapMap& path_overlap_map,
//                                                       const std::unordered_map<size_t, std::vector<size_t>>& ORF_path_map,
//                                                       const ORFVector& ORF_vector,
//                                                       const std::unordered_set<size_t>& target_ORFs,
//                                                       std::unordered_set<size_t>& unpaired_downstream,
//                                                       std::unordered_set<size_t> unpaired_upstream) {
//    // initilise vector to hold results
//    std::vector<std::pair<size_t, size_t>> connected_ORFs;
//
//    // iterate over target_ORFs
//    for (const auto &ORF_ID : target_ORFs) {
//        // check if ORF has already been paired upstream/downstream
//        if (unpaired_downstream.find(ORF_ID) == unpaired_downstream.end() && unpaired_upstream.find(ORF_ID) == unpaired_upstream.end())
//        {
//            continue;
//        }
//
////        // boolean to check if upstream and downstream paired
////        bool upstream_paired = false;
////        bool downstream_paired = false;
//
//        // find if there are any other ORFs in the same path
//        const auto &ORF_info = ORF_vector.at(ORF_ID);
//        const auto &ORF_path = std::get<6>(ORF_info);
//
//        auto ORF_list = ORF_path_map.at(ORF_path);
//        // if other ORFs traverse the path, order them
//        if (ORF_list.size() > 1)
//        {
//            order_ORFs_in_path(ORF_list, ORF_vector);
//            // work out if any of the other ORFs are target ORFs
//            for (int i = 0; i < ORF_list.size() - 1; i++)
//            {
//                connected_ORFs.push_back({ORF_list.at(i), ORF_list.at(i + 1)});
//                // remove all but first entry from unpaired_upstream
//                if (i != 0)
//                {
//                    auto it = unpaired_upstream.find(ORF_list.at(i));
//                    if (it != unpaired_upstream.end())
//                    {
//                        unpaired_upstream.erase(it);
//                    }
//                }
//                // do the same for the last entry as not reached
//                if (i == ORF_list.size() - 2)
//                {
//                    auto it = unpaired_upstream.find(ORF_list.at(i + 1));
//                    if (it != unpaired_upstream.end())
//                    {
//                        unpaired_upstream.erase(it);
//                    }
//                }
//
//                // remove all but last entry from unpaired_downstream (not iterating over so can ignore)
//                auto it = unpaired_downstream.find(ORF_list.at(i));
//                if (it != unpaired_downstream.end())
//                {
//                    unpaired_downstream.erase(it);
//                }
//            }
//            // traverse path_overlap_map upstream
//            if (unpaired_upstream.find(ORF_ID) == unpaired_upstream.end())
//            {
//                auto stream_ORFs = std::move(check_next_path(path_overlap_map, ORF_path, ORF_path_map, ORF_vector, -1));
//
//                // remove paired ORFs from unpaired_downstream
//                if (!stream_ORFs.empty())
//                {
//                    for (const auto& paired_ORF_ID : stream_ORFs)
//                    {
//                        connected_ORFs.push_back({paired_ORF_ID, ORF_ID});
//                        auto it = unpaired_downstream.find(paired_ORF_ID);
//                        if (it != unpaired_downstream.end())
//                        {
//                            unpaired_downstream.erase(it);
//                        }
//                    }
//
//                    // remove current ORF from unpaired_upstream
//                    auto it = unpaired_upstream.find(ORF_ID);
//                    if (it != unpaired_upstream.end())
//                    {
//                        unpaired_downstream.erase(it);
//                    }
//                }
//            }
//            // repeat for downstream
//            // traverse path_overlap_map upstream
//            if (unpaired_downstream.find(ORF_ID) == unpaired_downstream.end())
//            {
//                auto stream_ORFs = std::move(check_next_path(path_overlap_map, ORF_path, ORF_path_map, ORF_vector, 1));
//
//                // remove paired ORFs from unpaired_downstream
//                if (!stream_ORFs.empty())
//                {
//                    for (const auto& paired_ORF_ID : stream_ORFs)
//                    {
//                        connected_ORFs.push_back({ORF_ID, paired_ORF_ID});
//                        auto it = unpaired_upstream.find(paired_ORF_ID);
//                        if (it != unpaired_upstream.end())
//                        {
//                            unpaired_upstream.erase(it);
//                        }
//                    }
//
//                    // remove current ORF from unpaired_upstream
//                    auto it = unpaired_downstream.find(ORF_ID);
//                    if (it != unpaired_downstream.end())
//                    {
//                        unpaired_downstream.erase(it);
//                    }
//                }
//            }
//        }
//    }
//    return connected_ORFs;
//}
//
//void order_ORFs_in_path(std::vector<size_t>& ORF_list,
//                        const ORFVector& ORF_vector)
//{
//    std::vector<indexPair> ordered_ORFs;
//
//    // iterate over ORF_list and get the index of the ORF in the path
//    for (const auto& ORF_ID : ORF_list)
//    {
//        indexPair new_pair;
//
//        const auto& ORF_info = ORF_vector.at(ORF_ID);
////        const auto& ORF_path = std::get<6>(ORF_info);
//        const auto& path_indices = std::get<7>(ORF_info);
//        // need to find index of path in ORF_info
////        auto found = std::find(ORF_paths.begin(), ORF_paths.end(), path_ID);
////        size_t index = found - ORF_paths.begin();
//
//        // add to ordered_ORFs for sorting
//        new_pair.first = path_indices.first;
//        new_pair.second = ORF_ID;
//
//        ordered_ORFs.push_back(std::move(new_pair));
//    }
//
//    // sort based on first entry in ordered_ORFs
//    sort(ordered_ORFs.begin(), ordered_ORFs.end());
//
//    // insert into ORF_list in place
//    ORF_list.clear();
//    std::transform(ordered_ORFs.begin(), ordered_ORFs.end(), std::back_inserter(ORF_list),
//                   [] (auto const& pair) {return pair.second; });
//}
//
//std::vector<size_t> check_next_path (const PathOverlapMap& path_overlap_map,
//                                    const size_t& head_node,
//                                    const std::unordered_map<size_t, std::vector<size_t>>& ORF_path_map,
//                                    const ORFVector& ORF_vector,
//                                    const int& stream)
//{
//    // initialise vector to return
//    std::vector<size_t> stream_ORFs;
//
//    // generate stack of paths to traverse
//    std::stack<int> path_stack;
//
//    // convert head_node id to path_ID
//    int path_ID = head_node + 1 * stream;
//
//    // create path set for identification of repeats
//    std::unordered_set<int> path_set;
//    path_set.insert(path_ID);
//
//    // create first item in stack, multiply by stream - upstream (stream = -1) or downstream (stream = 1) and add stream ORF
//    path_stack.push(path_ID);
//
//    while(!path_stack.empty())
//    {
//        // pop node in stack
//        auto path_ID = path_stack.top();
//        path_stack.pop();
//
//        // get unitig_dict entry in graph_vector
//        const auto& path_dict = path_overlap_map.at(abs(path_ID) - 1);
//
//        // determine strand of unitig
//        const bool strand = (path_ID >= 0) ? true : false;
//
//        // iterate over neighbours
//        for (const auto& neighbour : path_dict.at(strand))
//        {
//            const auto& neighbour_id = std::get<0>(neighbour);
//
//            // check if path already iterated over to avoid cycles
//            if (path_set.find(neighbour_id) != path_set.end())
//            {
//                continue;
//            }
//
//            // add node to path_set
//            path_set.insert(neighbour_id);
//
//            // check if there are target ORFs in this path
//            if (ORF_path_map.find(abs(neighbour_id) - 1) != ORF_path_map.end())
//            {
//                auto ORF_list = ORF_path_map.at(abs(neighbour_id) - 1);
//                // if other ORFs traverse the path, order them
//                if (ORF_list.size() > 1)
//                {
//                    order_ORFs_in_path(ORF_list, ORF_vector);
//
//                    // if path is traversed in reverse, need to reverse the ORF_list ordering
//                    if (neighbour_id < 0)
//                    {
//                        std::reverse(ORF_list.begin(), ORF_list.end());
//                    }
//                }
//
//                // if upstream, need to add last ORF in list, if downstream add first
//                if (stream == -1)
//                {
//                    stream_ORFs.push_back(ORF_list.back());
//                } else
//                {
//                    stream_ORFs.push_back(ORF_list.at(0));
//                }
//            }
//            else
//            {
//                // if no ORFs found, continue traversal
//                path_stack.push(neighbour_id);
//            }
//        }
//    }
//    return stream_ORFs;
//}

void add_ORF_info (GraphVector& graph_vector,
                  const size_t& colour_ID,
                  const std::unordered_set<size_t>& target_ORFs,
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

//        // check if source and sink ORFs are the same, if not continue
//        if (source != sink)
//        {
//            // get graph information for source node
//            const auto & ORF_info = ORF_vector.at(sink);
//
//            // get start and end node IDs (-1 as graph is zero-based)
//            const auto source_node_id = abs(std::get<0>(ORF_info).at(0)) - 1;
//            const auto sink_node_id = abs(std::get<0>(ORF_info).back()) - 1;
//
//            // add ORF information to graph
//            graph_vector.at(source_node_id).set_ORFs(colour_ID, sink);
//            // check if ORF is present only on single node
//            if (source_node_id != sink_node_id)
//            {
//                graph_vector.at(sink_node_id).set_ORFs(colour_ID, sink);
//            }
//        }
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
                                                      std::unordered_set<int>& prev_node_set)
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
//            if (sink_node_ORFs.size() >= 2)
//            {
//                // if 2 or more ORFs overlapping, order and then add each connecting edge
//                const auto ordered_ORFs = order_ORFs_in_node(graph_vector, sink_node_ORFs, sink_node_id, ORF_vector);
//                for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
//                {
//                    ORF_edges.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
//                }
//            }

            // now traverse graph using DFS, finding next ORFs upstream and downstream of source and sink node respectively
            // check that sink hasn't been traversed in reverse (downstream) or source hasn't been traverse forward (downstream

//            // if going downstream (stream > 0), check that ORF is last, or upstream (stream < 0), check it is first
            if ((stream > 0 && source == ordered_ORFs.back()) || stream < 0 && source == ordered_ORFs.at(0))
            {
                auto next_ORFs = check_next_ORFs(graph_vector, start_node, source, colour_ID, stream, ORF_vector, max_ORF_path_length, prev_node_set);
                ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
            }

            prev_node_set.insert(start_node);
        }

//            // traverse downstream. Check that sink node hasn't been traversed in reverse, otherwise ORF has been paired downstream already.
//            if (prev_node_set.find(sink_node_id * -1) == prev_node_set.end())
//            {
//                auto next_ORFs = check_next_ORFs(graph_vector, sink_node_id, source, colour_ID, 1, ORF_vector, max_ORF_path_length, prev_node_set);
//                ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
//            }
//        // scope for second item in end_ORF pair
//        if (end_ORF_pair.first != end_ORF_pair.second)
//        {
//            const auto& start_ORF = end_ORF_pair.second;
//
//            // get ORF info
//            const auto& start_ORF_info = ORF_vector.at(start_ORF);
//
//            // get id for source_node and sink node for ORF
//            const auto& source_node_id = std::get<0>(start_ORF_info).at(0);
//            const auto& sink_node_id = std::get<0>(start_ORF_info).back();
//
//            // get source_node and sink_node info
//            const auto& source_node_info = graph_vector.at(abs(source_node_id) - 1);
//            const auto& sink_node_info = graph_vector.at(abs(sink_node_id) - 1);
//
//            // check if there are any overlapping ORFs on current node and order...
//            const auto& source_node_ORFs = source_node_info.get_ORFs(colour_ID);
//            const auto& sink_node_ORFs = sink_node_info.get_ORFs(colour_ID);
//
//
//            // check if there are overlapping ORFs on the source and sink node, as these will not be picked up
//            if (source_node_ORFs.size() >= 2)
//            {
//                // if 2 or more ORFs overlapping, order and then add each connecting edge
//                const auto ordered_ORFs = order_ORFs_in_node(graph_vector, source_node_ORFs, source_node_id, ORF_vector);
//                for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
//                {
//                    ORF_edges.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
//                }
//            }
//            if (sink_node_ORFs.size() >= 2)
//            {
//                // if 2 or more ORFs overlapping, order and then add each connecting edge
//                const auto ordered_ORFs = order_ORFs_in_node(graph_vector, sink_node_ORFs, sink_node_id, ORF_vector);
//                for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
//                {
//                    ORF_edges.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
//                }
//            }
//
//            // now traverse graph using DFS, finding next ORFs upstream and downstream of source and sink node respectively
//            // traverse upstream. Check that source node hasn't been traversed forward, otherwise ORF has been paired downstream already.
//            if (prev_node_set.find(source_node_id) == prev_node_set.end())
//            {
//                auto next_ORFs = check_next_ORFs(graph_vector, source_node_id, start_ORF, colour_ID, -1, ORF_vector, max_ORF_path_length, prev_node_set);
//                ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
//            }
//            // traverse downstream. Check that sink node hasn't been traversed in reverse, otherwise ORF has been paired downstream already.
//            if (prev_node_set.find(sink_node_id * -1) == prev_node_set.end())
//            {
//                auto next_ORFs = check_next_ORFs(graph_vector, sink_node_id, start_ORF, colour_ID, 1, ORF_vector, max_ORF_path_length, prev_node_set);
//                ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
//            }
//        }
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
                                                        std::unordered_set<int>& prev_node_set)
{
    // initialise return vector of upstream ORFs (vectors will break at branches in graph)
    std::vector<std::pair<size_t, size_t>> connected_ORFs;

    // set last traversed ORF
    size_t last_ORF = stream_source;

    // set expected next node to 0 (as 0 cannot be used as node descriptor as 1 indexed
    std::tuple<int, size_t, bool> expected_node = {0, 0, true};

    // get the colour array of the head node
    auto colour_arr = graph_vector.at(abs(head_node) - 1).full_colour();

    // generate path list, vector for path and the stack
    ORFStack ORF_stack;

    // create node set for identification of repeats
    std::unordered_set<int> node_set;
    node_set.insert(head_node * stream);

    // create first item in stack, multiply by stream - upstream (stream = -1) or downstream (stream = 1) and add stream ORF
    //ORF_stack.push(std::make_tuple(head_node * stream, last_ORF, expected_node, colour_arr, 0));
    ORF_stack.push(std::make_tuple(head_node * stream, colour_arr, 0));

    while(!ORF_stack.empty())
    {
        // pop node in stack
        auto path_tuple = ORF_stack.top();
        const auto& node_id = std::get<0>(path_tuple);
//        last_ORF = std::get<1>(path_tuple);
//        expected_node = std::get<2>(path_tuple);
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
            // make copy of last_ORF and expected_node
//            auto temp_last_ORF = last_ORF;
//            auto temp_expected_node = expected_node;
            auto temp_path_length = path_length;

            // parse neighbour information.
            const auto& neighbour_id = neighbour.first;

//            // check if unitig has already been traversed, and pass if repeat not specified
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

            // expected node not met, reset the temp_last_ORF to the stream_source
//            if (std::get<0>(temp_expected_node) != 0 && std::get<0>(temp_expected_node) != neighbour_id)
//            {
//                temp_last_ORF = stream_source;
//                temp_expected_node = {0, 0, true};
//            }

            // add node to temp_node_set
            node_set.insert(neighbour_id);

            // check if node is traversed by end of an ORF
            if (!neighbour_dict.ORFs_empty(current_colour))
            {
                // pull out next ORFs and order them
                const auto& next_ORFs = neighbour_dict.get_ORFs(current_colour);
                const auto ordered_ORFs = order_ORFs_in_node(graph_vector, next_ORFs, neighbour_id, ORF_vector);

                // add first entry and stream_source
                connected_ORFs.push_back({stream_source, ordered_ORFs.at(0)});

                // pair all entries
                if (ordered_ORFs.size() > 1)
                {
                    for (int i = 0; i < ordered_ORFs.size() - 1; i++)
                    {
                        connected_ORFs.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
                    }
                }
            }



                // bool for determining if paired ORFs are the same, and if they are already connected as end_ORFs
//                bool end_found = false;

//                // make sure that 5' and 3' of ORF is traversed if it is not uninode
//                if (temp_last_ORF == stream_source)
//                {
//                    // iterate over ordered_ORFs, if encounter uninode ORF then add and break. If not, then use final
//                    for (size_t i = 0; i < ordered_ORFs.size(); i++)
//                    {
//                        if (ordered_ORFs.at(i) != stream_source)
//                        {
//                            if (std::get<0>(ORF_vector.at(ordered_ORFs.at(i))).size() == 1)
//                            {
//                                connected_ORFs.push_back({stream_source, ordered_ORFs.at(i)});
//                                end_found = true;
//                                break;
//                            }
//                        }
//                    }
//                }
//                // if temp_last_ORF is not stream_source, and is present, means that both ends now encountered, so fully traversed
//                else if (std::find(ordered_ORFs.begin(), ordered_ORFs.end(), temp_last_ORF) != ordered_ORFs.end())
//                {
//                    connected_ORFs.push_back({stream_source, temp_last_ORF});
//                    end_found = true;
//                }

                // if end not found, then need to continue to find the next completely traversed ORF
//                if (!end_found)
//                {
//                    // need to assign next expected ORF
//                    bool forward = true;
//                    const auto& next_ORF_info = ORF_vector.at(ordered_ORFs.back());
//                    const auto& node_order = std::get<0>(next_ORF_info);
//                    auto it = std::find(node_order.begin(), node_order.end(), neighbour_id);
//
//                    if (it == node_order.end())
//                    {
//                        it = std::find(node_order.begin(), node_order.end(), neighbour_id * -1);
//                        forward = false;
//                    }
//
//                    // generate items for next iteration
//                    auto index = std::distance(node_order.begin(), it);
//                    int next_node;
//                    size_t next_ORF;
//
//                    // may be case that last entry is next_node but node is positive. If the node isn't start or end as expected, need to ignore.
//                    // (if end, should be negative, if front should be positive, else assume stream_source is next upstream ORF)
//                    if (forward && index == 0)
//                    {
//                        next_node = node_order.at(++index);
//                        next_ORF = ordered_ORFs.back();
//                    } else if (!forward && index == node_order.size() - 1)
//                    {
//                        next_node = node_order.at(--index) * -1;
//                        next_ORF = ordered_ORFs.back();
//                    } else
//                    {
//                        next_node = 0;
//                        next_ORF = stream_source;
//                    }

                    // add to stack if max_ORF_path_length not violated. If next node violates, still traverse to see if meet ORF
//                    if (temp_path_length <= max_ORF_path_length)
//                    {
//                        temp_path_length += neighbour_dict.size().second;
//                        ORF_stack.push({neighbour_id, next_ORF, {next_node, index, forward}, updated_colours_arr, temp_path_length});
//                    }
            else
            {
//                if (temp_last_ORF != stream_source)
//                {
//                    // calculate next expected node
//                    const auto& next_ORF_info = ORF_vector.at(temp_last_ORF);
//                    const auto& node_order = std::get<0>(next_ORF_info);
//
//                    auto index = std::get<1>(temp_expected_node);
//                    auto forward = std::get<2>(temp_expected_node);
//
//                    int next_node;
//
//                    if (forward)
//                    {
//                        next_node = node_order.at(++index);
//                    } else
//                    {
//                        next_node = node_order.at(--index) * -1;
//                    }
//
//                    temp_expected_node = {next_node, index, forward};
//                }
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

