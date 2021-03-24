// ggCaller header
#include "ggCaller_classes.h"
// define mutex for safe addition to robinhood_maps
//std::mutex mtx2;

std::vector<bool> negate_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2)
{
    std::vector<bool> output_array = array1;
    for (size_t i = 0; i < array1.size(); i++)
    {
        if (array1[i] == 1 && array2[i] == 0)
        {
            output_array[i] = 0;
        }
    }
    return output_array;
}

std::vector<bool> add_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2)
{
    std::vector<bool> output_array = array1;
    for (size_t i = 0; i < array1.size(); i++)
    {
        if (array1[i] == 0 && array2[i] == 1)
        {
            output_array[i] = 1;
        }
    }
    return output_array;
}

PathVector recur_nodes_binary (const unitigMap& graph_map,
                               const robin_hood::unordered_map<std::string, size_t>& head_kmer_map,
                               //const PathMap& previous_paths,
                               const std::vector<std::pair<size_t, bool>>& head_kmer_list,
                               const uint8_t& codon_arr,
                               std::vector<bool> colour_arr,
                               const std::set<std::pair<size_t, bool>>& kmer_set,
                               const size_t& length,
                               const size_t& length_max,
                               const bool& repeat,
                               const vector<bool>& empty_colour_arr)
{
    // generate path list, add head_kmer_list to it
    PathVector path_list;
    std::pair<std::vector<std::pair<size_t, bool>>, std::vector<bool>> path_pair (head_kmer_list, colour_arr);
    path_list.push_back(std::move(path_pair));

    // testing
//    std::set<std::pair<std::string, bool>> node_set;
//    std::vector<std::pair<std::string, bool>> node_vector;
//
//    //testing
//    for (const auto& node : head_kmer_list)
//    {
//        std::string node_str = graph_map.at(node.first).head_kmer;
//        std::pair<std::string, bool> node_pair(node_str, node.second);
//        node_set.insert(node_pair);
//        node_vector.push_back(node_pair);
//    }

//    int test = 0;
//    //std::pair<std::string, bool> head("GCTTGATATCCATTTCATCCCGCCCGGTCAT", true);
//    std::pair<std::string, bool> query1("GCTTGATATCCATTTCATCCCGCCCGGTCAT", false);
//    std::pair<std::string, bool> query2("TATCAAGTGCTTCCTTTGTTTTTAGAATGCG", false);
//    std::pair<std::string, bool> query3("TCATCTTATGCCTTTTTTATCTTCATTTGGG", false);
//    if (node_set.find(query1) != node_set.end() && node_set.find(query2) != node_set.end() && node_set.find(query3) != node_set.end() && colour_arr[45])
//    {
//        test = 1;
//    }

    // if path not complete, generate unitig object using head kmer
    if (codon_arr != 0)
    {
        // get last head-kmer pair in list, generate unitig map object
        const std::pair<size_t, bool>& iter_km_pair = head_kmer_list.back();

        // initilise next_unitigs as a reference
        const auto& next_unitigs = graph_map.at(iter_km_pair.first).neighbours;

        // iterate over neighbours, recurring through incomplete paths
        for (const auto& neighbour : next_unitigs.at(iter_km_pair.second))
        {
            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            size_t kmer_id = std::get<0>(neighbour);
            bool unitig_strand = std::get<1>(neighbour);
            const auto& frame = std::get<2>(neighbour);

            // make path entry for current kmer
            const auto kmer_pair = make_pair(kmer_id, unitig_strand);

            // check if unitig has already been traversed
            const bool is_in = kmer_set.find(kmer_pair) != kmer_set.end();

            // calculate modulus for path and updated cached path reading frame
            int modulus = length % 3;
            std::string cached_kmer_id = std::to_string(kmer_pair.first) + (kmer_pair.second ? "+" : "-");

            //bool cache_valid = 0;
            // don't need this, need to check if arrays AREN'T completed, this means caching can be used to help branching
            // calculate updated codon and colours array
//            if (codon_arr != graph_map.at(kmer_pair.first).part_codon.at(kmer_pair.second).at(modulus))
//            {
//                cache_valid = 1;
//            }

            // calculate condon array
            uint8_t updated_codon_arr = (codon_arr & ~(frame.at(modulus)));

            // update colour array with new information for head/tail k-mer of next unitig
            std::vector<bool> updated_colour_array;
            updated_colour_array = negate_colours_array(colour_arr, graph_map.at(kmer_pair.first).unitig_head_colour);

            // if node is end of contig (e.g. colour change, no connecting nodes) and path not completed
            // path to enable genes at beginning/end of contig to be returned
            if (graph_map.at(kmer_pair.first).end_contig)
            {
                // create temporary path to account for reaching end of contig
                std::vector<std::pair<size_t, bool>> temp_path = head_kmer_list;
                temp_path.push_back(kmer_pair);
                std::pair<std::vector<std::pair<size_t, bool>>, std::vector<bool>> path_pair (std::move(temp_path), colour_arr);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(path_pair));
            }

            // check that colours are viable
            if (updated_colour_array != empty_colour_arr)
            {
//                // get previous cached paths
//                PathVector cached_paths;
//                // check if cached path is present in previous paths
//                if (updated_codon_arr != 0)
//                {
//                    mtx2.lock();
//                    if (previous_paths.find(cached_kmer_id) != previous_paths.end())
//                    {
//                        cached_paths = previous_paths.at(cached_kmer_id);
//                    }
//                    mtx2.unlock();
//                }

                // if repeat is false and unitig has already been traversed ignore unitig
                if (!repeat && is_in)
                {
                    ;
                }
                // check if cached path exists, will only not be empty if next node path already complete and does not complete path
//                else if (!cached_paths.empty())
//                {
//                    // check if node has already been traversed and completed in same reading frame. If so, append completed paths to path list and break
//                    for (const auto& cached_path : cached_paths)
//                    {
//                        // create new_path variable which can be added to from exisiting paths
//                        std::vector<std::pair<size_t, bool>> new_path;
//                        new_path.insert(new_path.end(), head_kmer_list.begin(), head_kmer_list.end());
//
//                        // initialise bool to check if traversing nodes in cached path produces complete path
//                        bool path_complete;
//
//                        // iterate over new cached nodes to calculate length, checking if any new nodes end codon array
//                        std::set<std::pair<size_t, bool>> updated_kmer_set = kmer_set;
//                        size_t updated_length = length;
//                        std::vector<bool> cached_colour_array = colour_arr;
//                        uint8_t cached_codon_arr = codon_arr;
//                        for (const auto& cached_node : cached_path.first)
//                        {
//                            // if end of contig detected, add to path_list without worrying about codons or colours
//                            if (graph_map.at(kmer_pair.first).end_contig)
//                            {
//                                // create temporary path to account for reaching end of contig
//                                std::vector<std::pair<size_t, bool>> temp_path = new_path;
//                                temp_path.push_back(cached_node);
//                                std::pair<std::vector<std::pair<size_t, bool>>, std::vector<bool>> path_pair (std::move(temp_path), cached_colour_array);
//
//                                // update path_list to enable this to be returned
//                                path_list.push_back(std::move(path_pair));
//                            }
//
//                            // calculate updated statistics for new_path
//                            modulus = updated_length % 3;
//                            cached_codon_arr = (codon_arr & ~(graph_map.at(cached_node.first).part_codon.at(cached_node.second).at(modulus)));
//                            cached_colour_array = negate_colours_array(cached_colour_array, graph_map.at(cached_node.first).unitig_head_colour);
//
//                            // check if codon arrays are complete/colour arrays not empty
//                            if (cached_codon_arr != 0 && cached_colour_array != empty_colour_arr)
//                            {
//                                new_path.push_back(cached_node);
//                                updated_length += graph_map.at(cached_node.first).unitig_size.second;
//                                updated_kmer_set.insert(cached_node);
//                            } else if (cached_codon_arr == 0 && cached_colour_array != empty_colour_arr) {
//                                // if path is complete, update accordingly and add to path_list
//                                path_complete = true;
//                                new_path.push_back(cached_node);
//                                std::pair<std::vector<std::pair<size_t, bool>>, std::vector<bool>> path_pair (std::move(new_path), cached_colour_array);
//                                path_list.push_back(path_pair);
//                            } else {
//                                path_complete = true;
//                                break;
//                            }
//                        }
//
//                        // if path now complete, and path still valid (codon array and colour array not 0), need to recurse paths
//                        if (!path_complete && cached_codon_arr != 0 && cached_colour_array != empty_colour_arr)
//                        {
//                            for (const auto& iteration : recur_nodes_binary(graph_map, head_kmer_map, previous_paths, new_path, updated_codon_arr, cached_colour_array, updated_kmer_set, updated_length, length_max, repeat, empty_colour_arr))
//                            {
//                                path_list.push_back(iteration);
//                            }
//                        }
//                    }
//                }
                else {
                    // create new path object for iteration
                    std::vector<std::pair<size_t, bool>> path = head_kmer_list;

                    // add kmer pair to path and kmer_set
                    path.push_back(kmer_pair);
                    std::set<std::pair<size_t, bool>> updated_kmer_set = kmer_set;
                    updated_kmer_set.insert(kmer_pair);

                    // calculate updated length and modulus for codon array iteration
                    size_t updated_length = length;
                    updated_length += graph_map.at(kmer_pair.first).unitig_size.second;

                    // ensure max length has not been exceeded.
                    if (length <= length_max)
                    {
                        for (const auto& iteration : recur_nodes_binary(graph_map, head_kmer_map, path, updated_codon_arr, updated_colour_array, updated_kmer_set, updated_length, length_max, repeat, empty_colour_arr))
                        {
                            path_list.push_back(iteration);
                        }
                    }
                    // if updated length is greater than length_max, but only two nodes present, set codon_arr to 0 and let it be returned, as ORFs in long single nodes need to be returned.
                    else if (path.size() == 2)
                    {
                        updated_codon_arr = 0;
                        for (const auto& iteration : recur_nodes_binary(graph_map, head_kmer_map, path, updated_codon_arr, updated_colour_array, updated_kmer_set, updated_length, length_max, repeat, empty_colour_arr))
                        {
                            path_list.push_back(iteration);
                        }
                    }
                }
            }
        }
        // remove paths that aren't completed
        path_list.erase(path_list.begin());
    }
    return path_list;
}

PathPair traverse_graph(const GraphTuple& graph_tuple,
                         const bool repeat,
                         const vector<bool>& empty_colour_arr,
                         const size_t max_path_length)
{
    // define head_string_arr for parrallelisation
    std::vector<std::string> head_string_arr;

    // recur through nodes to file paths
    PathMap complete_paths;

    cout << "Traversing nodes in forward direction..." << endl;
    // traverse nodes in forward direction
    #pragma omp parallel
    {
        PathVector unitig_complete_paths;
        PathMap complete_paths_private;
        std::vector<std::string> head_string_arr_private;
        #pragma omp for nowait
        for (auto it = std::get<1>(graph_tuple).begin(); it < std::get<1>(graph_tuple).end(); it++)
        {
            const auto head_kmer_id = *it;
            // gather unitig information from graph_map
            const std::pair<size_t, bool> kmer_pair(head_kmer_id, true);
            const uint8_t codon_arr = std::get<0>(graph_tuple).at(head_kmer_id).full_codon.at(true).at(0);
            std::vector<bool> colour_arr = std::get<0>(graph_tuple).at(head_kmer_id).unitig_tail_colour;
            const size_t unitig_len = std::get<0>(graph_tuple).at(head_kmer_id).unitig_size.first;

            // generate vector and set for traversal
            std::vector<std::pair<size_t, bool>> head_kmer_list;
            std::set<std::pair<size_t, bool>> kmer_set;
            head_kmer_list.push_back(kmer_pair);
            kmer_set.insert(kmer_pair);

            // recur paths
            unitig_complete_paths = recur_nodes_binary(std::get<0>(graph_tuple), std::get<3>(graph_tuple), head_kmer_list, codon_arr, colour_arr, kmer_set, unitig_len, max_path_length, repeat, empty_colour_arr);

            // append to complete paths vector
            std::string head_id = std::to_string(head_kmer_id) + "+";
            if (!unitig_complete_paths.empty())
            {
//                //critical section, ensure data race is avoided when appending to/reading complete_paths
//                mtx2.lock();
//                complete_paths.emplace(head_id, std::move(unitig_complete_paths));
//                head_string_arr.push_back(head_id);
//                mtx2.unlock();
                complete_paths_private.emplace(head_id, std::move(unitig_complete_paths));
                head_string_arr_private.push_back(head_id);
            }
        }
        #pragma omp critical
        {
            complete_paths.insert(complete_paths_private.begin(), complete_paths_private.end());
            head_string_arr.insert(head_string_arr.end(), std::make_move_iterator(head_string_arr_private.begin()), std::make_move_iterator(head_string_arr_private.end()));
        }
    }

    cout << "Traversing nodes in reverse direction..." << endl;
    // traverse nodes in forward direction
    #pragma omp parallel
    {
        PathVector unitig_complete_paths;
        PathMap complete_paths_private;
        std::vector<std::string> head_string_arr_private;
        #pragma omp for nowait
        for (auto it = std::get<2>(graph_tuple).begin(); it < std::get<2>(graph_tuple).end(); it++)
        {
            const auto head_kmer_id = *it;
            // gather unitig information from graph_map
            const std::pair<size_t, bool> kmer_pair(head_kmer_id, false);
            const uint8_t codon_arr = std::get<0>(graph_tuple).at(head_kmer_id).full_codon.at(false).at(0);
            std::vector<bool> colour_arr = std::get<0>(graph_tuple).at(head_kmer_id).unitig_head_colour;
            const size_t unitig_len = std::get<0>(graph_tuple).at(head_kmer_id).unitig_size.first;

            // generate vector and set for traversal
            std::vector<std::pair<size_t, bool>> head_kmer_list;
            std::set<std::pair<size_t, bool>> kmer_set;
            head_kmer_list.push_back(kmer_pair);
            kmer_set.insert(kmer_pair);

            // recur paths
            unitig_complete_paths = recur_nodes_binary(std::get<0>(graph_tuple), std::get<3>(graph_tuple), head_kmer_list, codon_arr, colour_arr, kmer_set, unitig_len, max_path_length, repeat, empty_colour_arr);

            // append to complete paths vector
            std::string head_id = std::to_string(head_kmer_id) + "-";
            if (!unitig_complete_paths.empty())
            {
//                //critical section, ensure data race is avoided when appending to/reading complete_paths
//                mtx2.lock();
//                complete_paths.emplace(head_id, std::move(unitig_complete_paths));
//                head_string_arr.push_back(head_id);
//                mtx2.unlock();
                complete_paths_private.emplace(head_id, std::move(unitig_complete_paths));
                head_string_arr_private.push_back(head_id);
            }
        }
        #pragma omp critical
        {
            complete_paths.insert(complete_paths_private.begin(), complete_paths_private.end());
            head_string_arr.insert(head_string_arr.end(), std::make_move_iterator(head_string_arr_private.begin()), std::make_move_iterator(head_string_arr_private.end()));
        }
    }

    //std::sort(head_string_arr.begin(), head_string_arr.end());
    auto path_pair = std::make_pair(complete_paths, head_string_arr);
    return path_pair;
}