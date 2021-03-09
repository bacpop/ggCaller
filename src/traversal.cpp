// ggCaller header
#include "ggCaller_classes.h"
// define mutex for safe addition to robinhood_maps
std::mutex mtx2;

// compares two codon arrays, set's value to 0 if both are 1. Creates new array, does not edit either array
std::vector<bool> compare_codon_array(const std::vector<bool>& array1, const std::vector<bool>& array2)
{
    std::vector<bool> output_array = array1;
    for (size_t i = 0; i < 3; i++)
    {
        if (array1[i] == 1 && array2[i] == 1)
        {
            output_array[i] = 0;
        }
    }
    return output_array;
}

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

PathVector recur_nodes_binary (const ColoredCDBG<>& ccdbg,
                               const unitigMap& graph_map,
                               const robin_hood::unordered_map<std::string, size_t>& head_kmer_map,
                               const PathMap& previous_paths,
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

    // if path not complete, generate unitig object using head kmer
    if (codon_arr != 0)
    {
        // get last head-kmer pair in list, generate unitig map object
        const std::pair<size_t, bool>& iter_km_pair = head_kmer_list.back();

        //convert head_kmer_id to kmer sequence
        const std::string iter_km_str = graph_map.at(iter_km_pair.first).head_kmer;

        const Kmer iter_km = Kmer(iter_km_str.c_str());
        auto unitig = ccdbg.find(iter_km, true);


        // parse unitig strand
        unitig.strand = iter_km_pair.second;

        // generate neighbours, assume successors are correct
        const auto& next_unitigs = unitig.getSuccessors();

        // iterate over neighbours, recurring through incomplete paths
        for (const auto& um : next_unitigs)
        {
            // create new path object for iteration
            std::vector<std::pair<size_t, bool>> path = head_kmer_list;
            const std::string kmer_str = um.getUnitigHead().toString();
            size_t kmer_id = head_kmer_map.at(kmer_str);

            // make path entry
            bool unitig_strand = um.strand;
            const auto kmer_pair = make_pair(kmer_id, unitig_strand);

            // check if unitig has already been traversed
            const bool is_in = kmer_set.find(kmer_pair) != kmer_set.end();

            // calculate modulus for path and updated cached path reading frame
            int modulus = length % 3;
            bool cache_valid = 0;
            std::string cached_kmer_id = std::to_string(kmer_pair.first) + (um.strand ? "+" : "-");

            // calculate updated codon and colours array
            if (codon_arr == graph_map.at(kmer_pair.first).part_codon.at(um.strand).at(modulus))
            {
                cache_valid = 1;
            }

            // calculate condon array
            uint8_t updated_codon_arr = (codon_arr & ~(graph_map.at(kmer_pair.first).part_codon.at(um.strand).at(modulus)));

            // update colour array with new information for head/tail k-mer of next unitig
            std::vector<bool> updated_colour_array;
            updated_colour_array = negate_colours_array(colour_arr, graph_map.at(kmer_pair.first).unitig_head_colour);

            // if node is end of contig (e.g. colour change, no connecting nodes)
            // path to enable genes at beginning/end of contig to be returned
            if (graph_map.at(kmer_pair.first).end_contig)
            {
                // create temporary path to account for reaching end of contig
                std::vector<std::pair<size_t, bool>> temp_path = path;
                temp_path.push_back(kmer_pair);
                std::pair<std::vector<std::pair<size_t, bool>>, std::vector<bool>> path_pair (std::move(temp_path), colour_arr);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(path_pair));
            }

            // check that colours are viable
            if (updated_colour_array != empty_colour_arr)
            {
                // get previous cached paths
                PathVector cached_paths;
                // check if cached path is present in previous paths
                if (cache_valid)
                {
                    mtx2.lock();
                    if (previous_paths.find(cached_kmer_id) != previous_paths.end())
                    {
                        cached_paths = previous_paths.at(cached_kmer_id);
                    }
                    mtx2.unlock();
                }

                // if repeat is false and unitig has already been traversed ignore unitig
                if (!repeat && is_in)
                {
                    ;
                }
                // check if cached path exists, if following unitig does not complete path, and that cached path completes path
                else if (!cached_paths.empty() && (updated_codon_arr != 0) && (cache_valid))
                {
                    // check if node has already been traversed and completed in same reading frame. If so, append completed paths to path list and break
                    for (const auto& cached_path : cached_paths)
                    {
                        // check that path colours match
                        const auto cached_colour_array = negate_colours_array(colour_arr, cached_path.second);
                        if (cached_colour_array != empty_colour_arr)
                        {
                            // add cached paths to previously traversed nodes
                            std::pair<std::vector<std::pair<size_t, bool>>, std::vector<bool>> new_path;
                            new_path.first.insert(new_path.first.end(), path.begin(), path.end());
                            std::pair<std::vector<std::pair<size_t, bool>>, std::vector<bool>> original_path;
                            original_path.first.insert(original_path.first.end(), path.begin(), path.end());
                            new_path.first.insert(new_path.first.end(), cached_path.first.begin(), cached_path.first.end());
                            new_path.second = cached_colour_array;
                            path_list.push_back(std::move(new_path));
                        }
                    }
                }
                else {
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
                        for (const auto& iteration : recur_nodes_binary(ccdbg, graph_map, head_kmer_map, previous_paths, path, updated_codon_arr, updated_colour_array, updated_kmer_set, updated_length, length_max, repeat, empty_colour_arr))
                        {
                            path_list.push_back(iteration);
                        }
                    }
                    // if updated length is greater than length_max, but only two nodes present, set codon_arr to 0 and let it be returned, as ORFs in long single nodes need to be returned.
                    else if (path.size() == 2)
                    {
                        updated_codon_arr = 0;
                        for (const auto& iteration : recur_nodes_binary(ccdbg, graph_map, head_kmer_map, previous_paths, path, updated_codon_arr, updated_colour_array, updated_kmer_set, updated_length, length_max, repeat, empty_colour_arr))
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

PathPair traverse_graph(const ColoredCDBG<>& ccdbg,
                         const GraphTuple& graph_tuple,
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
            unitig_complete_paths = recur_nodes_binary(ccdbg, std::get<0>(graph_tuple), std::get<3>(graph_tuple), complete_paths, head_kmer_list, codon_arr, colour_arr, kmer_set, unitig_len, max_path_length, repeat, empty_colour_arr);

            // append to complete paths vector
            std::string head_id = std::to_string(head_kmer_id) + "+";
            if (!unitig_complete_paths.empty())
            {
                //critical section, ensure data race is avoided when appending to/reading complete_paths
                mtx2.lock();
                complete_paths.emplace(head_id, std::move(unitig_complete_paths));
                head_string_arr.push_back(head_id);
                mtx2.unlock();
            }
        }
    }

    cout << "Traversing nodes in reverse direction..." << endl;
    // traverse nodes in forward direction
    #pragma omp parallel
    {
        PathVector unitig_complete_paths;
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
            unitig_complete_paths = recur_nodes_binary(ccdbg, std::get<0>(graph_tuple), std::get<3>(graph_tuple), complete_paths, head_kmer_list, codon_arr, colour_arr, kmer_set, unitig_len, max_path_length, repeat, empty_colour_arr);

            // append to complete paths vector
            std::string head_id = std::to_string(head_kmer_id) + "-";
            if (!unitig_complete_paths.empty())
            {
                //critical section, ensure data race is avoided when appending to/reading complete_paths
                mtx2.lock();
                complete_paths.emplace(head_id, std::move(unitig_complete_paths));
                head_string_arr.push_back(head_id);
                mtx2.unlock();
            }
        }
    }

    //std::sort(head_string_arr.begin(), head_string_arr.end());
    auto path_pair = std::make_pair(complete_paths, head_string_arr);
    return path_pair;
}