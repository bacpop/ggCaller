// ggCaller header
#include "ggCaller_classes.h"
// define mutex for safe addition to robinhood_maps
std::mutex mtx;

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
                               const PathMap& previous_paths,
                               const std::vector<std::pair<std::string, bool>>& head_kmer_list,
                               const uint8_t& codon_arr,
                               const std::vector<bool>& colour_arr,
                               const std::set<std::pair<std::string, bool>>& kmer_set,
                               const size_t& length,
                               const bool& forward,
                               const size_t& length_max,
                               const bool& repeat,
                               const vector<bool>& empty_colour_arr)
{
    // generate path list, add head_kmer_list to it
    PathVector path_list;
    path_list.push_back(std::make_pair(head_kmer_list, colour_arr));

    // if path not complete, generate unitig object using head kmer
    if (codon_arr != 0)
    {
        // get last head-kmer pair in list, generate unitig map
        const std::pair<size_t, bool>& iter_kmer_pair = head_kmer_list.back();
        const iter_kmer_str = graph_map.at(iter_kmer_pair.first).head_kmer;
        const Kmer iter_km = Kmer(iter_km_str.c_str());
        auto unitig = ccdbg.find(iter_km, true);

        unitig.strand = iter_km_str.second;

        // generate neighbours, assume successors are correct
        const auto& next_unitigs = unitig.getSuccessors();

        // get predecessors if next node strand is reversed
        if ((forward && !unitig.strand) || (!forward && unitig.strand))
        {
            const auto& next_unitigs = unitig.getPredecessors();
        }

        // iterate over neighbours, recurring through incomplete paths
        for (const auto& um : next_unitigs)
        {

            // create new path object for iteration
            std::vector<std::pair<std::string, bool>> path = head_kmer_list;
            std::string kmer_id = um.getUnitigHead().toString();

            // make path entry
            auto kmer_pair = std::make_pair(kmer_id, um.strand);

            // check if unitig has already been traversed
            const bool is_in = kmer_set.find(kmer_pair) != kmer_set.end();

            // calculate modulus for path and updated cached path reading frame
            int modulus = length % 3;
            bool cache_valid = 0;
            std::string cached_kmer_id = kmer_id + (um.strand ? "+" : "-");

            // calculate updated codon and colours array
            if (codon_arr == graph_map.at(kmer_id).part_codon.at(um.strand).at(modulus))
            {
                cache_valid = 1;
            }
            //const auto updated_codon_arr = compare_codon_array(codon_arr, graph_map.at(kmer_id).part_codon.at(um.strand).at(modulus));

            const uint8_t updated_codon_arr = (codon_arr & ~(graph_map.at(kmer_id).part_codon.at(um.strand).at(modulus)));
            const auto updated_colour_array = negate_colours_array(colour_arr, graph_map.at(kmer_id).unitig_colour);

            // check that colours are viable
            if (updated_colour_array != empty_colour_arr)
            {
                // get previous cached paths
                PathVector cached_paths;
                // check if cached path is present in previous paths
                if (cache_valid)
                {
                    mtx.lock();
                    if (previous_paths.find(cached_kmer_id) != previous_paths.end())
                    {
                        cached_paths = previous_paths.at(cached_kmer_id);
                    }
                    mtx.unlock();
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
                            std::pair<std::vector<std::pair<std::string, bool>>, std::vector<bool>> new_path;
                            new_path.first.insert(new_path.first.end(), path.begin(), path.end());
                            std::pair<std::vector<std::pair<std::string, bool>>, std::vector<bool>> original_path;
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
                    std::set<std::pair<std::string, bool>> updated_kmer_set = kmer_set;
                    updated_kmer_set.insert(kmer_pair);

                    // calculate updated length and modulus for codon array iteration
                    int updated_length = length;
                    updated_length += graph_map.at(kmer_id).unitig_size.second;

                    // ensure max length has not been exceeded
                    if (length <= length_max)
                    {
                        for (const auto& iteration : recur_nodes_binary(ccdbg, graph_map, previous_paths, path, updated_codon_arr, updated_colour_array, updated_kmer_set, updated_length, forward, length_max, repeat, empty_colour_arr))
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
    bool is_forward = true;
    #pragma omp parallel
    {
        PathVector unitig_complete_paths;
        #pragma omp for nowait
        for (auto it = std::get<1>(graph_tuple).begin(); it < std::get<1>(graph_tuple).end(); it++)
        {
            const auto unitig_id = *it;
            // gather unitig information from graph_map
            const std::pair<size_t, bool> kmer_pair(unitig_id, is_forward);
            const uint8_t codon_arr = std::get<0>(graph_tuple).at(unitig_id).full_codon.at(is_forward).at(0);
            const vector<bool> colour_arr = std::get<0>(graph_tuple).at(unitig_id).unitig_colour;
            const size_t unitig_len = std::get<0>(graph_tuple).at(unitig_id).unitig_size.first;

            // generate vector and set for traversal
            std::vector<std::pair<size_t, bool>> head_kmer_list;
            std::set<std::pair<size_t, bool>> kmer_set;
            head_kmer_list.push_back(kmer_pair);
            kmer_set.insert(kmer_pair);

            // recur paths
            unitig_complete_paths = recur_nodes_binary(ccdbg, std::get<0>(graph_tuple), complete_paths, head_kmer_list, codon_arr, colour_arr, kmer_set, unitig_len, is_forward, max_path_length, repeat, empty_colour_arr);

            // append to complete paths vector
            std::string head_kmer_id = std::to_string(unitig_id) + "+";
            if (!unitig_complete_paths.empty())
            {
                //critical section, ensure data race is avoided when appending to/reading complete_paths
                mtx.lock();
                complete_paths.emplace(head_kmer_id, std::move(unitig_complete_paths));
                head_string_arr.push_back(head_kmer_id);
                mtx.unlock();
            }
        }
    }

    cout << "Traversing nodes in reverse direction..." << endl;
    // traverse nodes in reverse direction
    is_forward = false;
    #pragma omp parallel
    {
        PathVector unitig_complete_paths;
        #pragma omp for nowait
        for (auto it = std::get<2>(graph_tuple).begin(); it < std::get<2>(graph_tuple).end(); it++)
        {
            const auto unitig_id = *it;
            // gather unitig information from graph_map
            const std::pair<size_t, bool> kmer_pair(unitig_id, is_forward);
            const uint8_t codon_arr = std::get<0>(graph_tuple).at(unitig_id).full_codon.at(is_forward).at(0);
            const vector<bool> colour_arr = std::get<0>(graph_tuple).at(unitig_id).unitig_colour;
            const size_t unitig_len = std::get<0>(graph_tuple).at(unitig_id).unitig_size.first;

            // generate vector and set for traversal
            std::vector<std::pair<size_t, bool>> head_kmer_list;
            std::set<std::pair<size_t, bool>> kmer_set;
            head_kmer_list.push_back(kmer_pair);
            kmer_set.insert(kmer_pair);

            // recur paths
            unitig_complete_paths = recur_nodes_binary(ccdbg, std::get<0>(graph_tuple), complete_paths, head_kmer_list, codon_arr, colour_arr, kmer_set, unitig_len, is_forward, max_path_length, repeat, empty_colour_arr);

            // append to complete paths vector
            std::string head_kmer_id = head_kmer + "-";
            if (!unitig_complete_paths.empty())
            {
                // critical section, ensure data race is avoided when appending to/reading complete_paths
                mtx.lock();
                complete_paths.emplace(head_kmer_id, std::move(unitig_complete_paths));
                head_string_arr.push_back(head_kmer_id);
                mtx.unlock();
            }
        }
    }

    //std::sort(head_string_arr.begin(), head_string_arr.end());
    auto path_pair = std::make_pair(complete_paths, head_string_arr);
    return path_pair;
}