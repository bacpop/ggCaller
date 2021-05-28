// ggCaller header
#include "traversal.h"
// define mutex for safe addition to robinhood_maps
//std::mutex mtx2;

PathVector recur_nodes_binary (const UnitigVector& graph_vector,
                               const std::vector<int>& head_kmer_list,
                               const uint8_t& codon_arr,
                               const size_t& colour_ID,
                               const std::unordered_set<int>& kmer_set,
                               const size_t& length,
                               const size_t& length_max,
                               const bool& repeat)
{
    // generate path list, add head_kmer_list to it
    PathVector path_list;
    path_list.push_back(head_kmer_list);

    // testing
//    std::set<std::pair<std::string, bool>> node_set;
//    std::vector<std::pair<std::string, bool>> node_vector;
//
//    //testing
//    for (const auto& node : head_kmer_list)
//    {
//        std::string node_str = graph_vector.at(abs(node.first) - 1).head_kmer;
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
        // get last unitig id in list, this will be an int and needs to be converted to size_t
        const int& iter_id = head_kmer_list.back();

        // get unitig_dict entry in graph_vector
        const auto& unitig_dict = graph_vector.at(abs(iter_id) - 1);

        // determine strand of unitig
        bool strand = (iter_id >= 0) ? true : false;

        // iterate over neighbours, recurring through incomplete paths
        for (const auto& neighbour : unitig_dict.get_neighbours(strand))
        {
            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const auto& neighbour_id = neighbour.first;
            const auto& frame = neighbour.second;

            // get reference to unitig_dict object for neighbour
            const auto& neighbour_dict = graph_vector.at(abs(neighbour_id) - 1);

            // determine if neighbour is in same colour as iteration, if not pass
            if (!unitig_dict.full_colour().at(colour_ID))
            {
                continue;
            }

            // check if unitig has already been traversed
            const bool is_in = kmer_set.find(neighbour_id) != kmer_set.end();

            // calculate modulus for path and updated cached path reading frame
            int modulus = length % 3;

            // calculate condon array
            uint8_t updated_codon_arr = (codon_arr & ~(frame.at(modulus)));

            // if node is end of contig (e.g. colour change, no connecting nodes) and path not completed
            // path to enable genes at beginning/end of contig to be returned
            if (neighbour_dict.end_contig())
            {
                // create temporary path to account for reaching end of contig
                std::vector<int> temp_path = head_kmer_list;
                temp_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(temp_path));
            }

            // if repeat is false and unitig has already been traversed ignore unitig
            if (!repeat && is_in)
            {
                continue;
            }
            else {
                // create new path object for iteration
                std::vector<int> path = head_kmer_list;

                // add kmer pair to path and kmer_set
                path.push_back(neighbour_id);
                std::unordered_set<int> updated_kmer_set = kmer_set;
                updated_kmer_set.insert(neighbour_id);

                // calculate updated length and modulus for codon array iteration
                size_t updated_length = length;
                updated_length += neighbour_dict.size().second;

                // ensure max length has not been exceeded.
                if (length <= length_max)
                {
                    for (const auto& iteration : recur_nodes_binary(graph_vector, path, updated_codon_arr, colour_ID, updated_kmer_set, updated_length, length_max, repeat))
                    {
                        path_list.push_back(iteration);
                    }
                }
                // if updated length is greater than length_max, but only two nodes present, set codon_arr to 0 and let it be returned, as ORFs in long single nodes need to be returned.
                else if (path.size() == 2)
                {
                    updated_codon_arr = 0;
                    for (const auto& iteration : recur_nodes_binary(graph_vector, path, updated_codon_arr, colour_ID, updated_kmer_set, updated_length, length_max, repeat))
                    {
                        path_list.push_back(iteration);
                    }
                }
            }
        }
        // remove paths that aren't completed
        path_list.erase(path_list.begin());
    }
    return path_list;
}

AllPaths traverse_graph(const UnitigVector& graph_vector,
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

        // generate vector and set for traversal
        std::vector<int> head_kmer_list;
        std::unordered_set<int> kmer_set;
        head_kmer_list.push_back(head_id);
        kmer_set.insert(head_id);

        // recur paths
        PathVector unitig_complete_paths = recur_nodes_binary(graph_vector, head_kmer_list, codon_arr, colour_ID, kmer_set, unitig_len, max_path_length, repeat);

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

        // generate vector and set for traversal
        std::vector<int> head_kmer_list;
        std::unordered_set<int> kmer_set;
        head_kmer_list.push_back(head_id);
        kmer_set.insert(head_id);

        // recur paths
        PathVector unitig_complete_paths = recur_nodes_binary(graph_vector, head_kmer_list, codon_arr, colour_ID, kmer_set, unitig_len, max_path_length, repeat);

        if (!unitig_complete_paths.empty())
        {
            all_paths.push_back(std::move(unitig_complete_paths));
        }
    }

    return all_paths;
}