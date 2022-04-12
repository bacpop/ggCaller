#include "ORF_clustering.h"

ORFGroupPair group_ORFs(const ColourORFVectorMap& colour_ORF_map,
                         const std::vector<Kmer>& head_kmer_arr,
                         const size_t& overlap)
{
    // intialise Eigen Triplet
    std::vector<ET> tripletList;

    // generate a vector which maps ORF colour/ID to column in eigen matrix
//    ORFMatrixMap ORF_mat_map;

    // initialise length list for sorting of ORFs by length
    std::vector<std::pair<size_t, std::pair<size_t, size_t>>> ORF_length_list;

    // create a mapping between ORF_ID and hash to determine centroids that have same length
//    robin_hood::unordered_map<size_t, size_t> ID_hash_map;

    // initialise vector to hold all centroid IDs in, and dynamic bitset to determine which centroids have been assigned
    std::vector<std::tuple<int, int, size_t>> centroid_vector(head_kmer_arr.size(), {-1, -1, 0});

    // iterate over each ORF sequence with specific colours combination
    size_t ORF_ID = 0;
    for (const auto& colour : colour_ORF_map)
    {
        for (const auto& ORF_entry : colour.second)
        {
            // add to ORF_mat_map, initialising empty vector
//            ORF_mat_map[ORF_ID] = {colour.first, ORF_entry.first};

            // extract info from ORF_entry
            const auto& ORF_info = ORF_entry.second;

            // iterate over nodes traversed by ORF
            const auto& ORF_nodes = std::get<0>(ORF_info);

            ORF_length_list.push_back({std::get<2>(ORF_info), {colour.first, ORF_entry.first}});

            // generate a hash based on kmer string
            std::string ORF_kmer;

            // hold indices that are complete unitigs i.e. start/end not in overlap
            std::vector<size_t> complete_nodes;

            // traverse to generate hash
            for (int i = 0; i < ORF_nodes.size(); i++)
            {
                const auto& node_indices = std::get<1>(ORF_info).at(i);
                // determine if at end and in overlap region of unitigs, if so pass
                if (node_indices.second - node_indices.first < overlap)
                {
                    // if no entries added, then at start so break
                    if (complete_nodes.empty())
                    {
                       continue;
                    }
                    // else at end, so break
                    else
                    {
                        break;
                    }
                }
                // convert node_traversed to size_t, minus 1 as unitigs are one-based, needs to be zero based
                const int& node_traversed = ORF_nodes.at(i);
                const size_t node_ID = abs(node_traversed) - 1;
                ORF_kmer += head_kmer_arr.at(node_ID).toString();
                complete_nodes.push_back(node_ID);
            }

            // get hash to hash map
            size_t ORF_hash = hasher{}(ORF_kmer);

            // get length of ORF entry
            const size_t& ORF_length = std::get<2>(ORF_info);

            // go over again, this time determining centroid status
            for (const auto& node_ID : complete_nodes)
            {
                // determine the current size of the ORF and the centroid if previous assigned
                if (std::get<0>(centroid_vector.at(node_ID)) != -1)
                {
                    auto& centroid_ID = centroid_vector[node_ID];
                    const auto& centroid_info = colour_ORF_map.at(std::get<0>(centroid_ID)).at(std::get<1>(centroid_ID));
//                    const auto& centroid_entry = ORF_mat_map.at(centroid_ID);
                    const auto& centroid_length = std::get<2>(centroid_info);

                    // determine if current ORF is larger than current centroid
                    // if equal, then add lowest ORF index of two as centroid
                    if (ORF_length >= centroid_length)
                    {
                        // if equal, then add lowest ORF index of two as centroid
                        if (ORF_length == centroid_length)
                        {
                            const auto& centroid_hash = std::get<2>(centroid_ID);

                            if (ORF_hash < centroid_hash)
                            {
                                centroid_ID = {colour.first, ORF_entry.first, ORF_hash};
                            } else if (ORF_hash == centroid_hash)
                            {
                                if (colour.first < std::get<0>(centroid_ID))
                                {
                                    centroid_ID = {colour.first, ORF_entry.first, ORF_hash};
                                }
                            }
                        } else
                        {
                            centroid_ID = {colour.first, ORF_entry.first, ORF_hash};
                        }
                    }
                } else
                {
                    centroid_vector[node_ID] = {colour.first, ORF_entry.first, ORF_hash};
                }
            }
            // increment ORF_ID
            ORF_ID++;
        }
    }

    // sort the ORFs in descending order
    std::stable_sort(ORF_length_list.rbegin(), ORF_length_list.rend());

    //return group information
    auto return_pair = std::pair(centroid_vector, ORF_length_list);
    return return_pair;
}

ORFClusterMap produce_clusters(const ColourORFVectorMap& colour_ORF_map,
                               const ColoredCDBG<MyUnitigMap>& ccdbg,
                               const std::vector<Kmer>& head_kmer_arr,
                               const size_t& DBG_overlap,
                               ORFGroupPair& ORF_group_pair,
                               const double& id_cutoff,
                               const double& len_diff_cutoff)
{
    // unpack ORF_group_pair
    auto& centroid_vector = ORF_group_pair.first;
    auto& ORF_length_list = ORF_group_pair.second;

    // initialise temporary cluster map, first item is centroid, second is set of nodes in that cluster
    tbb::concurrent_unordered_map<std::string, tbb::concurrent_unordered_set<size_t>> cluster_map;

    // iterate over all ORFs, checking if all centroids have been identified
    std::vector<size_t> group_indices(ORF_length_list.size());
    std::iota(group_indices.begin(), group_indices.end(), 0);

    // iterate through group_map, ensuring all ORFs have been checked against all potential centroids
    while (true)
    {
        // create set to keep track of assigned ORFs to avoid ORFs being added to other groups
        tbb::concurrent_unordered_set<std::string> encountered_set;

        //create a copy of group_indices for updating from within the for loop
        std::vector<size_t> group_indices_iter;

        // create vector to hold centroid sizes
        std::vector<std::vector<std::pair<size_t, size_t>>> centroid_vector_sizes(centroid_vector.size());

        #pragma omp parallel
        {
            std::vector<std::vector<std::pair<size_t, size_t>>> centroid_vector_sizes_private(centroid_vector.size());
            std::vector<size_t> group_indices_private;

            #pragma omp parallel for
            for (int index = 0; index < group_indices.size(); index++)
            {
                // get entry information for the current ORF
                const size_t& ORF_ID = group_indices.at(index);

                const auto& ORF2_entry = ORF_length_list.at(ORF_ID).second;

                std::string ORF2_entry_str = std::to_string(ORF2_entry.first) + "_" + std::to_string(ORF2_entry.second);

                // check if ORF encountered previously
                if (encountered_set.find(ORF2_entry_str) != encountered_set.end())
                {
                    continue;
                }

                const auto& ORF2_info = colour_ORF_map.at(ORF2_entry.first).at(ORF2_entry.second);
                const auto& ORF2_len = std::get<2>(ORF2_info);

                // create a set for all the centroids for the groups the ORF belongs to, holding hash of ORFs for ordering
                std::set<std::pair<size_t, size_t>> centroid_set;

                // keep track of whether ORF itself is centroid
                bool is_centroid = false;

                // go through the group_IDs
                for (const auto& node_ID : std::get<0>(ORF2_info))
                {
                    size_t group_ID = abs(node_ID) - 1;
                    if (std::get<0>(centroid_vector.at(group_ID)) != -1)
                    {
                        const auto& centroid_ID = centroid_vector.at(group_ID);

                        // check if entry is centroid itself
                        if (std::get<0>(centroid_ID) == ORF2_entry.first && std::get<1>(centroid_ID) == ORF2_entry.second)
                        {
                            is_centroid = true;
                            continue;
                        }

                        // add all centroids to centroid set for shared k-mer groups
                        centroid_set.insert({std::get<2>(centroid_ID), group_ID});
                    }
                }

                // go through and align to each centroid, determine what is the highest scoring assignment
                std::pair<size_t, size_t> assigned_centroid = ORF2_entry;
                double assigned_perc_id = 0;

                for (const auto& centroid_pair : centroid_set)
                {
                    // get ORF_ID
                    const auto& group_ID = centroid_pair.second;
                    const auto& ORF1_entry = centroid_vector.at(group_ID);

                    const auto& ORF1_info = colour_ORF_map.at(std::get<0>(ORF1_entry)).at(std::get<1>(ORF1_entry));
                    const auto& ORF1_len = std::get<2>(ORF1_info);

                    // check that perc_len_diff is greater than cut-off, otherwise pass
                    double perc_len_diff = (double)ORF2_len / (double)ORF1_len;

                    if (perc_len_diff < len_diff_cutoff)
                    {
                        continue;
                    }

                    // calculate the perc_id between the current centroid and the ORF
                    double current_perc_id = align_seqs(ORF1_info, ORF2_info, ccdbg, head_kmer_arr, DBG_overlap);

                    // assign to first centroid with current_perc_id >= id_cutoff
                    if (current_perc_id >= id_cutoff)
                    {
                        assigned_perc_id = current_perc_id;
                        assigned_centroid = {std::get<0>(ORF1_entry), std::get<1>(ORF1_entry)};
                        break;
                    }
                }

                // check if centroid perc_id is greater than cut-off, if not then assign the cluster to it's own centroid
                if (assigned_perc_id > 0 || is_centroid)
                {
                    std::string assigned_centroid_str = std::to_string(assigned_centroid.first) + "_" + std::to_string(assigned_centroid.second);
                    cluster_map[assigned_centroid_str].insert(ORF_ID);
                    encountered_set.insert(assigned_centroid_str);
                }
                else
                {
                    // if singleton and not centroid, add the current ORF as a centroid for all groups it is part of
                    // enabling searching in next iteration
                    for (const auto& node_ID : std::get<0>(ORF2_info))
                    {
                        // add ORF as new centroid of group
                        centroid_vector_sizes_private[abs(node_ID) - 1].push_back({ORF2_len, ORF_ID});
                    }
                    // add entry to group_indices_private to re-cluster in next pass
                    group_indices_private.push_back(ORF_ID);
                }
                // add to encountered_set
                encountered_set.insert(ORF2_entry_str);
            }

            #pragma omp critical
            {
                // clear centroid vector and reset with all new centroids
                for (int i = 0; i < centroid_vector_sizes_private.size(); i++)
                {
                    centroid_vector_sizes[i].insert(centroid_vector_sizes[i].end(), std::make_move_iterator(centroid_vector_sizes_private[i].begin()), std::make_move_iterator(centroid_vector_sizes_private[i].end()));
                }

                // update group_indices
                group_indices_iter.insert(group_indices_iter.end(), std::make_move_iterator(group_indices_private.begin()), std::make_move_iterator(group_indices_private.end()));
            }
        }

        // update group_indices with group_indices_iter
        group_indices = std::move(group_indices_iter);
        if (!group_indices.size())
        {
            break;
        } else
        {
            // reset centroid vector
            for (int i = 0; i < centroid_vector.size(); i++)
            {
                centroid_vector[i] = {-1, -1, 0};
            }

            // add only longest sequences to centroid_vector
            for (int i = 0; i < centroid_vector_sizes.size(); i++)
            {
                // sort each entry by size, take largest
                auto& group_entry = centroid_vector_sizes[i];
                // if group_entry not empty, sort and then take largest
                if (!group_entry.empty())
                {
                    std::stable_sort(group_entry.rbegin(), group_entry.rend());
                    auto& largest = group_entry[0];
                    auto& largest_entry = ORF_length_list.at(largest.second).second;

                    // generate a hash for largest entry
                    size_t largest_hash;
                    {
                        const auto& ORF_info = colour_ORF_map.at(largest_entry.first).at(largest_entry.second);

                        // generate a hash based on kmer string
                        std::string ORF_kmer;

                        // hold indices that are complete unitigs i.e. start/end not in overlap
                        std::vector<size_t> complete_nodes;

                        // traverse to generate hash
                        for (int k = 0; k < std::get<0>(ORF_info).size(); k++)
                        {
                            const auto& node_indices = std::get<1>(ORF_info).at(k);
                            // determine if at end and in overlap region of unitigs, if so pass
                            if (node_indices.second - node_indices.first < DBG_overlap)
                            {
                                // if no entries added, then at start so break
                                if (complete_nodes.empty())
                                {
                                    continue;
                                }
                                    // else at end, so break
                                else
                                {
                                    break;
                                }
                            }
                            // convert node_traversed to size_t, minus 1 as unitigs are one-based, needs to be zero based
                            const int& node_traversed = std::get<0>(ORF_info).at(k);
                            const size_t node_ID = abs(node_traversed) - 1;
                            ORF_kmer += head_kmer_arr.at(node_ID).toString();
                            complete_nodes.push_back(node_ID);
                        }

                        // get hash to hash map
                        largest_hash = hasher{}(ORF_kmer);
                    }

                    for (int j = 1; j < group_entry.size(); j++)
                    {
                        if (group_entry[j].first < largest.first)
                        {
                            break;
                        } else if (group_entry[j].first == largest.first)
                        {
                            // if size is equal, then take lowest hash as centroid
                            const auto& new_entry = ORF_length_list.at(group_entry[j].second).second;

                            // generate new hash for ORF
                            size_t new_entry_hash;
                            {
                                const auto& ORF_info = colour_ORF_map.at(new_entry.first).at(new_entry.second);

                                // generate a hash based on kmer string
                                std::string ORF_kmer;

                                // hold indices that are complete unitigs i.e. start/end not in overlap
                                std::vector<size_t> complete_nodes;

                                // traverse to generate hash
                                for (int k = 0; k < std::get<0>(ORF_info).size(); k++)
                                {
                                    const auto& node_indices = std::get<1>(ORF_info).at(k);
                                    // determine if at end and in overlap region of unitigs, if so pass
                                    if (node_indices.second - node_indices.first < DBG_overlap)
                                    {
                                        // if no entries added, then at start so break
                                        if (complete_nodes.empty())
                                        {
                                            continue;
                                        }
                                            // else at end, so break
                                        else
                                        {
                                            break;
                                        }
                                    }
                                    // convert node_traversed to size_t, minus 1 as unitigs are one-based, needs to be zero based
                                    const int& node_traversed = std::get<0>(ORF_info).at(k);
                                    const size_t node_ID = abs(node_traversed) - 1;
                                    ORF_kmer += head_kmer_arr.at(node_ID).toString();
                                    complete_nodes.push_back(node_ID);
                                }

                                // get hash to hash map
                                new_entry_hash = hasher{}(ORF_kmer);
                            }

                            if (new_entry_hash < largest_hash)
                            {
                                largest = group_entry[j];
                                largest_hash = new_entry_hash;
                            } else if (new_entry_hash == largest_hash)
                            {
                                if (new_entry.first < largest_entry.first)
                                {
                                    largest = group_entry[j];
                                    largest_hash = new_entry_hash;
                                }
                            }
                        }
                    }

                    // assign ID to centroid vector
                    centroid_vector[i] = {largest_entry.first, largest_entry.second, largest_hash};

                    // add the entry to cluster_map
                    std::string largest_str = std::to_string(largest_entry.first) + "_" + std::to_string(largest_entry.second);
                    cluster_map[largest_str].insert(largest.second);
                }
            }
        }
    }

    // generate set to determine which ORFs have already been assigned clusters
    std::set<std::pair<size_t, size_t>> cluster_assigned;

    // initialise map as intermediate to hold cluster IDs
    ORFClusterMap final_clusters;

    // iterate over ORF_length_list, pulling out centroids and their assigned clustered ORFs
    size_t cluster_id = 0;
    for (size_t i = 0; i < ORF_length_list.size(); i++)
    {
        const auto& ORF_ID = ORF_length_list.at(i).second;

        if (cluster_assigned.find(ORF_ID) != cluster_assigned.end())
        {
            continue;
        }

        std::string ORF_ID_str = std::to_string(ORF_ID.first) + "_" + std::to_string(ORF_ID.second);

        // check if entry is a centroid, if so add to final_clusters along with all of its attached ORFs that are unassigned
        if (cluster_map.find(ORF_ID_str) != cluster_map.end())
        {
            // add current ORF to centroid entry
            final_clusters[cluster_id].push_back(ORF_ID);

            // add the centroid to the cluster_assigned set
            cluster_assigned.insert(ORF_ID);

            // add rest of homologs to centroid entry
            for (const auto& homolog_alias : cluster_map.at(ORF_ID_str))
            {
                const auto& homolog_ID = ORF_length_list.at(homolog_alias).second;

                std::string homolog_ID_str = std::to_string(homolog_ID.first) + "_" + std::to_string(homolog_ID.second);

                // if homolog is a centroid, do not add
                if (cluster_map.find(homolog_ID_str) != cluster_map.end())
                {
                    continue;
                }

                // if the homolog is not already assigned to a cluster, assign and add to cluster_assigned
                if (cluster_assigned.find(homolog_ID) == cluster_assigned.end())
                {
                    final_clusters[cluster_id].push_back(homolog_ID);
                    cluster_assigned.insert(homolog_ID);
                }
                // if homolog already assigned to a cluster, skip
                else
                {
                    continue;
                }
            }
        }
        cluster_id++;
    }

    return final_clusters;
}

double align_seqs(const ORFNodeVector& ORF1_info,
                  const ORFNodeVector& ORF2_info,
                  const ColoredCDBG<MyUnitigMap>& ccdbg,
                  const std::vector<Kmer>& head_kmer_arr,
                  const size_t& DBG_overlap)
{
    std::string ORF1_aa;
    std::string ORF2_aa;

    // scope for DNA sequences
    {
        const auto ORF1_sequence = generate_sequence_nm(std::get<0>(ORF1_info), std::get<1>(ORF1_info), DBG_overlap, ccdbg, head_kmer_arr);
        const auto ORF2_sequence = generate_sequence_nm(std::get<0>(ORF2_info), std::get<1>(ORF2_info), DBG_overlap, ccdbg, head_kmer_arr);

        // translate DNA sequence
        ORF1_aa = (translate(ORF1_sequence)).aa();
        ORF2_aa = (translate(ORF2_sequence)).aa();
    }

    // convert sequence to const char * and calculate edit distance
    const char * seq1 = ORF1_aa.c_str();
    const char * seq2 = ORF2_aa.c_str();
    EdlibAlignResult result = edlibAlign(seq1, ORF1_aa.size(), seq2, ORF1_aa.size(), edlibDefaultAlignConfig());
    size_t edit_distance = result.editDistance;
    edlibFreeAlignResult(result);

    // convert edit distance into percent identity
    double perc_id = 1 - ((double)edit_distance / (double)ORF2_aa.size());
    return perc_id;
}