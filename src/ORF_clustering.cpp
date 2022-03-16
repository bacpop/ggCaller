#include "ORF_clustering.h"

ORFGroupTuple group_ORFs(const ColourORFVectorMap& colour_ORF_map,
                         const std::vector<Kmer>& head_kmer_arr)
{
    // intialise Eigen Triplet
    std::vector<ET> tripletList;

    // generate a vector which maps ORF colour/ID to column in eigen matrix
    ORFMatrixVector ORF_mat_vector;

    // create a mapping between ORF_ID and hash to determine centroids that have same length
    robin_hood::unordered_map<size_t, size_t> ID_hash_map;

    // iterate over each ORF sequence with specific colours combination
    size_t ORF_ID = 0;
    for (const auto& colour : colour_ORF_map)
    {
        for (int i = 0; i < colour.second.size(); i++)
        {
            // add to ORF_mat_map, initialising empty vector
            ORF_mat_vector.push_back({colour.first, i});

            // iterate over nodes traversed by ORF
            const auto& ORF_nodes = std::get<0>(colour.second.at(i));

            // generate a hash based on start and end node Kmer
            {
                std::string ORF_kmer = head_kmer_arr.at(abs(ORF_nodes.at(0)) - 1).toString();
                ORF_kmer += head_kmer_arr.at(abs(ORF_nodes.back()) - 1).toString();
                ID_hash_map[ORF_ID] = hasher{}(ORF_kmer);
            }

            for (const auto& node_traversed : ORF_nodes)
            {
                // add to triplet list, with temp_ORF_ID (row), node id (column) and set value as 1 (true)
                // convert node_traversed to size_t, minus 1 as unitigs are one-based, needs to be zero based
                const size_t abs_node_id = abs(node_traversed) - 1;
                tripletList.push_back(ET(ORF_ID, abs_node_id, 1));
            }
            // increment ORF_ID
            ORF_ID++;
        }
    }

    // initialise map which holds which nodes map to which ORFs
    std::vector<std::unordered_set<size_t>> ORF_group_vector(ORF_mat_vector.size());

    // initialise vector to hold all centroid IDs in, and dynamic bitset to determine which centroids have been assigned
    std::vector<size_t> centroid_vector(head_kmer_arr.size());
    boost::dynamic_bitset<> centroid_assigned(head_kmer_arr.size());

    // initialise sparse matrix
    Eigen::SparseMatrix<double> mat(ORF_mat_vector.size(), head_kmer_arr.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    // iterate over non-zero entries in matrix and calculate overlaps
    for (int outit = 0; outit < mat.outerSize(); ++outit)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator init(mat, outit); init; ++init)
        {
            // get column id (node) and row ID (ORF)
            const auto& node_ID = init.col();
            const auto& ORF_ID = init.row();

            // add the ORF to the group
            ORF_group_vector[ORF_ID].insert(node_ID);

            // get reference to entry in ORF_mat_vector
            const auto& ORF_entry = ORF_mat_vector.at(ORF_ID);

            // get length of ORF entry
            const size_t& ORF_length = std::get<2>(colour_ORF_map.at(ORF_entry.first).at(ORF_entry.second));

            // determine the current size of the ORF and the centroid if previous assigned
            if ((bool)centroid_assigned[node_ID])
            {
                auto& centroid_ID = centroid_vector[node_ID];
                const auto& centroid_entry = ORF_mat_vector.at(centroid_ID);
                const auto& centroid_length = std::get<2>(colour_ORF_map.at(centroid_entry.first).at(centroid_entry.second));

                // determine if current ORF is larger than current centroid
                // if equal, then add lowest ORF index of two as centroid
                if (ORF_length >= centroid_length)
                {
                    // if equal, then add lowest ORF index of two as centroid
                    if (ORF_length == centroid_length)
                    {
                        if (ID_hash_map.at(ORF_ID) < ID_hash_map.at(centroid_ID))
                        {
                            centroid_ID = ORF_ID;
                        }
                    } else
                    {
                        centroid_ID = ORF_ID;
                    }
                }
            } else
            {
                centroid_vector[node_ID] = ORF_ID;
                centroid_assigned[node_ID] = 1;
            }
        }
    }

    //return group information
    auto return_tuple = std::make_tuple(ORF_mat_vector, ORF_group_vector, centroid_vector, ID_hash_map);
    return return_tuple;
}

ORFClusterMap produce_clusters(const ColourORFVectorMap& colour_ORF_map,
                               const ColoredCDBG<MyUnitigMap>& ccdbg,
                               const std::vector<Kmer>& head_kmer_arr,
                               const size_t& DBG_overlap,
                               const ORFMatrixVector& ORF_mat_vector,
                               const std::vector<std::unordered_set<size_t>>& ORF_group_vector,
                               std::vector<size_t>& centroid_vector,
                               const robin_hood::unordered_map<size_t, size_t>& ID_hash_map,
                               const double& id_cutoff,
                               const double& len_diff_cutoff)
{
    // initialise temporary cluster map, first item is centroid, second is set of nodes in that cluster
    tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_set<size_t>> cluster_map;

    // initialise final cluster map for return
    ORFClusterMap final_clusters;

    // initialise length list for sorting of ORFs by length
    std::vector<std::pair<size_t, size_t>> ORF_length_list;

    // go through all comparisons, identifying which ORFs do not find good alignments to current centroids, then repeat for only those ORFs
    // not assigned to centroid (i.e. repeat twice, should cluster all ORFs)

    // iterate over all ORFs, checking if all centroids have been identified
    std::vector<size_t> group_indices(ORF_group_vector.size());
    std::iota(group_indices.begin(), group_indices.end(), 0);

    // keep track of assigned bitsets
    boost::dynamic_bitset<> centroid_assigned(centroid_vector.size());
    centroid_assigned.set();

    // iterate through group_map, ensuring all ORFs have been checked against all potential centroids
    while (true)
    {
        // create set to keep track of assigned ORFs to avoid ORFs being added to other groups
        tbb::concurrent_unordered_set<size_t> encountered_set;

        //create a copy of group_indices for updating from within the for loop
        std::vector<size_t> group_indices_iter;

        // create vector to hold centroid sizes
        std::vector<std::vector<std::pair<size_t, size_t>>> centroid_vector_sizes(centroid_vector.size());

        #pragma omp parallel
        {
            std::vector<std::pair<size_t, size_t>> ORF_length_list_private;
            std::vector<std::vector<std::pair<size_t, size_t>>> centroid_vector_sizes_private(centroid_vector.size());
            std::vector<size_t> group_indices_private;

            #pragma omp parallel for schedule(dynamic)
            for (size_t index = 0; index < group_indices.size(); index++)
            {
                // get entry information for the current ORF
                const size_t& ORF_ID = group_indices.at(index);
                const auto& ORF2_entry = ORF_mat_vector.at(ORF_ID);
                const auto& ORF2_info = colour_ORF_map.at(ORF2_entry.first).at(ORF2_entry.second);
                const auto& ORF2_len = std::get<2>(ORF2_info);

                // check if ORF encountered previously
                if (encountered_set.find(ORF_ID) != encountered_set.end())
                {
                    continue;
                }

                // create a set for all the centroids for the groups the ORF belongs to
                std::unordered_set<size_t> centroid_set;

                // go through the group_IDs
                for (const auto& group_ID : ORF_group_vector.at(ORF_ID))
                {
                    if ((bool)centroid_assigned[group_ID])
                    {
                        // add all centroids to centroid set for shared k-mer groups
                        centroid_set.insert(centroid_vector.at(group_ID));
                    }
                }

                // keep track of whether ORF itself is centroid
                bool is_centroid = false;

                // go through and align to each centroid, determine what is the highest scoring assignment
                size_t assigned_centroid = ORF_ID;
                double assigned_perc_id = 0;
                size_t centroid_len = 0;

                for (const auto& centroid : centroid_set)
                {
                    if (centroid != ORF_ID)
                    {
                        const auto& ORF1_entry = ORF_mat_vector.at(centroid);
                        const auto& ORF1_info = colour_ORF_map.at(ORF1_entry.first).at(ORF1_entry.second);
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
                            assigned_centroid = centroid;
                            centroid_len = ORF1_len;
                            break;
                        }
                    } else
                    {
                        is_centroid = true;
                    }
                }

                // check if centroid perc_id is greater than cut-off, if not then assign the cluster to it's own centroid
                if (assigned_perc_id)
                {
                    cluster_map[assigned_centroid].insert(ORF_ID);
                    cluster_map[assigned_centroid].insert(assigned_centroid);
                    encountered_set.insert(ORF_ID);

                    // assign ORF_ID and length, as well as for centroid, if not previously assigned
                    ORF_length_list_private.push_back({ORF2_len, ORF_ID});
                    if (encountered_set.find(assigned_centroid) == encountered_set.end())
                    {
                        ORF_length_list_private.push_back({centroid_len, assigned_centroid});
                        encountered_set.insert(assigned_centroid);
                    }
                }
                else
                {
                    // if is a singleton centroid, add to own cluster
                    if (is_centroid)
                    {
                        cluster_map[ORF_ID].insert(ORF_ID);
                        // assign ORF_ID and length
                        ORF_length_list_private.push_back({ORF2_len, ORF_ID});
                    } else
                    {
                        // if singleton and not centroid, add the current ORF as a centroid for all groups it is part of
                        // enabling searching in next iteration
                        for (const auto& group_ID : ORF_group_vector.at(ORF_ID))
                        {
                            // add ORF as new centroid of group
                            centroid_vector_sizes_private[group_ID].push_back({ORF2_len, ORF_ID});
                        }
                        // add entry to group_indices_private to re-cluster in next pass
                        group_indices_private.push_back(ORF_ID);
                    }

                    // add to encountered_set
                    encountered_set.insert(ORF_ID);
                }
            }

            #pragma omp critical
            {
                // update node_colour_vector with calculated colours
                ORF_length_list.insert(ORF_length_list.end(), std::make_move_iterator(ORF_length_list_private.begin()), std::make_move_iterator(ORF_length_list_private.end()));

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
            // clear centroid vector
            centroid_vector.clear();
            centroid_vector.resize(centroid_vector_sizes.size());
            centroid_assigned.reset();

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

                    for (int j = 1; j < group_entry.size(); j++)
                    {
                        if (group_entry[j].first < largest.first)
                        {
                            break;
                        } else if (group_entry[j].first == largest.first)
                        {
                            // if size is equal, then take lowest hash as centroid
                            if (ID_hash_map.at(group_entry[j].second) < ID_hash_map.at(largest.second))
                            {
                                largest = group_entry[j];
                            }
                        }
                    }

                    // assign ID to centroid vector
                    centroid_vector[i] = largest.second;
                    centroid_assigned[i] = 1;
                }
            }
        }
    }

    // sort the ORFs in descending order
    std::stable_sort(ORF_length_list.rbegin(), ORF_length_list.rend());

    // generate set to determine which ORFs have already been assigned clusters
    robin_hood::unordered_set<size_t> cluster_assigned;

    // iterate over ORF_length_list, pulling out centroids and their assigned clustered ORFs
    for (size_t i = 0; i < ORF_length_list.size(); i++)
    {
        const size_t& ORF_ID = ORF_length_list.at(i).second;

        if (cluster_assigned.find(ORF_ID) != cluster_assigned.end())
        {
            continue;
        }

        // check if entry is a centroid, if so add to final_clusters along with all of its attached ORFs that are unassigned
        if (cluster_map.find(ORF_ID) != cluster_map.end())
        {
            // add current ORF to centroid entry
            final_clusters[ORF_ID].push_back(ORF_ID);

            // add the centroid to the cluster_assigned set
            cluster_assigned.insert(ORF_ID);

            // add rest of homologs to centroid entry
            for (const auto& homolog_ID : cluster_map.at(ORF_ID))
            {
                // if the homolog is not already assigned to a cluster, assign and add to cluster_assigned
                if (cluster_assigned.find(homolog_ID) == cluster_assigned.end())
                {
                    final_clusters[ORF_ID].push_back(homolog_ID);
                    cluster_assigned.insert(homolog_ID);
                }
                // if homolog already assigned to a cluster, skip
                else
                {
                    continue;
                }

                // check whether homolog is also a centroid, and that it is not the current ORF being incremented
                if (cluster_map.find(homolog_ID) != cluster_map.end() && homolog_ID != ORF_ID)
                {
                    // if entry is centroid, pull all groups it is centroid for into final cluster find next largest entry and assign that as centroid
                    for (const auto& homolog2_ID : cluster_map.at(homolog_ID))
                    {
                        // if homolog2 cluster already assigned, skip
                        if (cluster_assigned.find(homolog2_ID) != cluster_assigned.end())
                        {
                            continue;
                        }

                        // else put it in the cluster
                        final_clusters[ORF_ID].push_back(homolog2_ID);
                        cluster_assigned.insert(homolog2_ID);

                        // check if this homolog is a centroid. If it is, go over each group it is centroid for and assign highest as centroid
                        if (cluster_map.find(homolog2_ID) != cluster_map.end() && homolog2_ID != homolog_ID)
                        {
                            size_t centroid_id = 0;
                            size_t centroid_len = 0;
                            tbb::concurrent_unordered_set<size_t> homolog_set;

                            for (const auto& homolog3_ID : cluster_map.at(homolog2_ID))
                            {
                                // if homolog3_ID already assigned, then pass
                                if (cluster_assigned.find(homolog3_ID) != cluster_assigned.end())
                                {
                                    continue;
                                }

                                // check that homolog is not the same as current iteration, and that this
                                // homolog is not in the current set for the current centroid
                                if (homolog3_ID != homolog2_ID && cluster_map.at(ORF_ID).find(homolog3_ID) == cluster_map.at(homolog_ID).end())
                                {
                                    const auto& ORF_entry = ORF_mat_vector.at(homolog3_ID);
                                    const auto& ORF_info = colour_ORF_map.at(ORF_entry.first).at(ORF_entry.second);
                                    const auto& homolog_len = std::get<2>(ORF_info);

                                    // assign new centroid ID if the homolog is longer or equal and has a lower ID number
                                    if (homolog_len > centroid_len || (homolog_len == centroid_len && homolog3_ID < centroid_id))
                                    {
                                        centroid_id = homolog3_ID;
                                        centroid_len = homolog_len;
                                    }

                                    homolog_set.insert(homolog3_ID);
                                }
                            }
                            // assign a new centroid with the unassigned ORFs
                            if (!homolog_set.empty())
                            {
                                // check if the new centroid already has entry in cluster_map, if so, append
                                if (cluster_map.find(centroid_id) != cluster_map.end())
                                {
                                    cluster_map[centroid_id].insert(homolog_set.begin(), homolog_set.end());
                                } else
                                {
                                    cluster_map[centroid_id] = homolog_set;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    auto num_ORFs = ORF_mat_vector.size();
    auto num_clusters = final_clusters.size();

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