#include "ORF_clustering.h"

ORFGroupPair group_ORFs(const std::map<size_t, std::string>& ORF_file_paths,
                        const ColoredCDBG<MyUnitigMap>& ccdbg,
                        const std::vector<Kmer>& head_kmer_arr,
                        const size_t& overlap)
{
    // intialise Eigen Triplet
    std::vector<ET> tripletList;

    // initialise length list for sorting of ORFs by length
    std::vector<std::pair<size_t, std::pair<size_t, size_t>>> ORF_length_list;

    // initialise vector to hold all centroid IDs in, and dynamic bitset to determine which centroids have been assigned
    std::shared_ptr<std::string> placeholder = std::make_shared<std::string>("NA");
    std::vector<std::tuple<int, int, size_t, size_t, std::shared_ptr<std::string>>> centroid_vector(head_kmer_arr.size(), {-1, -1, 0, 0, placeholder});

    // iterate over each ORF sequence with specific colours combination
    for (const auto& colour : ORF_file_paths)
    {
        ORFNodeMap ORF_map;
        // read in ORF_map file
        {
            std::ifstream ifs(colour.second);
            boost::archive::text_iarchive ia(ifs);
            ia >> ORF_map;
        }
        
        for (const auto& ORF_entry : ORF_map)
        {
            // extract info from ORF_entry
            const auto& ORF_info = ORF_entry.second;

            // iterate over nodes traversed by ORF
            const auto& ORF_nodes = std::get<0>(ORF_info);

            ORF_length_list.push_back({std::get<2>(ORF_info), {colour.first, ORF_entry.first}});

            // hold indices that are complete unitigs i.e. start/end not in overlap
            std::vector<size_t> complete_nodes;

            for (int i = 0; i < ORF_nodes.size(); i++)
            {
                const auto& node_indices = std::get<1>(ORF_info).at(i);
                // determine if at end and in overlap region of unitigs, if so pass
                if ((node_indices.second - node_indices.first) < overlap)
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
                complete_nodes.push_back(abs(ORF_nodes.at(i)) - 1);
            }
            // get hash to hash map
            std::string ORF_seq = generate_sequence_nm(std::get<0>(ORF_info), std::get<1>(ORF_info), overlap, ccdbg, head_kmer_arr);
            size_t ORF_hash = hasher{}(ORF_seq);
            // generate compact amino acid sequence
            ORF_seq = translate(ORF_seq);

            // get length of ORF entry
            const size_t& ORF_length = std::get<2>(ORF_info);

            // go over again, this time determining centroid status
            for (const auto& node_ID : complete_nodes)
            {
                // determine the current size of the ORF and the centroid if previous assigned
                if (std::get<0>(centroid_vector.at(node_ID)) != -1)
                {
                    auto& centroid_ID = centroid_vector[node_ID];
                    const auto& centroid_length = std::get<3>(centroid_ID);

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
                                centroid_ID = {colour.first, ORF_entry.first, ORF_hash, ORF_length, std::make_shared<std::string>(ORF_seq)};
                            } else if (ORF_hash == centroid_hash)
                            {
                                if (colour.first < std::get<0>(centroid_ID))
                                {
                                    centroid_ID = {colour.first, ORF_entry.first, ORF_hash, ORF_length, std::make_shared<std::string>(ORF_seq)};
                                }
                            }
                        } else
                        {
                            centroid_ID = {colour.first, ORF_entry.first, ORF_hash, ORF_length, std::make_shared<std::string>(ORF_seq)};
                        }
                    }
                } else
                {
                    centroid_vector[node_ID] = {colour.first, ORF_entry.first, ORF_hash, ORF_length, std::make_shared<std::string>(ORF_seq)};
                }
            }
        }
    }

    // sort the ORFs in descending order
    std::stable_sort(ORF_length_list.rbegin(), ORF_length_list.rend());

    //return group information
    auto return_pair = std::pair(centroid_vector, ORF_length_list);
    return return_pair;
}

// why not store all centroid sequences in a list, then iterate through all genes as before
// compare them to all centroids and calculate pairwise scores, cluster greedily
ORFClusterMap produce_clusters(const std::map<size_t, std::string>& ORF_file_paths,
                               const ColoredCDBG<MyUnitigMap>& ccdbg,
                               const std::vector<Kmer>& head_kmer_arr,
                               const size_t& overlap,
                               ORFGroupPair& ORF_group_pair,
                               const double& id_cutoff,
                               const double& len_diff_cutoff)
{
    // unpack ORF_group_pair
    auto& centroid_vector = ORF_group_pair.first;
    auto& ORF_length_list = ORF_group_pair.second;

    robin_hood::unordered_map<std::string, std::vector<std::pair<size_t, size_t>>> CentroidToORFMap;


    // iterate over each ORF sequence with specific colours combination
    for (const auto& colour : ORF_file_paths)
    {
        ORFNodeMap ORF_map;
        // read in ORF_map file
        {
            std::ifstream ifs(colour.second);
            boost::archive::text_iarchive ia(ifs);
            ia >> ORF_map;
        }
        
        for (const auto& ORF_entry : ORF_map)
        {
            const auto& ORF_info = ORF_entry.second;
            std::string ORF_seq = translate(generate_sequence_nm(std::get<0>(ORF_info), std::get<1>(ORF_info), overlap, ccdbg, head_kmer_arr));

            const auto& ORF_nodes = std::get<0>(ORF_info);
            const auto& ORF_len = std::get<2>(ORF_info);

            bool is_centroid = false;

            for (int i = 0; i < std::get<0>(ORF_info).size(); i++)
            {
                const auto& node_indices = std::get<1>(ORF_info).at(i);
                // determine if at end and in overlap region of unitigs, if so pass
                if ((node_indices.second - node_indices.first) < overlap)
                {
                    // if ORF doesn't span with node overlap, ignore
                    continue;
                }
                // convert node_traversed to size_t, minus 1 as unitigs are one-based, needs to be zero based
                const size_t current_node = abs(ORF_nodes.at(i)) - 1;

                const auto& centroid_ID = centroid_vector.at(current_node);
                std::string ORF_ID_string = std::to_string(colour.first) + "_" + std::to_string(ORF_entry.first);

                // check if entry is centroid itself, if already set then ignore as already added
                if (std::get<0>(centroid_ID) == colour.first && std::get<1>(centroid_ID) == ORF_entry.first)
                {
                    if (!is_centroid)
                    {
                        CentroidToORFMap[ORF_ID_string].push_back({colour.first, ORF_entry.first});
                        is_centroid = true;
                    }
                } else 
                {
                    const auto& centroid_len = std::get<3>(centroid_ID);

                    // check that perc_len_diff is greater than cut-off, otherwise pass
                    double perc_len_diff = (double)ORF_len / (double)centroid_len;

                    if (perc_len_diff < len_diff_cutoff)
                    {
                        continue;
                    }

                    const auto centroid_seq = *std::get<4>(centroid_ID);

                    // calculate the perc_id between the current centroid and the ORF
                    double current_perc_id = align_seqs(centroid_seq, ORF_seq);

                    // assign to all centroids with current_perc_id >= id_cutoff

                    if (current_perc_id >= id_cutoff)
                    {
                        std::string Centroid_ID_string = std::to_string(std::get<0>(centroid_ID)) + "_" + std::to_string(std::get<1>(centroid_ID));
                        CentroidToORFMap[Centroid_ID_string].push_back({colour.first, ORF_entry.first});
                    }
                }
            }
        }
    }

    // generate set to determine which ORFs have already been assigned clusters
    robin_hood::unordered_set<std::string> cluster_assigned;

    // initialise map as intermediate to hold cluster IDs
    ORFClusterMap final_clusters;

    // iterate over ORF_length_list, pulling out centroids and their assigned clustered ORFs
    for (size_t i = 0; i < ORF_length_list.size(); i++)
    {
        const auto& ORF_ID = ORF_length_list.at(i).second;

        std::string ORF_ID_str = std::to_string(ORF_ID.first) + "_" + std::to_string(ORF_ID.second);

        if (cluster_assigned.find(ORF_ID_str) != cluster_assigned.end())
        {
            continue;
        }

        // check if entry is a centroid, if so add to final_clusters along with all of its attached ORFs that are unassigned
        if (CentroidToORFMap.find(ORF_ID_str) != CentroidToORFMap.end())
        {
            // add current ORF to centroid entry
            final_clusters[ORF_ID_str].push_back(ORF_ID);

            // add the centroid to the cluster_assigned set
            cluster_assigned.insert(ORF_ID_str);

            // add rest of homologs to centroid entry
            for (const auto& homolog_ID : CentroidToORFMap.at(ORF_ID_str))
            {
                std::string homolog_ID_str = std::to_string(homolog_ID.first) + "_" + std::to_string(homolog_ID.second);

                // if the homolog is not already assigned to a cluster, assign and add to cluster_assigned
                if (cluster_assigned.find(homolog_ID_str) == cluster_assigned.end())
                {
                    final_clusters[ORF_ID_str].push_back(homolog_ID);
                    cluster_assigned.insert(homolog_ID_str);
                }
                // if homolog already assigned to a cluster, skip
                else
                {
                    continue;
                }
            }
            // erase from cluster_map
            CentroidToORFMap.erase(ORF_ID_str);
        } else
        {
            // if not assigned and not a centroid, then add to it's own cluster as singleton
            final_clusters[ORF_ID_str].push_back(ORF_ID);

            cluster_assigned.insert(ORF_ID_str);
        }
    }

    return final_clusters;
}

double align_seqs(const std::string& ORF1_aa,
                  const std::string& ORF2_aa)
{
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