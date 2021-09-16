//
// Created by sth19 on 01/09/2021.
//

#include "ORF_clustering.h"

std::tuple<ORFMatrixVector, std::vector<std::unordered_set<size_t>>, std::vector<std::pair<size_t, size_t>>> group_ORFs(const std::unordered_map<size_t, ORFNodeMap>& colour_ORF_map,
                                                               const GraphVector& graph_vector)
{
    // initilaise IndentityVector
    //IndentityVector id_vector;

    // calculate expected % matching for given id_cutoff and k-mer length
//    const double exp_matching = pow(id_cutoff, (double) (DBG_overlap + 1));

    // intialise Eigen Triplet
    std::vector<ET> tripletList;

    // generate a vector which maps ORF colour/ID to column in eigen matrix
    ORFMatrixVector ORF_mat_vector;

    // iterate over each ORF sequence with specific colours combination
    size_t ORF_ID = 0;
    for (const auto& colour : colour_ORF_map)
    {
        ORF_mat_vector.reserve(ORF_mat_vector.size() + colour.second.size());
        for (const auto& ORF_map : colour.second)
        {
            // add to ORF_mat_map, initialising empty vector
//            ORF_mat_vector.push_back({{colour.first, ORF_map.first}, {}});
            ORF_mat_vector.push_back({colour.first, ORF_map.first});

            // iterate over nodes traversed by ORF
            const auto& ORF_nodes = std::get<0>(ORF_map.second);
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
    std::vector<std::pair<size_t, size_t>> centroid_vector(graph_vector.size());

    // initialise sparse matrix
    Eigen::SparseMatrix<double> mat(ORF_mat_vector.size(), graph_vector.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());


    // conduct transposition + matrix multiplication to calculate ORFs sharing nodes
    //auto ORF_overlap_mat = ((mat * mat.transpose()).pruned()).eval();

    // iterate over non-zero entries in matrix and calculate overlaps
    for (int outit = 0; outit < mat.outerSize(); ++outit)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator init(mat, outit); init; ++init)
        {
//            // iterate through the bottom half of the symmetric matrix and ignore line of symmetry, pass if row <= col
//            if(init.row() <= init.col())
//            {
//                continue;
//            }

            // get column id (node) and row ID (ORF)
            const auto& node_id = init.col();
            const auto& ORF_id = init.row();

            // add the ORF to the group
            ORF_group_vector[ORF_id].insert(node_id);

            // get reference to entry in ORF_mat_vector
            const auto& ORF_entry = ORF_mat_vector.at(ORF_id);

            // get length of ORF entry
            const size_t& ORF_length = std::get<2>(colour_ORF_map.at(ORF_entry.first).at(ORF_entry.second));

            // determine the current size of the ORF and the centroid. Initialiser value will be 0 for unassigned centroid
            if (centroid_vector.at(node_id).second != 0)
            {
                auto& centroid = centroid_vector.at(node_id);

                // determine if current ORF is larger than current centroid
                if (ORF_length > centroid.second)
                {
                    centroid = {ORF_id, ORF_length};
                }
            } else
            {
                centroid_vector[node_id] = {ORF_id, ORF_length};
            }

//            // Assign temporary values for ORF1 and ORF2, not sorted by traversed node vector.
//            const auto& temp_ORF1_ID = ORF_mat_vector.at(init.col()).first;
//            const auto& temp_ORF2_ID = ORF_mat_vector.at(init.row()).first;
//
//
//            // Get nodes traversed by genes. Order by length of the node vector; ORF1 is the longer of the two vectors. Need to copy ORF2, as may be reversed
//            // check if ORF1 traverse more nodes than ORF2
//            const bool temp_ORF1_longer = ((std::get<2>(colour_ORF_map.at(temp_ORF1_ID.first).at(temp_ORF1_ID.second))
//                    >= std::get<2>(colour_ORF_map.at(temp_ORF2_ID.first).at(temp_ORF2_ID.second))) ? true : false);
//
//            // get respective ORF IDs
//            const std::pair<size_t, size_t>& ORF1_ID = (temp_ORF1_longer ? temp_ORF1_ID : temp_ORF2_ID);
//            const std::pair<size_t, size_t>& ORF2_ID = (temp_ORF1_longer ? temp_ORF2_ID : temp_ORF1_ID);
//
//            // assign matrix IDs for ORF_mat_vector
//            const auto& ORF1_mat_ID = (temp_ORF1_longer ? init.col() : init.row());
//            const auto& ORF2_mat_ID = (temp_ORF1_longer ? init.row() : init.col());
//
//            // get reference to ORF1_info information and unpack
//            const auto& ORF1_info = colour_ORF_map.at(ORF1_ID.first).at(ORF1_ID.second);
//            const auto& ORF1_node_ids = std::get<0>(ORF1_info);
//            const auto& ORF1_node_coords = std::get<1>(ORF1_info);
//            const auto& ORF1_len = std::get<2>(ORF1_info);
//
//            // unpack ORF2_info
//            const auto& ORF2_info = colour_ORF_map.at(ORF2_ID.first).at(ORF2_ID.second);
//            const auto& ORF2_node_ids = std::get<0>(ORF2_info);
//            const auto& ORF2_node_coords = std::get<1>(ORF2_info);
//            const auto& ORF2_len = std::get<2>(ORF2_info);
//
//
//            //testing
//            const auto ORF1_sequence = generate_sequence_private(std::get<0>(ORF1_info), std::get<1>(ORF1_info), DBG_overlap, graph_vector);
//            const auto ORF2_sequence = generate_sequence_private(std::get<0>(ORF2_info), std::get<1>(ORF2_info), DBG_overlap, graph_vector);
////
////            int test = 0;
////            if (ORF1_sequence == "TTGTTAGAAATCGATTTGACTGTCCTGATCGATTTGTCATGTTCTTATTTCATTTTACTATATTTTTGGTTCGCGGGAAGTCTACTAAGATACTTAAAGATGCAGATAGTAAAAAAATGTAGACATTACCGTAAAAAAGTGATATAA"
////            && ORF2_sequence == "ATGTTCTTATTTCATTTTACTATATTTTTGTTTCGCGGGAAGTCTACTAAGATACTTAAAGATGCAGATAGTAAAAAAAATGTAGACATTACCGTAAAAAAGTGA")
////            {
////                test = 1;
////            }
//
//
//            // check that perc_len_diff is greater than cut-off, otherwise pass
//            double perc_len_diff = (double)ORF2_len / (double)ORF1_len;
//
//            if (perc_len_diff < len_diff_cutoff)
//            {
//                continue;
//            }
//
//            // make sets of ORF2_node_IDs and ORF1_node_IDs and calculate % identity
//            std::set<int> ORF1_node_set(ORF1_node_ids.begin(), ORF1_node_ids.end());
//            std::set<int> ORF2_node_set(ORF2_node_ids.begin(), ORF2_node_ids.end());
//
//            std::vector<int> set_union;
//            std::set_intersection(ORF1_node_set.begin(), ORF1_node_set.end(), ORF2_node_set.begin(), ORF2_node_set.end(),
//                                  std::back_inserter(set_union));
//
//            // if set_union is empty, ORFs share no nodes, so pass
//            if (set_union.empty())
//            {
//                continue;
//            }
//
//            // go through set_union, identify the nodes in each ORF set and determine overlap
//            size_t abs_overlap = 0;
//            for (const auto& node : set_union)
//            {
//                // get index of overlapping node in ORF1_node_ids/ORF2_node_ids
//                auto ORF1_it = std::find(ORF1_node_ids.begin(), ORF1_node_ids.end(), node);
//                auto ORF2_it = std::find(ORF2_node_ids.begin(), ORF2_node_ids.end(), node);
//                size_t ORF1_index = ORF1_it - ORF1_node_ids.begin();
//                size_t ORF2_index = ORF2_it - ORF2_node_ids.begin();
//
//                // compare coordinates, calculate overlap
//                const auto& ORF1_coords = ORF1_node_coords.at(ORF1_index);
//                const auto& ORF2_coords = ORF2_node_coords.at(ORF2_index);
//                const size_t& ORF1_start = std::get<0>(ORF1_coords);
//                const size_t& ORF1_end = std::get<1>(ORF1_coords);
//                const size_t& ORF2_start = std::get<0>(ORF2_coords);
//                const size_t& ORF2_end = std::get<1>(ORF2_coords);
//
//                const auto node_seq = graph_vector.at(abs(node) - 1).seq();
//
//                // check there is an intersection
//                size_t node_overlap = 0;
//                if ((ORF1_start <= ORF2_end) && (ORF2_start <= ORF1_end))
//                {
//                    size_t low_index = std::max(ORF1_start, ORF2_start);
//                    size_t high_index = std::min(ORF1_end, ORF2_end);
//                    node_overlap += (high_index - low_index) + 1;
//
//                    // need to check if previous node is in set, if so, then determine overlapping region from previous node
//                    if (ORF1_index != 0 && ORF2_index != 0)
//                    {
//                        // see if node prior to this node is present in the overlap
//                        if (ORF1_node_ids.at(ORF1_index - 1) == ORF2_node_ids.at(ORF2_index - 1))
//                        {
//                            size_t prev_overlap = (DBG_overlap - low_index);
//                            //if index
//                            if (high_index >= DBG_overlap)
//                            {
//                                // negate overlap from previous node
//                                node_overlap -= prev_overlap;
//                            } else
//                            {
//                                node_overlap = 0;
//                            }
//                        }
//                        // need to think about insertions that will not allow easy calculation of node overlap
//                        else if (ORF1_index != 1 && ORF2_index != 1)
//                        {
//                            // look two nodes back. If insertion present, then node in between will be less than 2 * DBG_overlap + 1
//                            if (ORF1_node_ids.at(ORF1_index - 2) == ORF2_node_ids.at(ORF2_index - 2))
//                            {
//                                // get lengths of previous nodes if they don't match
//                                const auto& ORF1_prev_len = graph_vector.at(abs(ORF1_node_ids.at(ORF1_index - 1)) - 1).size().first;
//                                const auto& ORF2_prev_len = graph_vector.at(abs(ORF2_node_ids.at(ORF2_index - 1)) - 1).size().first;
//
//                                // calculate the smaller of the two
//                                size_t min_len = std::min(ORF1_prev_len, ORF2_prev_len);
//
//                                // if the length is less than 2* the DBG overlap + 1, then need to account for prev node overlap
//                                if (min_len < 2 * DBG_overlap)
//                                {
//                                    size_t prev_overlap = ((2 * DBG_overlap) - min_len);
//                                    if (high_index >= DBG_overlap)
//                                    {
//                                        // negate overlap from previous node
//                                        node_overlap -= prev_overlap;
//                                    } else
//                                    {
//                                        node_overlap = 0;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                abs_overlap += node_overlap;
//            }
//
//            // initialise bool to determine if perc_identity is accurate
//            bool id_accurate = false;
//
//            // calculate percentage identity based on k-mers
//            double perc_id = (double)abs_overlap / (double)ORF2_len;
//
////            if (perc_id < 1 && abs_overlap >= DBG_overlap + 1)
////            {
////                size_t kmers_matching = abs_overlap - DBG_overlap + 2;
////                double perc_matching = (double)kmers_matching / (double)(ORF2_len - DBG_overlap + 2);
////
////                if (perc_matching >= exp_matching)
////                {
////                    perc_id = align_seqs(ORF1_info, ORF2_info, graph_vector, DBG_overlap);
////                    id_accurate = true;
////                }
////            } else
////            {
////                id_accurate = true;
////            }
//
//            // if perc_id / kmer length is greater than expected matching but less than 1, need to align to
//            // calculate perc_id accurately. Otherwise can ignore.
//
//
//            int test2 = 0;
//            if (perc_id > 1)
//            {
//                test2 = 1;
//            }
//
//            // if cut-offs satisfied and perc_id is accurate, add to id_vector
////            if (perc_id >= id_cutoff && id_accurate)
////            {
//                ORF_mat_vector[ORF1_mat_ID].second[ORF2_mat_ID] = {perc_id, perc_len_diff};
//                ORF_mat_vector[ORF2_mat_ID].second[ORF1_mat_ID] = {perc_id, perc_len_diff};
//            }
        }
    }

    //std::pair<IndentityVector, ORFMatrixVector> ORF_cluster_pair = std::make_pair(id_vector, ORF_mat_vector);
    //return ORF_cluster_pair;
    auto return_tuple = std::make_tuple(ORF_mat_vector, ORF_group_vector, centroid_vector);
    return return_tuple;
}

std::unordered_map<size_t, std::vector<size_t>> produce_clusters(const std::unordered_map<size_t, ORFNodeMap>& colour_ORF_map,
                               const GraphVector& graph_vector,
                               const size_t& DBG_overlap,
                               const ORFMatrixVector& ORF_mat_vector,
                               const std::vector<std::unordered_set<size_t>>& ORF_group_vector,
                               const std::vector<std::pair<size_t, size_t>>& centroid_vector,
                               const double& id_cutoff,
                               const double& len_diff_cutoff)
{
    // initialise temporary cluster map, first item is centroid, second is set of nodes in that cluster
    std::unordered_map<size_t, std::unordered_set<size_t>> cluster_map;

    // initialise final cluster map for return
    std::unordered_map<size_t, std::vector<size_t>> final_clusters;

    // initialise length list for sorting of ORFs by length
    std::vector<std::pair<size_t, size_t>> ORF_length_list;
    ORF_length_list.reserve(ORF_group_vector.size());

    // go through ORF_mat_vector, assign each ORF to a group is any of its homologues are present in a group, or a new group if non of its current
    // homologues are in a group. This can be done for all the homologues in the vector.
    // then need to identify the longest ORF in each group, assign this as a single group (merging all groups it is part of
    // and then assign each ORF to the group it has highest identity with (need to keep identity measure to to this!)
//    std::unordered_set<size_t> ORF_index_set;
//
//    // mapping of node IDs to the groups it belongs to
//    std::unordered_map<size_t, std::unordered_set<size_t>> ORF_to_group_map;
//
//
//    // mapping of group to ORFs contained within
//    std::unordered_map<size_t, std::unordered_set<size_t>> group_map;

    // create set to keep track of assigned ORFs to avoid ORFs being added to other groups
    std::unordered_set<size_t> encountered_set;

    int number_alignments = 0;

    // iterate through group_map
    for (size_t ORF_ID = 0; ORF_ID < ORF_group_vector.size(); ORF_ID++)
    {
        // get entry information for the current ORF
        const auto& ORF2_entry = ORF_mat_vector.at(ORF_ID);
        const auto& ORF2_info = colour_ORF_map.at(ORF2_entry.first).at(ORF2_entry.second);
        const auto& ORF2_len = std::get<2>(ORF2_info);

        // assign ORF_ID and length
        ORF_length_list.push_back({ORF2_len, ORF_ID});

        // check if ORF encountered previously
        if (encountered_set.find(ORF_ID) != encountered_set.end())
        {
            continue;
        }

        // create a set for all the centroids for the groups the ORF belongs to
        std::unordered_set<size_t> centroid_set;

        // create vector that holds identical sequences
        std::vector<size_t> identical_ORFs = {ORF_ID};

        // go through the group_IDs
        for (const auto& group_ID : ORF_group_vector.at(ORF_ID))
        {
            centroid_set.insert(centroid_vector.at(group_ID).first);
        }


        // go through and align to each centroid, determine what is the highest scoring assignment
        size_t assigned_centroid = ORF_ID;
        double assigned_perc_id = 0;
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
                double current_perc_id = align_seqs(ORF1_info, ORF2_info, graph_vector, DBG_overlap);
                number_alignments++;

                if (current_perc_id > assigned_perc_id)
                {
                    assigned_perc_id = current_perc_id;
                    assigned_centroid = centroid;
                }
                // when current_perc_id == assigned_perc_id, assign assigned_centroid to lowest index
                else if (current_perc_id == assigned_perc_id)
                {
                    assigned_centroid = std::min(centroid, assigned_centroid);
                }

                // if the two sequences are identical add to identical_ORFs
                if (current_perc_id == 1)
                {
                    identical_ORFs.push_back(centroid);
                }
            }
        }

        // chec if there are identical ORFs present
        if (identical_ORFs.size() != 1)
        {
            // add all entries to first entry key in cluster_map
            std::sort(identical_ORFs.begin(), identical_ORFs.end());
            cluster_map[identical_ORFs[0]].insert(identical_ORFs.begin(), identical_ORFs.end());
            encountered_set.insert(identical_ORFs.begin(), identical_ORFs.end());
        }
        // check if centroid perc_id is greater than cut-off, if not then assign the cluster to it's own centroid
        else if (assigned_perc_id >= id_cutoff)
        {
            cluster_map[assigned_centroid].insert(ORF_ID);
            cluster_map[assigned_centroid].insert(assigned_centroid);
            encountered_set.insert(ORF_ID);
            encountered_set.insert(assigned_centroid);
        }
        else
        {
            cluster_map[ORF_ID].insert(ORF_ID);
            encountered_set.insert(ORF_ID);
        }
    }

//    for (size_t ORF_index = 0; ORF_index < ORF_mat_vector.size(); ORF_index++)
//    {
//        // check if ORF_index already assigned to group
////        if (ORF_index_set.find(ORF_index) != ORF_index_set.end())
////        {
////            continue;
////        }
//
//        // create unordered_set for groups for current ORF
//        std::unordered_set<size_t> ORF_group_set;
//
//        // add current entry to ORF_group_set
//        ORF_group_set.insert(ORF_index);
//
//        // get reference to entry in ORF_mat_vector
//        const auto& ORF_entry = ORF_mat_vector.at(ORF_index);
//
//        // get length of ORF entry
//        size_t centroid_length = std::get<2>(colour_ORF_map.at(ORF_entry.first.first).at(ORF_entry.first.second));
//
//        // set longest ORF in comparison, this will be used as the centre ORF
//        size_t centroid = ORF_index;
//
//        // go through ORF_entry.second, check if any other ORFs have been assigned to a group. If not, then create new group
//        bool homolog_present = false;
//        for (const auto& homolog : ORF_entry.second)
//        {
//            const auto& homolog_ID = homolog.first;
//
//            // check if homolog has already been added to a group
////            if (ORF_index_set.find(homolog_ID) != ORF_index_set.end())
////            {
////                continue;
////            }
//
//            // multiply identity and perc_len_diff to get identity score
//            double identity_score = homolog.second.first * homolog.second.second;
//
//            const auto& homolog_entry = ORF_mat_vector.at(homolog_ID);
//
//            const size_t& ORF_length = std::get<2>(colour_ORF_map.at(homolog_entry.first.first).at(homolog_entry.first.second));
//
//            if (ORF_length > centroid_length)
//            {
//                centroid = homolog_ID;
//                centroid_length = ORF_length;
//            }
//
//            // add to the group set
//            ORF_group_set.insert(homolog_ID);
////            ORF_index_set.insert(homolog_ID);
//
//            // keep track of which groups the current ORF has been added to
//            ORF_to_group_map[homolog_ID].insert(group_ID)
//
////            if (identity_score < 1)
////            {
////                // check entry for homolog to see if there are any higher scoring pairs, if so ignore
////                // get reference to entry in ORF_mat_vector
////                const auto& homolog_entry = ORF_mat_vector.at(homolog_ID);
////
////                bool not_in_group = false;
////
////                for (const auto& homolog2 : homolog_entry.second)
////                {
////                    double identity_score2 = homolog2.second.first * homolog2.second.second;
////
////                    if (identity_score2 > identity_score)
////                    {
////                        not_in_group = true;
////                        break;
////                    }
////                }
////
////                // if no higher scoring homolog found, add to the group
////                if (!not_in_group)
////                {
////                    ORF_group_set.insert(homolog_ID);
////                    ORF_index_set.insert(homolog_ID);
////                }
////
////            } else
////            {
////                ORF_group_set.insert(homolog_ID);
////                ORF_index_set.insert(homolog_ID);
////            }
//
//
////            // go through and add the groups other homologs belong to
////            if (ORF_2_group_map.find(homolog_ID) != ORF_2_group_map.end())
////            {
////                homolog_present = true;
////                // merge ORF_group_set for current ORF and homolog
////                ORF_group_set.insert(ORF_2_group_map.at(homolog_ID).begin(), ORF_2_group_map.at(homolog_ID).end());
////            }
//            centroid_map[group_ID] = {centroid, centroid_length};
//        }
//
//        ORF_index_set.insert(ORF_index);
//
//        group_map[group_ID++] = std::move(ORF_group_set);
//
////        // if no homolog present, create new group and add to it and iterate group_ID
////        if (!homolog_present)
////        {
////            ORF_group_set.insert(group_ID++);
////        }
////
////        // for each group in ORF_group_set, determine if current ORF is the longest and update group_2_centre_map
////        for (const auto& group : ORF_group_set)
////        {
////            if (group_2_centre_map.find(group) != group_2_centre_map.end())
////            {
////                auto& current_ORF = group_2_centre_map[group];
////                // if the new ORF is longer than the current ORF, the update entry
////                if (ORF_length > current_ORF.second)
////                {
////                    current_ORF = {ORF_index, ORF_length};
////                }
////            } else
////            {
////                // if group not present, then add to map
////                group_2_centre_map[group] = {ORF_index, ORF_length};
////            }
////        }
////
////        ORF_2_group_map[ORF_index] = std::move(ORF_group_set);
//    }

    // sort the ORFs in descending order
    sort(ORF_length_list.rbegin(), ORF_length_list.rend());

    // generate set to determine which ORFs have already been assigned clusters
    std::set<size_t> cluster_assigned;

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

                // check whether homolog is also a centroid, and that it is not the current ORF being incremented
                if (cluster_map.find(homolog_ID) != cluster_map.end() && homolog_ID != ORF_ID)
                {
                    // if entry is centroid, find next largest entry and assign that as centroid
                    size_t centroid_id;
                    size_t centroid_length = 0;
                    std::unordered_set<size_t> homolog_set;

                    for (const auto& homolog2_ID : cluster_map.at(homolog_ID))
                    {
                        // check that homolog is not the same as current iteration, and that this
                        // homolog is not in the current set for the current centroid
                        if (homolog2_ID != homolog_ID && cluster_map.at(ORF_ID).find(homolog2_ID) == cluster_map.at(ORF_ID).end())
                        {
                            const auto& ORF_entry = ORF_mat_vector.at(homolog2_ID);
                            const auto& ORF_info = colour_ORF_map.at(ORF_entry.first).at(ORF_entry.second);
                            const auto& homolog_len = std::get<2>(ORF_info);

                            if (homolog_len > centroid_length)
                            {
                                centroid_id = homolog2_ID;
                                centroid_length = homolog_len;
                            }

                            homolog_set.insert(homolog2_ID);
                        }
                    }

                    // assign a new centroid with the unassigned ORFs
                    if (!homolog_set.empty())
                    {
                        cluster_map[centroid_id] = homolog_set;
                    }
                }
            }
        }
    }

    // sanity check to ensure no ORFs are being assigned to more than one cluster
    std::unordered_set<size_t> check_set1;
    int test2 = 0;
    for (const auto& group : cluster_map)
    {
        for (const auto& node : group.second)
        {
            if (check_set1.find(node) != check_set1.end())
            {
                test2 = 1;
            } else
            {
                check_set1.insert(node);
            }
        }
    }
    std::unordered_set<size_t> check_set2;
    test2 = 0;
    for (const auto& group : final_clusters)
    {
        for (const auto& node : group.second)
        {
            if (check_set2.find(node) != check_set2.end())
            {
                test2 = 1;
            } else
            {
                check_set2.insert(node);
            }
        }
    }

    // check all ORFs are accounted for
    size_t all_ORFs = ORF_mat_vector.size();
    size_t grouped_ORFs_pre = check_set1.size();
    size_t grouped_ORFs_post = check_set2.size();

    return final_clusters;
}

std::string generate_sequence_private(const std::vector<int>& nodelist,
                                      const std::vector<indexPair>& node_coords,
                                      const size_t& overlap,
                                      const GraphVector& _GraphVector)
{
    std::string sequence;
    for (size_t i = 0; i < nodelist.size(); i++)
    {
        // initialise sequence items
        std::string unitig_seq;
        std::string substring;

        // parse information
        const auto& id = nodelist[i];
        const auto& coords = node_coords[i];
        bool strand = (id >= 0) ? true : false;

        if (strand)
        {
            unitig_seq = _GraphVector.at(abs(id) - 1).seq();
        } else {
            unitig_seq = reverse_complement(_GraphVector.at(abs(id) - 1).seq());
        }

        if (sequence.empty())
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - std::get<0>(coords)) + 1;
            substring = unitig_seq.substr(std::get<0>(coords), node_seq_len);
        } else
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - overlap) + 1;
            // need to account for overlap, if overlap is greater than the end of the node, sequence already accounted for
            if (node_seq_len > 0)
            {
                substring = unitig_seq.substr(overlap, node_seq_len);
            } else
            {
                break;
            }
        }
        sequence += substring;
    }
    return sequence;
}

double align_seqs(const ORFNodeVector& ORF1_info,
                  const ORFNodeVector& ORF2_info,
                  const GraphVector& graph_vector,
                  const size_t& DBG_overlap)
{
    std::string ORF1_aa;
    std::string ORF2_aa;

    // scope for DNA sequences
    {
        const auto ORF1_sequence = generate_sequence_private(std::get<0>(ORF1_info), std::get<1>(ORF1_info), DBG_overlap, graph_vector);
        const auto ORF2_sequence = generate_sequence_private(std::get<0>(ORF2_info), std::get<1>(ORF2_info), DBG_overlap, graph_vector);

        // reserve memory for amino acids
        ORF1_aa.reserve(ORF1_sequence.size() / 3);
        ORF2_aa.reserve(ORF2_sequence.size() / 3);

        const auto ORF1_dna_vect = ORF1_sequence | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;
        const auto ORF2_dna_vect = ORF2_sequence | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;

        const auto ORF1_aa_vect = ORF1_dna_vect | seqan3::views::translate_single(seqan3::translation_frames::FWD_FRAME_0);
        const auto ORF2_aa_vect = ORF2_dna_vect | seqan3::views::translate_single(seqan3::translation_frames::FWD_FRAME_0);

        for (auto && residue : ORF1_aa_vect)
        {
            ORF1_aa += residue.to_char();
        }

        for (auto && residue : ORF2_aa_vect)
        {
            ORF2_aa += residue.to_char();
        }
    }

    // Configure the alignment kernel, using hamming distance
    // and -1 penalisation for gaps. Enables easy calculation of perc identity
//    auto config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
//                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
//                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
//                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}} |
//                  seqan3::align_cfg::scoring_scheme{
//                          seqan3::aminoacid_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{0}}} |
//                  seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{0}, seqan3::gap_open_score{0}}};

//    auto config = seqan3::align_cfg::method_global{} |
//                  seqan3::align_cfg::scoring_scheme{
//                          seqan3::aminoacid_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{0}}} |
//                  seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{0}, seqan3::gap_open_score{0}}};

//    auto config = seqan3::align_cfg::method_local |
//                  seqan3::align_cfg::scoring_scheme{
//                          seqan3::aminoacid_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{0}}} |
//                  seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{0}, seqan3::gap_open_score{0}}};

    // calculate score of alignment
//    auto results = seqan3::align_pairwise(std::tie(ORF1_aa, ORF2_aa), config);
//    auto & res = *results.begin();

//    double perc_id = (double)res.score() / (double)ORF1_aa.size();

    const char * seq1 = ORF1_aa.c_str();
    const char * seq2 = ORF2_aa.c_str();
    EdlibAlignResult result = edlibAlign(seq1, ORF1_aa.size(), seq2, ORF1_aa.size(), edlibDefaultAlignConfig());
    size_t edit_distance = result.editDistance;
    double perc_id = 1 - ((double)edit_distance / (double)ORF2_aa.size());
    return perc_id;
}