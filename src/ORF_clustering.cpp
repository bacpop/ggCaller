//
// Created by sth19 on 01/09/2021.
//

#include "ORF_clustering.h"

ORFGroupTuple group_ORFs(const ColourORFMap& colour_ORF_map,
                         const GraphVector& graph_vector)
{
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
    std::vector<std::vector<std::pair<size_t, size_t>>> centroid_vector(graph_vector.size());

    // initialise sparse matrix
    Eigen::SparseMatrix<double> mat(ORF_mat_vector.size(), graph_vector.size());
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

            // determine the current size of the ORF and the centroid. Initialiser value will be 0 for unassigned centroid
            if (!centroid_vector.at(node_ID).empty())
            {
                auto& centroid = centroid_vector[node_ID][0];

                // determine if current ORF is larger than current centroid
                // if equal, then add lowest ORF index of two as centroid
                if (ORF_length > centroid.second || (ORF_length == centroid.second && ORF_ID < centroid.first))
                {
                    centroid = {ORF_ID, ORF_length};
                }
            } else
            {
                centroid_vector[node_ID].push_back({ORF_ID, ORF_length});
            }
        }
    }

    //return group information
    auto return_tuple = std::make_tuple(ORF_mat_vector, ORF_group_vector, centroid_vector);
    return return_tuple;
}

ORFClusterMap produce_clusters(const ColourORFMap& colour_ORF_map,
                               const GraphVector& graph_vector,
                               const size_t& DBG_overlap,
                               const ORFMatrixVector& ORF_mat_vector,
                               const std::vector<std::unordered_set<size_t>>& ORF_group_vector,
                               std::vector<std::vector<std::pair<size_t, size_t>>>& centroid_vector,
                               const double& id_cutoff,
                               const double& len_diff_cutoff)
{
    // initialise temporary cluster map, first item is centroid, second is set of nodes in that cluster
    std::unordered_map<size_t, std::unordered_set<size_t>> cluster_map;

    // initialise final cluster map for return
    ORFClusterMap final_clusters;

    // initialise length list for sorting of ORFs by length
    std::vector<std::pair<size_t, size_t>> ORF_length_list;
    ORF_length_list.reserve(ORF_group_vector.size());

    // initialise mapping from group_id to potential singleton ORFs
    std::unordered_map<size_t, std::unordered_set<size_t>> group_singleton_map;

    // create set to keep track of assigned ORFs to avoid ORFs being added to other groups
    std::unordered_set<size_t> encountered_set;

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

        // go through the group_IDs
        for (const auto& group_ID : ORF_group_vector.at(ORF_ID))
        {
            // iterate over all centroids for the given k-mer group
            for (const auto& centroid_ID : centroid_vector.at(group_ID))
            {
                centroid_set.insert(centroid_ID.first);
            }
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
            }
        }

        // check if centroid perc_id is greater than cut-off, if not then assign the cluster to it's own centroid
        if (assigned_perc_id >= id_cutoff)
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

            // if singleton, add the current ORF as a centroid for all groups it is part of
            for (const auto& group_ID : ORF_group_vector.at(ORF_ID))
            {
                // check that ORF is not already the centroid for the group
                if (centroid_vector.at(group_ID).at(0).first != ORF_ID)
                {
                    centroid_vector[group_ID].push_back({ORF_ID, ORF2_len});
                }
            }
        }
    }

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
                    size_t centroid_id = 0;
                    size_t centroid_len = 0;
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

                            // assign new centroid ID if the homolog is longer or equal and has a lower ID number
                            if (homolog_len > centroid_len || (homolog_len == centroid_len && homolog2_ID < centroid_id))
                            {
                                centroid_id = homolog2_ID;
                                centroid_len = homolog_len;
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

        ORF1_aa = (translate(ORF1_sequence)).aa();
        ORF2_aa = (translate(ORF2_sequence)).aa();

//        const auto ORF1_dna_vect = ORF1_sequence | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;
//        const auto ORF2_dna_vect = ORF2_sequence | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;
//
//        const auto ORF1_aa_vect = ORF1_dna_vect | seqan3::views::translate_single(seqan3::translation_frames::forward_frame0);
//        const auto ORF2_aa_vect = ORF2_dna_vect | seqan3::views::translate_single(seqan3::translation_frames::forward_frame0);
//
//        for (auto && residue : ORF1_aa_vect)
//        {
//            ORF1_aa += residue.to_char();
//        }
//
//        for (auto && residue : ORF2_aa_vect)
//        {
//            ORF2_aa += residue.to_char();
//        }
    }

    // convert sequence to const char * and calculate edit distance
    const char * seq1 = ORF1_aa.c_str();
    const char * seq2 = ORF2_aa.c_str();
    EdlibAlignResult result = edlibAlign(seq1, ORF1_aa.size(), seq2, ORF1_aa.size(), edlibDefaultAlignConfig());
    size_t edit_distance = result.editDistance;

    // convert edit distance into percent identity
    double perc_id = 1 - ((double)edit_distance / (double)ORF2_aa.size());
    return perc_id;
}