#include "gene_overlap.h"
#include "ORF_clustering.h"

inline std::pair<std::vector<int>, std::vector<indexPair>> combine_nodes(const std::vector<int>& ORF_nodes,
                                                                          const std::vector<indexPair>& ORF_coords,
                                                                          const std::vector<int>& TIS_nodes,
                                                                          const std::vector<indexPair>& TIS_coords)
{
    // initialise vectors to return
    std::vector<int> combined_nodes;
    std::vector<indexPair> combined_coords;

    // iterate over TIS_nodes, checking if they are the first entry in ORF_nodes
    int i = 0;
    for (; i < TIS_nodes.size(); i++)
    {
        if (TIS_nodes.at(i) != ORF_nodes.at(0))
        {
            combined_nodes.push_back(TIS_nodes.at(i));
            combined_coords.push_back(TIS_coords.at(i));
        } else
        {
            break;
        }
    }

    // as now aligned across TIS and ORF nodes, now combine
    int j = 0;
    for (; i < TIS_nodes.size(); i++)
    {
        combined_nodes.push_back(TIS_nodes.at(i));
        combined_coords.push_back(TIS_coords.at(i));
        combined_coords.back().second = ORF_coords.at(j).second;
        j++;
    }

    // add the remaining entries for ORF_nodes
    combined_nodes.insert(combined_nodes.end(), ORF_nodes.begin() + j, ORF_nodes.end());
    combined_coords.insert(combined_coords.end(), ORF_coords.begin() + j, ORF_coords.end());

    return {combined_nodes, combined_coords};
}

void reverse_ORFNodeVector(const ColoredCDBG<MyUnitigMap>& ccdbg,
                           const std::vector<Kmer>& head_kmer_arr,
                           std::pair<std::vector<int>, std::vector<indexPair>>& ORF2_nodes,
                           int& ORF2_start_node,
                           int& ORF2_end_node,
                           std::pair<int, size_t>& ORF2_5p,
                           std::pair<int, size_t>& ORF2_3p)
{
    // Reverse ORF2 node vector
    std::reverse(ORF2_nodes.first.begin(), ORF2_nodes.first.end());
    // Reverse ORF2 node coordinate vector
    std::reverse(ORF2_nodes.second.begin(), ORF2_nodes.second.end());

    // reverse sign on each ID in ORF2_nodes
    for (auto & node_id : ORF2_nodes.first)
    {
        node_id = node_id * -1;
    }

    // iterate over ORF2_nodes coordinate vector, reversing the coordinates relative to the end index of the node
    for (int i = 0; i < ORF2_nodes.second.size(); i++)
    {
        // get a reference to the unitig map object
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, ORF2_nodes.first.at(i));
        auto& um = um_pair.first;

        // get absolute last node index (same as unitig length minus 1 as zero indexed)
        size_t node_end = um.size - 1;
        // get difference from original start to absolute last node index
        size_t reversed_end = node_end - std::get<0>(ORF2_nodes.second.at(i));
        // get difference from original end to absolute last node index
        size_t reversed_start = node_end - std::get<1>(ORF2_nodes.second.at(i));
        // reassigned the entry in-place in ORF2_nodes.second
        ORF2_nodes.second[i] = std::make_pair(reversed_start, reversed_end);
    }

    // correct ORF2 5p and 3p positions, reversing the coordinates
    ORF2_5p.second = std::get<1>(ORF2_nodes.second.back());
    ORF2_3p.second = std::get<0>(ORF2_nodes.second[0]);

    // correct node labels
    ORF2_start_node = ORF2_nodes.first[0];
    ORF2_end_node = ORF2_nodes.first.back();
    ORF2_5p.first *= -1;
    ORF2_3p.first *= -1;
}

std::tuple<bool, std::vector<size_t>, std::vector<size_t>> slice_ORFNodeVector(const std::pair<std::vector<int>, std::vector<indexPair>>& ORF1_nodes,
                                                                               const std::pair<std::vector<int>, std::vector<indexPair>>& ORF2_nodes,
                                                                               const int& ORF2_start_node,
                                                                               const int& ORF2_end_node,
                                                                               const bool is_ref,
                                                                               const fm_index_coll& fm_idx)
{
    bool overlap_complete = false;

    // initialise vectors to capture overlapping node indices for ORFs
    std::vector<size_t> ORF_1_overlap_node_index;
    std::vector<size_t> ORF_2_overlap_node_index;

    // initialise start and end index lists for ORF2 within ORF1
    std::vector<size_t> start_index_list;
    std::vector<size_t> end_index_list;

    // look for presence of ORF2 start node in ORF1 nodes
    auto iter = std::find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_start_node);
    while(iter != std::get<0>(ORF1_nodes).end())
    {
        size_t index = iter - std::get<0>(ORF1_nodes).begin();
        start_index_list.push_back(std::move(index));
        iter = std::find(iter + 1, std::get<0>(ORF1_nodes).end(), ORF2_start_node);
    }

    // look for presence of ORF2 end node in ORF1 nodes
    iter = std::find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_end_node);
    while (iter != std::get<0>(ORF1_nodes).end())
    {
        size_t index = iter - std::get<0>(ORF1_nodes).begin();
        end_index_list.push_back(std::move(index));
        iter = std::find(iter + 1, std::get<0>(ORF1_nodes).end(), ORF2_end_node);
    }

    // if both start and end indexes are found, check if ORF2 sits within ORF1
    if (!start_index_list.empty() && !end_index_list.empty())
    {
        for (const auto& start_index : start_index_list)
        {
            for (const auto& end_index : end_index_list)
            {
                // check moving in forward direction
                std::vector<int> ORF1_nodes_sliced;
                if (start_index <= end_index)
                {
                    // slice ORF1 node vector from first entry to last
                    ORF1_nodes_sliced = std::vector<int> (std::get<0>(ORF1_nodes).begin() + start_index, std::get<0>(ORF1_nodes).begin() + end_index + 1);
                }

                // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
                if (ORF1_nodes_sliced == ORF2_nodes.first)
                {
                    for (size_t i1 = start_index; i1 < end_index + 1; i1++)
                    {
                        ORF_1_overlap_node_index.push_back(i1);
                    }
                    for (size_t i2 = 0; i2 < ORF2_nodes.first.size(); i2++)
                    {
                        ORF_2_overlap_node_index.push_back(i2);
                    }

                    overlap_complete = true;

                    break;
                }
            }
            // if chain is complete, don't iterate over any more start codons
            if (overlap_complete)
            {
                break;
            }
        }
    }

    // if overlap not found, check if beginning of ORF2 overlaps with end of ORF1. May be case that start and end nodes are found, but not complete match across ORF1
    if (!overlap_complete && !start_index_list.empty())
    {
        for (const auto& start_index : start_index_list)
        {
            // check moving in forward direction
            std::vector<int> ORF1_nodes_sliced;
            std::vector<int> ORF2_nodes_sliced;

            // slice ORF1 node vector from first entry to last
            ORF1_nodes_sliced = std::vector<int> (std::get<0>(ORF1_nodes).begin() + start_index, std::get<0>(ORF1_nodes).end());
            // check if ORF1_slice is too large to slice ORF1 (meaning likely ORF2 is reversed)
            if (ORF1_nodes_sliced.size() <= ORF2_nodes.first.size())
            {
                // slice ORF1 node vector from first entry to the equivalent length of ORF2_nodes_sliced
                ORF2_nodes_sliced = std::vector<int>(ORF2_nodes.first.begin(), ORF2_nodes.first.begin() + ORF1_nodes_sliced.size());
            }

            // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
            if (ORF1_nodes_sliced == ORF2_nodes_sliced)
            {
                for (size_t i1 = start_index; i1 < std::get<0>(ORF1_nodes).size(); i1++)
                {
                    ORF_1_overlap_node_index.push_back(i1);
                }
                for (size_t i2 = 0; i2 < ORF2_nodes_sliced.size(); i2++)
                {
                    ORF_2_overlap_node_index.push_back(i2);
                }

                // check if overlap is correct
                if (is_ref)
                {
                    std::vector<int> search_vector = ORF1_nodes.first;
                    search_vector.insert(search_vector.end(), ORF2_nodes.first.begin() + ORF_2_overlap_node_index.back() + 1, ORF2_nodes.first.end());

                    const auto present = path_search(search_vector, fm_idx);
                    if (!present.first)
                    {
                        continue;
                    }
                }

                overlap_complete = true;

                break;
            }
        }
    }
    // if overlap not found, check if end of ORF2 overlaps with beginning of ORF1.
    if (!overlap_complete && !end_index_list.empty())
    {
        for (const auto& end_index : end_index_list)
        {
            // check moving in forward direction
            std::vector<int> ORF1_nodes_sliced;
            std::vector<int> ORF2_nodes_sliced;

            // slice ORF1 node vector from first entry to last
            ORF1_nodes_sliced = std::vector<int> (std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).begin() + end_index + 1);
            // check that slice is small enough to correctly slice ORF2_nodes
            if (ORF1_nodes_sliced.size() <= ORF2_nodes.first.size())
            {
                // slice ORF2 node vector from last entry - end_index to the equivalent length of ORF2_nodes_sliced
                ORF2_nodes_sliced = std::vector<int>(ORF2_nodes.first.end() - end_index - 1, ORF2_nodes.first.end());
            }

            // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
            if (ORF1_nodes_sliced == ORF2_nodes_sliced)
            {
                for (size_t i1 = 0; i1 < ORF1_nodes_sliced.size(); i1++)
                {
                    ORF_1_overlap_node_index.push_back(i1);
                }
                for (size_t i2 = ORF2_nodes.first.size() - end_index - 1; i2 < ORF2_nodes.first.size(); i2++)
                {
                    ORF_2_overlap_node_index.push_back(i2);
                }
                // check if overlap is correct
                if (is_ref)
                {
                    std::vector<int> search_vector = ORF2_nodes.first;
                    search_vector.insert(search_vector.end(), ORF1_nodes.first.begin() + ORF_1_overlap_node_index.back() + 1, ORF1_nodes.first.end());

                    const auto present = path_search(search_vector, fm_idx);
                    if (!present.first)
                    {
                        continue;
                    }
                }

                overlap_complete = true;

                break;
            }
        }
    }
    return {overlap_complete, ORF_1_overlap_node_index, ORF_2_overlap_node_index};
}

ORFOverlapMap calculate_overlaps(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                 const std::vector<Kmer>& head_kmer_arr,
                                 const ORFVector& ORF_vector,
                                 const int DBG_overlap,
                                 const size_t max_overlap,
                                 const bool is_ref,
                                 const fm_index_coll& fm_idx)
{
    // initialise overlap map for each ORF per colour (each first ORF is the first ORF on positive strand etc.)
    ORFOverlapMap ORF_overlap_map;

    // intialise Eigen Triplet
    std::vector<ET> tripletList;

    // iterate over each ORF sequence with specific colours combination
    for (size_t ORF_index = 0; ORF_index < ORF_vector.size(); ORF_index++)
    {
        // iterate over nodes traversed by ORF and its TIS
        const auto& ORF_info = ORF_vector.at(ORF_index);
        std::unordered_set<int> ORF_nodes;
        std::copy(std::get<3>(ORF_info).begin(), std::get<3>(ORF_info).end(), std::inserter(ORF_nodes, ORF_nodes.end()));
        std::copy(std::get<0>(ORF_info).begin(), std::get<0>(ORF_info).end(), std::inserter(ORF_nodes, ORF_nodes.end()));
        for (const auto& node_traversed : ORF_nodes)
        {
            // add to triplet list, with temp_ORF_ID (row), node id (column) and set value as 1 (true)
            // convert node_traversed to size_t, minus 1 as unitigs are one-based, needs to be zero based
            const size_t abs_node_id = abs(node_traversed) - 1;
            tripletList.push_back(ET(ORF_index, abs_node_id, 1));
        }
    }

    // initialise sparse matrix
    Eigen::SparseMatrix<double> mat(ORF_vector.size(), head_kmer_arr.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    // conduct transposition + matrix multiplication to calculate ORFs sharing nodes
    auto ORF_overlap_mat = ((mat * mat.transpose()).pruned()).eval();

    // iterate over non-zero entries in matrix and calculate overlaps
    for (int outit = 0; outit < ORF_overlap_mat.outerSize(); ++outit)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator init(ORF_overlap_mat, outit); init; ++init)
        {
            // iterate through the bottom half of the symmetric matrix and ignore line of symmetry, pass if row <= col
            if(init.row() <= init.col())
            {
                continue;
            }

            // Assign temporary values for ORF1 and ORF2, not sorted by traversed node vector.
            auto temp_ORF1_ID = init.col();
            auto temp_ORF2_ID = init.row();

            // Get nodes traversed by genes. Order by length of the node vector; ORF1 is the longer of the two vectors. Need to copy ORF2, as may be reversed
            // check if ORF1 traverse more nodes than ORF2
            const bool temp_ORF1_longer = ((std::get<0>(ORF_vector.at(temp_ORF1_ID)).size() >= std::get<0>(ORF_vector.at(temp_ORF2_ID)).size()) ? true : false);

            // get respective ORF IDs
            const size_t& ORF1_ID = (temp_ORF1_longer ? temp_ORF1_ID : temp_ORF2_ID);
            const size_t& ORF2_ID = (temp_ORF1_longer ? temp_ORF2_ID : temp_ORF1_ID);

            // get reference to ORF1_info information and unpack
            const auto& ORF1_info = ORF_vector.at(ORF1_ID);
            const auto& ORF1_node_ids = std::get<0>(ORF1_info);
            const auto& ORF1_node_coords = std::get<1>(ORF1_info);
            const auto& ORF1_TIS_ids = std::get<3>(ORF1_info);
            const auto& ORF1_TIS_coords = std::get<4>(ORF1_info);
            const auto& ORF1_strand = std::get<5>(ORF1_info);
            const auto& ORF1_len = std::get<2>(ORF1_info);
            const bool ORF1_TIS = (!ORF1_TIS_ids.empty()) ? true : false;

            // unpack ORF2_info
            const auto& ORF2_info = ORF_vector.at(ORF2_ID);
            const auto& ORF2_node_ids = std::get<0>(ORF2_info);
            const auto& ORF2_node_coords = std::get<1>(ORF2_info);
            const auto& ORF2_TIS_ids = std::get<3>(ORF2_info);
            const auto& ORF2_TIS_coords = std::get<4>(ORF2_info);
            const auto& ORF2_strand = std::get<5>(ORF2_info);
            const auto& ORF2_len = std::get<2>(ORF2_info);
            const bool ORF2_TIS = (!ORF2_TIS_ids.empty()) ? true : false;

            // merge ORF and TIS nodes for ORF1 and ORF2
            const auto ORF1_nodes = std::move(combine_nodes(ORF1_node_ids, ORF1_node_coords, ORF1_TIS_ids, ORF1_TIS_coords));
            auto ORF2_nodes = std::move(combine_nodes(ORF2_node_ids, ORF2_node_coords, ORF2_TIS_ids, ORF2_TIS_coords));

            // initialise overlap type
            // n = no overlap
            // u = unidirectional overlap (3' of first overlaps with 5' of second ->->)
            // c = convergent overlap (3' of first overlaps with 3' of second -><-)
            // d = divergent overlap (5' of first overlaps with 5' of second <-->)
            // w = ORF lies completely within another
            // i = incompatible overlap (i.e. if overlap is greater than maximum or same stop codon shared in same frame)
            char overlap_type = 'n';

            // set reversed variable to determine if ORFs are in same strand,
            // or bidirectional if both orientations of ORF2 must be checked
            bool reversed = false;
            bool bidirectional = false;

            // set first ORF in relative ordering on strand
            int first_ORF = 1;

            // initialise absolute overlap
            size_t abs_overlap = 0;

            // initialise overlap_complete check
            bool overlap_complete = false;

            // work out strand of 3p node
            const bool ORF1_3p_strand = (std::get<0>(ORF1_nodes).back() >= 0) ? true : false;
            const bool ORF2_3p_strand = (ORF2_nodes.first.back() >= 0) ? true : false;

            //get ORF1 and ORF2 5' and 3' ends in pair <node_head_kmer, position>. ORF1 can be references, ORF2 cannot as can change
            const std::pair<int, size_t> ORF1_5p(std::get<0>(ORF1_nodes)[0], std::get<0>(std::get<1>(ORF1_nodes)[0]));
            const std::pair<int, size_t> ORF1_3p(std::get<0>(ORF1_nodes).back(), std::get<1>(std::get<1>(ORF1_nodes).back()));
            std::pair<int, size_t> ORF2_5p(ORF2_nodes.first[0], std::get<0>(ORF2_nodes.second[0]));
            std::pair<int, size_t> ORF2_3p(ORF2_nodes.first.back(), std::get<1>(ORF2_nodes.second.back()));

            // initialise values for start and end nodes. Again ORF1 can be references, ORF2 cannot as can change
            const int& ORF1_start_node = std::get<0>(ORF1_nodes)[0];
            const int& ORF1_end_node = std::get<0>(ORF1_nodes).back();
            int ORF2_start_node = ORF2_nodes.first[0];
            int ORF2_end_node = ORF2_nodes.first.back();

            // work out if node 1 is negative by checking strand
            bool negative = ORF1_strand;

            // determine if ORF2 should be reversed or not by looking at orientation of component nodes
            auto start_iter = std::find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_start_node);
            auto end_iter = std::find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_end_node);
            auto start_iter_rev = std::find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_start_node * -1);
            auto end_iter_rev = std::find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_end_node * -1);

            // determine if ORFs are reversed, or if ORFs are bidirectional,
            // meaning both orientations of ORF2 must be investigated
            if ((start_iter != std::get<0>(ORF1_nodes).end() && start_iter_rev != std::get<0>(ORF1_nodes).end()) ||
                (end_iter != std::get<0>(ORF1_nodes).end() && end_iter_rev != std::get<0>(ORF1_nodes).end()))
            {
                // if forward and reverse nodes present, then need to check both orientations
                bidirectional = true;
            }
                // else if reversed nodes present but forward are not, then reverse
            else if (start_iter == std::get<0>(ORF1_nodes).end() && end_iter == std::get<0>(ORF1_nodes).end())
            {
                reversed = true;
            }

            // check if 3p matches between the two ORFs and node is in the same strand without reversal.
            // If so, set as incompatible, and the overlap as the shorter of the two ORFS
            if (ORF1_3p == ORF2_3p && ORF1_3p_strand == ORF2_3p_strand && !reversed)
            {
                overlap_type = 'i';
                overlap_complete = true;
                abs_overlap = (ORF1_len >= ORF2_len) ? ORF2_len : ORF1_len;

                // set the first ORF to be the longest, regardless of strand
                first_ORF = (ORF1_len >= ORF2_len) ? 1 : 2;
            }

            // initialise overlap node_indices
            std::vector<size_t> ORF_1_overlap_node_index;
            std::vector<size_t> ORF_2_overlap_node_index;

            if (!overlap_complete)
            {
                // if not bidirectional, run as normal
                if (!bidirectional)
                {
                    if (reversed)
                    {
                        reverse_ORFNodeVector(ccdbg, head_kmer_arr, ORF2_nodes, ORF2_start_node, ORF2_end_node, ORF2_5p, ORF2_3p);
                    }
                    auto return_tuple = std::move(slice_ORFNodeVector(ORF1_nodes, ORF2_nodes, ORF2_start_node, ORF2_end_node, is_ref, fm_idx));
                    overlap_complete = std::get<0>(return_tuple);
                    ORF_1_overlap_node_index = std::get<1>(return_tuple);
                    ORF_2_overlap_node_index = std::get<2>(return_tuple);
                }
                // else, need to try in both directions, testing if overlap_complete
                else
                {
                    auto return_tuple = std::move(slice_ORFNodeVector(ORF1_nodes, ORF2_nodes, ORF2_start_node, ORF2_end_node, is_ref, fm_idx));
                    overlap_complete = std::get<0>(return_tuple);
                    ORF_1_overlap_node_index = std::get<1>(return_tuple);
                    ORF_2_overlap_node_index = std::get<2>(return_tuple);

                    if (!overlap_complete)
                    {
                        reverse_ORFNodeVector(ccdbg, head_kmer_arr, ORF2_nodes, ORF2_start_node, ORF2_end_node, ORF2_5p, ORF2_3p);
                        return_tuple = std::move(slice_ORFNodeVector(ORF1_nodes, ORF2_nodes, ORF2_start_node, ORF2_end_node, is_ref, fm_idx));
                        overlap_complete = std::get<0>(return_tuple);
                        ORF_1_overlap_node_index = std::get<1>(return_tuple);
                        ORF_2_overlap_node_index = std::get<2>(return_tuple);

                        // set reversed as true
                        reversed = true;
                    }
                }
            }

            // if overlap complete, calculate overlap if not already calculated
            // check that an overlapping region has been found. If so, calculate absolute overlap in base-pairs
            if (overlap_complete)
            {
                if (abs_overlap == 0)
                {
                    // initialise options for computing which type of overlap has occurred
                    size_t overlap_start = 0;
                    size_t overlap_end = 0;

                    for (size_t i = 0; i < ORF_1_overlap_node_index.size(); i++)
                    {
                        // get the index of the corresponding overlapping nodes between ORF1 and ORF2
                        const size_t& ORF1_overlap_node = ORF_1_overlap_node_index[i];
                        const size_t& ORF2_overlap_node = ORF_2_overlap_node_index[i];

                        // get the node coordinates traversed by each ORF
                        const size_t& ORF1_start = std::get<0>(std::get<1>(ORF1_nodes)[ORF1_overlap_node]);
                        const size_t& ORF1_end = std::get<1>(std::get<1>(ORF1_nodes)[ORF1_overlap_node]);
                        const size_t& ORF2_start = std::get<0>(ORF2_nodes.second[ORF2_overlap_node]);
                        const size_t& ORF2_end = std::get<1>(ORF2_nodes.second[ORF2_overlap_node]);

                        // get the first node involved in the overlap
                        const int overlap_node = std::get<0>(ORF1_nodes)[ORF1_overlap_node];

                        // check there is an intersection
                        if ((ORF1_start <= ORF2_end) && (ORF2_start <= ORF1_end))
                        {
                            size_t low_index = std::max(ORF1_start, ORF2_start);
                            size_t high_index = std::min(ORF1_end, ORF2_end);


                            // calculate overlap for first node, and work out overlap orientation
                            if (i == 0)
                            {
                                // calculate overlap. Add 1 as positions are zero-indexed, do not negate node_overlap
                                abs_overlap += (high_index - low_index) + 1;

                                // check what are the beginning and end indexes of the overlap, if they match a certain prime end, then can calculate the type of overlap
                                if (!reversed)
                                {
                                    if (overlap_node == ORF1_5p.first && low_index == ORF1_5p.second)
                                    {
                                        overlap_start = 15;
                                    } else if (overlap_node == ORF2_5p.first && low_index == ORF2_5p.second)
                                    {
                                        overlap_start = 25;
                                    }
                                } else
                                {
                                    if (overlap_node == ORF1_5p.first && low_index == ORF1_5p.second)
                                    {
                                        overlap_start = 15;
                                    } else if (overlap_node == ORF2_5p.first && low_index == ORF2_5p.second)
                                    {
                                        overlap_start = 25;
                                    } else if (overlap_node == ORF1_3p.first && low_index == ORF1_3p.second)
                                    {
                                        overlap_start = 13;
                                    } else if (overlap_node == ORF2_3p.first && low_index == ORF2_3p.second)
                                    {
                                        overlap_start = 23;
                                    }
                                }
                            } else
                            {
                                // calculate node overlap from DBG structure.
                                size_t node_overlap = (DBG_overlap - low_index);

                                //if index
                                if (high_index >= DBG_overlap)
                                {
                                    // calculate overlap, negating node overlap. Add one as zero indexed.
                                    abs_overlap += ((high_index - low_index) - node_overlap) + 1;
                                }
                            }

                            // if the node is the last in the overlap, check what end matches the overlap
                            if (i == ORF_1_overlap_node_index.size() - 1)
                            {
                                if (!reversed)
                                {
                                    if (overlap_node == ORF1_3p.first && high_index == ORF1_3p.second)
                                    {
                                        overlap_end = 13;
                                    } else if (overlap_node == ORF2_3p.first && high_index == ORF2_3p.second)
                                    {
                                        overlap_end = 23;
                                    }
                                } else
                                {
                                    if (overlap_node == ORF1_5p.first && high_index == ORF1_5p.second)
                                    {
                                        overlap_end = 15;
                                    } else if (overlap_node == ORF2_5p.first && high_index == ORF2_5p.second)
                                    {
                                        overlap_end = 25;
                                    } else if (overlap_node == ORF1_3p.first && high_index == ORF1_3p.second)
                                    {
                                        overlap_end = 13;
                                    } else if (overlap_node == ORF2_3p.first && high_index == ORF2_3p.second)
                                    {
                                        overlap_end = 23;
                                    }
                                }
                            }
                        }
                            // if no intersection detection, determine how ORFs are ordered
                        else if ((ORF1_start > ORF2_end && !negative) || (ORF2_start > ORF1_end && negative))
                        {
                            first_ORF = 2;
                        }
                    }

                    // check if TIS present for specific overlaps and negate
                    if ((overlap_start == 25 && ORF2_TIS) || (overlap_start == 15 && ORF1_TIS))
                    {
                        if (abs_overlap < 16)
                        {
                            abs_overlap = 0;
                        } else
                        {
                            abs_overlap -= 16;
                        }
                    }

                    // convert overlap_start and overlap_end to string to enable creation of ID for switch
                    // check if any overlap detected.
                    std::string overlap_ID_str = std::to_string(overlap_start) + std::to_string(overlap_end);
                    int overlap_ID = std::stoi(overlap_ID_str);

                    // go over combinations of overlap start and overlap end to determine type of overlap
                    switch (overlap_ID)
                    {
                        // unidirectional
                        case 1523:
                            overlap_type = 'u';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            if (!negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                        case 2513:
                            overlap_type = 'u';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            if (negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                            // ORF lies completely within another
                        case 1513:
                            overlap_type = 'w';
                            // ORF1 sits fully in ORF2
                            first_ORF = 2;
                            // think about case where ORF1 and ORF2 are reverse complements of eachother
                            if (ORF1_5p == ORF2_3p && ORF1_3p == ORF2_5p && !negative)
                            {
                                first_ORF = 1;
                            }
                            break;
                        case 2523:
                            overlap_type = 'w';
                            // ORF2 sits fully in ORF1, leave first_ORF = 1
                            break;
                        case 1315:
                            overlap_type = 'w';
                            // ORF1 sits fully in ORF2
                            first_ORF = 2;
                            break;
                        case 2325:
                            overlap_type = 'w';
                            // ORF2 sits fully in ORF1, leave first_ORF = 1
                            break;
                            // convergent
                        case 2313:
                            overlap_type = 'c';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            if (negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                        case 1323:
                            overlap_type = 'c';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            if (!negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                        case 1313:
                            overlap_type = 'c';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            if (negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                            // divergent
                        case 2515:
                            overlap_type = 'd';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            if (negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                        case 1525:
                            overlap_type = 'd';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            // check for case where 5p of ORF1 and 3p of ORF2 match, if they do keep ORF as first ORF
                            if (ORF1_5p != ORF2_3p && !negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                        case 1515:
                            overlap_type = 'd';
                            // adjust for negativity of ORF1 from default (first_ORF = 1)
                            if (!negative)
                            {
                                first_ORF = 2;
                            }
                            break;
                    }
                }
                // if overlap is greater than max_overlap, set as incompatible
                if (abs_overlap > max_overlap)
                {
                    overlap_type = 'i';
                } else if (abs_overlap == 0)
                {
                    overlap_type = 'n';
                }

                // add overlap type to map, where the first ORF on the positive strand is the second key,
                // and the second ORF is the first key (edge weights are generated for the sink node based on Balrog DAG)

                // overlap_tuple contains the overlap_type, abs_overlap, and the strand of the first and second ORF in the order they appear in the map
                if (first_ORF == 1)
                {
                    std::pair<char, size_t> overlap_tuple(overlap_type, abs_overlap);
                    ORF_overlap_map[ORF2_ID][ORF1_ID] = std::move(overlap_tuple);

                } else {
                    std::pair<char, size_t> overlap_tuple(overlap_type, abs_overlap);
                    ORF_overlap_map[ORF1_ID][ORF2_ID] = std::move(overlap_tuple);
                }
            }
        }
    }

    return ORF_overlap_map;
}
