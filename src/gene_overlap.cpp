#include "ggCaller_classes.h"

// return robin_hood::unordered_map<std::string, std::string> in future, which pairs each ORF to its overlapping ORFs
void calculate_overlaps(const unitigMap& unitig_map,
                        const StrandORFNodeMap& ORF_node_paths,
                        const std::unordered_map<std::string, std::vector<std::string>> ORF_colours_map,
                        const int& DBG_overlap,
                        const size_t& max_overlap)
{

    // iterate over each colour combination in ORF_colours_map
    for (auto it = ORF_colours_map.begin(); it != ORF_colours_map.end(); it++)
    {
        cout << "Colour: " << it->first << endl;
        // intialise Eigen Triplet
        std::vector<ET> tripletList;

        // initialise ORF ID count;
        size_t ORF_ID = 0;

        // map to determine the sequence and strand of a gene
        robin_hood::unordered_map<size_t, std::pair<std::string, std::string>> ORF_map;

        // iterate over each ORF sequence with specific colours combination
        for (const auto& ORF_seq : it->second)
        {
            // search for the sequence within the positive ORF paths
            if (ORF_node_paths.at("+").find(ORF_seq) != ORF_node_paths.at("+").end())
            {
                std::pair seq_pair(ORF_seq, "+");
                ORF_map[ORF_ID] = std::move(seq_pair);
                // iterate over nodes traversed by ORF
                for (const auto& node_traversed : (ORF_node_paths.at("+").at(ORF_seq).first))
                {
                    // add to triplet list, with ORRF_ID (row), node id (column) and set value as 1 (true)
                    tripletList.push_back(ET(ORF_ID, unitig_map.at(node_traversed).unitig_id, 1));
                }
                ORF_ID++;
            }

            // repeat for reverse node paths
            if (ORF_node_paths.at("-").find(ORF_seq) != ORF_node_paths.at("-").end())
            {
                std::pair seq_pair(ORF_seq, "-");
                ORF_map[ORF_ID] = std::move(seq_pair);
                // iterate over nodes traversed by ORF
                for (const auto& node_traversed : (ORF_node_paths.at("-").at(ORF_seq).first))
                {
                    // add to triplet list, with ORF_ID (row), node id (column) and set value as 1 (true)
                    tripletList.push_back(ET(ORF_ID, unitig_map.at(node_traversed).unitig_id, 1));
                }
                ORF_ID++;
            }
        }
        // initialise sparse matrix
        Eigen::SparseMatrix<double> mat(ORF_ID, unitig_map.size());
        mat.setFromTriplets(tripletList.begin(), tripletList.end());

        // conduct transposition + matrix multiplication to calculate ORFs sharing nodes
        auto ORF_overlap_mat = ((mat * mat.transpose()).pruned()).eval();

        // iterate over non-zero entries in matrix
        for (int k = 0; k < ORF_overlap_mat.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(ORF_overlap_mat, k); it; ++it)
            {
                // iterate through the bottom half of the symmetric matrix and ignore line of symmetry, pass if row <= col
                if(it.row() <= it.col())
                {
                    continue;
                }

                // gather information on overlapping ORFs. Assigned ORF1 as the longer ORF, ORF2 as the shorter ORF
                auto ORF1 = (ORF_map.at(it.col()).first.size() >= ORF_map.at(it.row()).first.size()) ? ORF_map.at(it.col()) : ORF_map.at(it.row());
                auto ORF2 = (ORF_map.at(it.col()).first.size() >= ORF_map.at(it.row()).first.size()) ? ORF_map.at(it.row()) : ORF_map.at(it.col());

                // gather nodes traversed by genes
                auto ORF1_nodes = ORF_node_paths.at(ORF1.second).at(ORF1.first);
                auto ORF2_nodes = ORF_node_paths.at(ORF2.second).at(ORF2.first);


                // check if both strands are the same. If not, reverse the nodes of ORF2 and their within-node coordinates
                if (ORF1.second != ORF2.second)
                {
                    // Reverse ORF2 node vector
                    std::reverse(ORF2_nodes.first.begin(), ORF2_nodes.first.end());
                    // Reverse ORF2 node coordinate vector
                    std::reverse(ORF2_nodes.second.begin(), ORF2_nodes.second.end());

                    // iterate over ORF2_nodes coordinate vector, reversing the coordinates relative to the end index of the node
                    for (auto & node_coord : ORF2_nodes.second)
                    {
                        // get absolute last node index
                        size_t node_end = std::get<2>(node_coord);
                        // get difference from original start to absolute last node index
                        size_t reversed_end = node_end - std::get<0>(node_coord);
                        // get difference from original end to absolute last node index
                        size_t reversed_start = node_end - std::get<1>(node_coord);
                        // reassigned the entry in-place in ORF2_nodes.second
                        node_coord = make_tuple(reversed_start, reversed_end, node_end);
                    }
                }

                // initialise vectors to capture overlapping node indices for ORFs
                std::vector<size_t> ORF_1_overlap_node_index;
                std::vector<size_t> ORF_2_overlap_node_index;

                // initialise absolute overlap
                size_t abs_overlap = 0;

                // generate a set for ORF2 nodes to enable fast querying of ORF1 node elements
                std::unordered_set<std::string> ORF2_node_set(ORF2_nodes.first.begin(), ORF2_nodes.first.end());

                std::string ORF1_start_node = ORF1_nodes.first[0];
                std::string ORF1_end_node = ORF1_nodes.first.back();
                std::string ORF2_start_node = ORF2_nodes.first[0];
                std::string ORF2_end_node = ORF2_nodes.first.back();

                // set overlap_complete as false, enable correct overlap to be detected if nodes are repeated
                bool overlap_complete = false;

                // define iterator for detection of matching nodes in find methods
                std::vector<std::string>::iterator iter;

                // beginning of ORF1 overlaps with end of ORF2. If first node in ORF1 is in ORF2, ensure the nodes match up until the end of ORF2
                if (ORF2_node_set.find(ORF1_start_node) != ORF2_node_set.end())
                {
                    // find all entries in ORF2 that matches the first entry in ORF1
                    std::vector<std::size_t> index_list;
                    iter = find(ORF2_nodes.first.begin(), ORF2_nodes.first.end(), ORF1_start_node);
                    while(iter != ORF2_nodes.first.end())
                    {
                        size_t index = iter - ORF2_nodes.first.begin();
                        index_list.push_back(index);
                        iter = find(iter + 1, ORF2_nodes.first.end(), ORF1_start_node);
                    }

                    // iterate over all indexes of the start node of ORF1
                    for (auto start_index : index_list)
                    {
                        // slice ORF2 node vector from the last entry to the end
                        std::vector<std::string> ORF2_nodes_sliced(ORF2_nodes.first.begin() + start_index, ORF2_nodes.first.end());
                        // slice ORF1 node vector from first entry to the equivalent length of ORF2_nodes_sliced
                        std::vector<std::string> ORF1_nodes_sliced(ORF1_nodes.first.begin(), ORF1_nodes.first.begin() + ORF2_nodes_sliced.size());
                        // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
                        if (ORF1_nodes_sliced == ORF2_nodes_sliced)
                        {
                            for (size_t i1 = 0; i1 < ORF1_nodes_sliced.size(); i1++)
                            {
                                ORF_1_overlap_node_index.push_back(i1);
                            }
                            for (size_t i2 = start_index; i2 < ORF2_nodes.first.size(); i2++)
                            {
                                ORF_2_overlap_node_index.push_back(i2);
                            }
                            overlap_complete = true;
                            break;
                        }
                    }
                }
                // end of ORF1 overlaps with beginning of ORF2. If first node in ORF1 is in ORF2, ensure the nodes match up until the end of ORF2
                if (!overlap_complete)
                {
                    if (!overlap_complete && ORF2_node_set.find(ORF1_end_node) != ORF2_node_set.end())
                    {
                        // find all entries in ORF2 that matches the last entry in ORF1
                        std::vector<std::size_t> index_list;
                        iter = find(ORF2_nodes.first.begin(), ORF2_nodes.first.end(), ORF1_end_node);
                        while(iter != ORF2_nodes.first.end())
                        {
                            size_t index = iter - ORF2_nodes.first.begin();
                            // add one to index to enable correct vector slicing from end of vector
                            index_list.push_back(index + 1);
                            iter = find(iter + 1, ORF2_nodes.first.end(), ORF1_end_node);
                        }

                        // iterate over all indexes of the start node of ORF1
                        for (auto end_index : index_list)
                        {
                            // slice ORF2 node vector from the beginning to the entry
                            std::vector<std::string> ORF2_nodes_sliced(ORF2_nodes.first.begin(), ORF2_nodes.first.begin() + end_index);
                            // slice ORF1 node vector from equivalent distance from end for ORF1 beginning to the start_index to the end
                            std::vector<std::string> ORF1_nodes_sliced(ORF1_nodes.first.end() - end_index, ORF1_nodes.first.end());
                            // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes and break
                            if (ORF1_nodes_sliced == ORF2_nodes_sliced)
                            {
                                for (size_t i1 = ORF1_nodes.first.size() - end_index; i1 < ORF1_nodes.first.size(); i1++)
                                {
                                    ORF_1_overlap_node_index.push_back(i1);
                                }
                                for (size_t i2 = 0; i2 < ORF2_nodes_sliced.size(); i2++)
                                {
                                    ORF_2_overlap_node_index.push_back(i2);
                                }
                                overlap_complete = true;
                                break;
                            }
                        }
                    }
                }
                if (!overlap_complete)
                {
                    // ORF2 is within ORF1. Look for first and last entries of ORF2 in ORF1 and slice
                    // find all entries of ORF2 start node in ORF1
                    std::vector<std::size_t> start_index_list;
                    iter = find(ORF1_nodes.first.begin(), ORF1_nodes.first.end(), ORF2_start_node);
                    while(iter != ORF1_nodes.first.end())
                    {
                        size_t index = iter - ORF1_nodes.first.begin();
                        start_index_list.push_back(index);
                        iter = find(iter + 1, ORF1_nodes.first.end(), ORF2_start_node);
                    }

                    // find all entries of ORF2 end node in ORF1
                    std::vector<std::size_t> end_index_list;
                    iter = find(ORF1_nodes.first.begin(), ORF1_nodes.first.end(), ORF2_end_node);
                    while(iter != ORF1_nodes.first.end())
                    {
                        size_t index = iter - ORF1_nodes.first.begin();
                        // add one to index to enable correct vector slicing from end of vector
                        end_index_list.push_back(index + 1);
                        iter = find(iter + 1, ORF1_nodes.first.end(), ORF2_end_node);
                    }

                    // iterate over all indexes of the start node and end node of ORF1. Ensure ORF2 is completely traversed (chain_complete)
                    bool chain_complete = false;
                    for (auto start_index : start_index_list)
                    {
                        for (auto end_index : end_index_list)
                        {
                            // slice ORF1 node vector from the beginning to the entry
                            std::vector<std::string> ORF1_nodes_sliced(ORF1_nodes.first.begin() + start_index, ORF1_nodes.first.begin() + end_index);
                            // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes and break
                            if (ORF1_nodes_sliced == ORF2_nodes.first)
                            {
                                for (size_t i1 = start_index; i1 < end_index; i1++)
                                {
                                    ORF_1_overlap_node_index.push_back(i1);
                                }
                                for (size_t i2 = 0; i2 < ORF2_nodes.first.size(); i2++)
                                {
                                    ORF_2_overlap_node_index.push_back(i2);
                                }
                                chain_complete = true;
                                overlap_complete = true;

                                // if all nodes of ORF2 match internally within ORF1,
                                // the overlap between the two is equivalent to the length of ORF2
                                abs_overlap = ORF2.first.size();
                                break;
                            }
                        }
                        // Ensure ORF2 is completely traversed (chain_complete)
                        if (chain_complete)
                        {
                            break;
                        }
                    }
                }

                // check that an overlapping region has been found. If so, calculate absolute overlap in base-pairs
                if (overlap_complete && abs_overlap == 0)
                {

                    // get the index of the corresponding overlapping nodes between ORF1 and ORF2
                    size_t ORF1_overlap_node = ORF_1_overlap_node_index[0];
                    size_t ORF2_overlap_node = ORF_2_overlap_node_index[0];


                    // calculate overlap of first node and get overlap of it's end
                    // get the node coordinates traversed by each ORF
                    size_t ORF1_start = std::get<0>(ORF1_nodes.second[ORF1_overlap_node]);
                    size_t ORF1_end = std::get<1>(ORF1_nodes.second[ORF1_overlap_node]);
                    size_t ORF2_start = std::get<0>(ORF2_nodes.second[ORF2_overlap_node]);
                    size_t ORF2_end = std::get<1>(ORF2_nodes.second[ORF2_overlap_node]);


                    // check there is an intersection
                    if ((ORF1_end >= ORF2_start) || (ORF2_end >= ORF1_start))
                    {
                        size_t low_index = std::max(ORF1_start, ORF2_start);
                        size_t high_index = std::min(ORF1_end, ORF2_end);

                        // calculate overlap. Add 1 as positions are zero-indexed
                        abs_overlap += (high_index - low_index) + 1;
                    }

                    // iterate over remaining matching indexes across the two ORFs. Negate node_overlap due to DBG structure.
                    for (size_t i = 1; i < ORF_1_overlap_node_index.size(); i++)
                    {
                        // get the index of the corresponding overlapping nodes between ORF1 and ORF2
                        ORF1_overlap_node = ORF_1_overlap_node_index[i];
                        ORF2_overlap_node = ORF_2_overlap_node_index[i];

                        // get the node coordinates traversed by each ORF
                        ORF1_start = std::get<0>(ORF1_nodes.second[ORF1_overlap_node]);
                        ORF1_end = std::get<1>(ORF1_nodes.second[ORF1_overlap_node]);
                        ORF2_start = std::get<0>(ORF2_nodes.second[ORF2_overlap_node]);
                        ORF2_end = std::get<1>(ORF2_nodes.second[ORF2_overlap_node]);

                        // check there is an intersection
                        if ((ORF1_end >= ORF2_start) || (ORF2_end >= ORF1_start))
                        {
                            size_t low_index = std::max(ORF1_start, ORF2_start);
                            size_t high_index = std::min(ORF1_end, ORF2_end);

                            // check that the ORF overlap high index is greater than the DBG_overlap,
                            // otherwise pass as this sequence has already been counted
                            if (high_index >= DBG_overlap)
                            {
                                // calculate node overlap from DBG structure.
                                size_t node_overlap = (DBG_overlap - low_index);

                                // calculate overlap, negating node overlap
                                abs_overlap += ((high_index - low_index) - node_overlap) + 1;
                            }
                        }
                    }
                }
            }
        }
        int j = 1;
    }
}
