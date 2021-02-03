#include "ggCaller_classes.h"
// return robin_hood::unordered_map<std::string, std::string> in future, which pairs each ORF to its overlapping ORFs
std::pair<ORFOverlapMap, FullORFMap> calculate_overlaps(const unitigMap& unitig_map,
                                                                     const StrandORFNodeMap& ORF_node_paths,
                                                                     const std::pair<ORFColoursMap, std::vector<std::string>>& ORF_colours_tuple,
                                                                     const int& DBG_overlap,
                                                                     const size_t& max_overlap)
{
    // initialise overlap group map for each colour
    //robin_hood::unordered_map<std::string, std::vector<std::unordered_set<std::string>>> overlap_group_map;

    // initialse an ORF set to be returned with overlap Map for easy iteration of gene scoring
    //std::unordered_set<std::string> ORF_set;

    // initialise full ORF map for easy iteration of gene scoring and gene removal
    FullORFMap full_ORF_map;

    // initialise overlap map for each ORF per colour (each first ORF is the first ORF on positive strand etc.)
    ORFOverlapMap ORF_overlap_map;

    // iterate over each colour combination in ORF_colours_map
    for (auto colit = ORF_colours_tuple.second.begin(); colit < ORF_colours_tuple.second.end(); colit++)
    {
        // initialise overlap group vector
        std::vector<std::unordered_set<std::string>> overlap_group_vector;

        // intialise Eigen Triplet
        std::vector<ET> tripletList;

        // initialise ORF ID count;
        size_t ORF_ID = 0;

        // map to determine the sequence and strand of a gene
        robin_hood::unordered_map<size_t, std::pair<std::string, std::string>> ORF_ID_map;

        // iterate over each ORF sequence with specific colours combination
        for (const auto& ORF_seq : ORF_colours_tuple.first.at(*colit))
        {
            // search for the sequence within the positive ORF paths
            if (ORF_node_paths.at("+").find(ORF_seq) != ORF_node_paths.at("+").end())
            {
                // add ORF sequence to full_ORF_map
                full_ORF_map[*colit][ORF_seq] = "+";

                std::pair seq_pair(ORF_seq, "+");
                ORF_ID_map[ORF_ID] = std::move(seq_pair);
                // iterate over nodes traversed by ORF
                for (const auto& node_traversed : (ORF_node_paths.at("+").at(ORF_seq).first))
                {
                    // add to triplet list, with ORRF_ID (row), node id (column) and set value as 1 (true)
                    tripletList.push_back(ET(ORF_ID, unitig_map.at(node_traversed).unitig_id, 1));
                }
                ORF_ID++;
            }

            // repeat for reverse node paths
            else if (ORF_node_paths.at("-").find(ORF_seq) != ORF_node_paths.at("-").end())
            {
                // add ORF sequence to full_ORF_map
                full_ORF_map[*colit][ORF_seq] = "-";

                std::pair seq_pair(ORF_seq, "-");
                ORF_ID_map[ORF_ID] = std::move(seq_pair);
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
                auto temp_ORF1 = ORF_ID_map.at(init.col());
                auto temp_ORF2 = ORF_ID_map.at(init.row());

                // Get nodes traversed by genes. Order by length of the node vector; ORF1 is the longer of the two vectors
                auto ORF1_nodes = ((ORF_node_paths.at(temp_ORF1.second).at(temp_ORF1.first).first.size() >= ORF_node_paths.at(temp_ORF2.second).at(temp_ORF2.first).first.size()) ? ORF_node_paths.at(temp_ORF1.second).at(temp_ORF1.first) : ORF_node_paths.at(temp_ORF2.second).at(temp_ORF2.first));
                auto ORF2_nodes = ((ORF_node_paths.at(temp_ORF1.second).at(temp_ORF1.first).first.size() >= ORF_node_paths.at(temp_ORF2.second).at(temp_ORF2.first).first.size()) ? ORF_node_paths.at(temp_ORF2.second).at(temp_ORF2.first) : ORF_node_paths.at(temp_ORF1.second).at(temp_ORF1.first));
                auto ORF1 = ((ORF_node_paths.at(temp_ORF1.second).at(temp_ORF1.first).first.size() >= ORF_node_paths.at(temp_ORF2.second).at(temp_ORF2.first).first.size()) ? temp_ORF1 : temp_ORF2);
                auto ORF2 = ((ORF_node_paths.at(temp_ORF1.second).at(temp_ORF1.first).first.size() >= ORF_node_paths.at(temp_ORF2.second).at(temp_ORF2.first).first.size()) ? temp_ORF2 : temp_ORF1);

                // for debugging
                auto row = init.row();
                auto col = init.col();

                if ((ORF1.first == "ATTGGTAATTTTGTATATGCCTACTTTTAGGGATGATGCTTCAACGAAGGCGTATAATCTTGATTACGATAAGGTAATAAATTCATTTCAAGATTTTTATAACAGAAAAGTTAAAGTATTGATTCGTTTTCATCCAAATGTAGATAATACATTTTTTAATAATACTGATAAAAGATTAATTAATGTGACAGATTATCCTAATCCGCAGGATTTAATGTTTGTGGCTGATATTATGATTTCAGACTATTCGTCAGCACCCATAGATTTTTTGTTATTAAATCGAGTAGTCTTTCTGTATCTACCAGATTTTAAAGAATATCAGAGCGATAAAAATCCGTTTTTTGAAGTTTTCAAAGTTTCGAAAACCAAAGGCATTGCGCTTGATCCGTTTGATGAGATTATTGGTCGCTTCCAGTTTGGCGTTAGAATAGTGTAG" || ORF2.first == "ATTGGTAATTTTGTATATGCCTACTTTTAGGGATGATGCTTCAACGAAGGCGTATAATCTTGATTACGATAAGGTAATAAATTCATTTCAAGATTTTTATAACAGAAAAGTTAAAGTATTGATTCGTTTTCATCCAAATGTAGATAATACATTTTTTAATAATACTGATAAAAGATTAATTAATGTGACAGATTATCCTAATCCGCAGGATTTAATGTTTGTGGCTGATATTATGATTTCAGACTATTCGTCAGCACCCATAGATTTTTTGTTATTAAATCGAGTAGTCTTTCTGTATCTACCAGATTTTAAAGAATATCAGAGCGATAAAAATCCGTTTTTTGAAGTTTTCAAAGTTTCGAAAACCAAAGGCATTGCGCTTGATCCGTTTGATGAGATTATTGGTCGCTTCCAGTTTGGCGTTAGAATAGTGTAG")
                && (ORF1.first == "TTTTGTCCTTTCTTTTTTGATGTTCAAAGCGATAAAAATCCGTTTTTTGAAGTTTTCAAAGTTTCGAAAACCAAAGGCATTGCGCTTGATAAGTTTGATGAGATTATTGGTCGCTTCCAGTTTGGCATTAGAATAG" || ORF2.first == "TTTTGTCCTTTCTTTTTTGATGTTCAAAGCGATAAAAATCCGTTTTTTGAAGTTTTCAAAGTTTCGAAAACCAAAGGCATTGCGCTTGATAAGTTTGATGAGATTATTGGTCGCTTCCAGTTTGGCATTAGAATAG"))
                {
                    int test = 1;
                }


                // get ORF1 and ORF2 5' and 3' ends in pair <node_head_kmer, position>
                std::pair<std::string, size_t> ORF1_5p(ORF1_nodes.first[0], std::get<0>(ORF1_nodes.second[0]));
                std::pair<std::string, size_t> ORF1_3p(ORF1_nodes.first.back(), std::get<1>(ORF1_nodes.second.back()));
                std::pair<std::string, size_t> ORF2_5p(ORF2_nodes.first[0], std::get<0>(ORF2_nodes.second[0]));
                std::pair<std::string, size_t> ORF2_3p(ORF2_nodes.first.back(), std::get<1>(ORF2_nodes.second.back()));

                // set reversed variable to determine if ORFs are in same strand
                bool reversed = false;

                // set negative variable to determine if both strands are negative
                bool negative = false;

                // initialise which ORF is first positionally in positive strand. Default is 1.
                int first_ORF = 1;

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

                    // set reversed to true
                    reversed = true;

                    // correct ORF2 5p and 3p positions, reversing the coordinates
                    ORF2_5p.second = std::get<1>(ORF2_nodes.second.back());
                    ORF2_3p.second = std::get<0>(ORF2_nodes.second[0]);
                }
                // check if ORF1 is negative, enabling determination of positions of ORFs on positive strand
                if (ORF1.second == "-")
                {
                    negative = true;
                }

                // initialise vectors to capture overlapping node indices for ORFs
                std::vector<size_t> ORF_1_overlap_node_index;
                std::vector<size_t> ORF_2_overlap_node_index;

                // generate a set for ORF2 nodes to enable fast querying of ORF1 node elements
                std::unordered_set<std::string> ORF2_node_set(ORF2_nodes.first.begin(), ORF2_nodes.first.end());

                // initialise absolute overlap
                size_t abs_overlap = 0;

                // initialise values for start and end nodes
                std::string ORF1_start_node = ORF1_nodes.first[0];
                std::string ORF1_end_node = ORF1_nodes.first.back();
                std::string ORF2_start_node = ORF2_nodes.first[0];
                std::string ORF2_end_node = ORF2_nodes.first.back();

                // set overlap_complete as false, enable correct overlap to be detected if nodes are repeated
                bool overlap_complete = false;

                // initialise overlap type
                // n = no overlap
                // u = unidirectional overlap (3' of first overlaps with 5' of second ->->)
                // c = convergent overlap (3' of first overlaps with 3' of second -><-)
                // d = divergent overlap (5' of first overlaps with 5' of second <-->)
                // w = ORF lies completely within another
                // i = incompatible overlap (i.e. if overlap is greater than maximum or same stop codon shared in same frame)
                char overlap_type = 'n';

                // check if genes are unidirectional and share same stop codon. If so, genes are incompatible and overlap will be equal to the shorter of the two ORFs
                if (!reversed && ORF1_3p == ORF2_3p)
                {
                    overlap_type = 'i';
                    overlap_complete = true;
                    abs_overlap = ((ORF1.first.size() >= ORF2.first.size()) ? ORF2.first.size() : ORF1.first.size()) - 16;
                    // set the first ORF to be the longest if not negative. If is negative, both ORFs start at same position on positive strand, so leave as default
                    if (!negative)
                    {
                        first_ORF = ((ORF1.first.size() >= ORF2.first.size()) ? 1 : 2);
                    }
                }

                // check overlap not already completed by matching stop codons
                if (!overlap_complete)
                {
                    // beginning of ORF1 overlaps with end of ORF2. If first node in ORF1 is in ORF2, ensure the nodes match up until the end of ORF2
                    if (ORF2_node_set.find(ORF1_start_node) != ORF2_node_set.end())
                    {
                        // find all entries in ORF2 that matches the first entry in ORF1
                        std::vector<std::size_t> index_list;

                        // define iterator for detection of matching nodes in find methods
                        auto iter = find(ORF2_nodes.first.begin(), ORF2_nodes.first.end(), ORF1_start_node);
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
                }
                // end of ORF1 overlaps with beginning of ORF2. If first node in ORF1 is in ORF2, ensure the nodes match up until the end of ORF2
                if (!overlap_complete)
                {
                    if (!overlap_complete && ORF2_node_set.find(ORF1_end_node) != ORF2_node_set.end())
                    {
                        // find all entries in ORF2 that matches the last entry in ORF1
                        std::vector<std::size_t> index_list;
                        auto iter = find(ORF2_nodes.first.begin(), ORF2_nodes.first.end(), ORF1_end_node);
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
                    auto iter = find(ORF1_nodes.first.begin(), ORF1_nodes.first.end(), ORF2_start_node);
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
                            // ensure indexes are correctly orientated (may be reversed if one gene is artifically generated)
                            if (start_index > end_index)
                            {
                                continue;
                            }
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
                                // the overlap between the two is equivalent to the length of ORF2, minus upstream region
                                abs_overlap = ORF2.first.size() - 16;

                                // gene pair is set overlap_types as 'w'
                                overlap_type = 'w';
                                // First ORF default is ORF1 as it overlaps ORF2 completely, leave as is.
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
                if (overlap_complete)
                {
                    if (abs_overlap == 0)
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

                        // get the first node involved in the overlap
                        std::string overlap_node = ORF1_nodes.first[ORF1_overlap_node];

                        // initialise options for computing which type of overlap has occurred
                        size_t overlap_start;
                        size_t overlap_end;

                        // check there is an intersection
                        if ((ORF1_start <= ORF2_end) && (ORF2_start <= ORF1_end))
                        {
                            size_t low_index = std::max(ORF1_start, ORF2_start);
                            size_t high_index = std::min(ORF1_end, ORF2_end);

                            // calculate overlap. Add 1 as positions are zero-indexed
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

                            // if overlap is only single node, calculate overlap_end also
                            if (ORF_1_overlap_node_index.size() == 1)
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

                            // get the first node involved in the overlap
                            overlap_node = ORF1_nodes.first[ORF1_overlap_node];

                            // check there is an intersection
                            if ((ORF1_start <= ORF2_end) && (ORF2_start <= ORF1_end))
                            {
                                size_t low_index = std::max(ORF1_start, ORF2_start);
                                size_t high_index = std::min(ORF1_end, ORF2_end);

                                // calculate node overlap from DBG structure.
                                size_t node_overlap = (DBG_overlap - low_index);

                                // calculate overlap, negating node overlap. Add one as zero indexed.
                                abs_overlap += ((high_index - low_index) - node_overlap) + 1;

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
                        }
                        // convert overlap_start and overlap_end to string to enable creation of ID for switch
                        // check if any overlap detected.
                        if (abs_overlap > 0)
                        {
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
                    }
                    // if overlap is greater than max_overlap, set as incompatible
                    if (abs_overlap > max_overlap)
                    {
                        overlap_type = 'i';
                    }

                    // add overlap type to map, where the first ORF on the positive strand is the second key,
                    // and the second ORF is the first key (edge weights are generated for the sink node based on Balrog DAG)

                    // overlap_tuple contains the overlap_type, abs_overlap, and the strand of the first and second ORF in the order they appear in the map
                    if (first_ORF == 1)
                    {
                        std::pair<char, size_t> overlap_tuple(overlap_type, abs_overlap);
                        ORF_overlap_map[*colit][ORF2.first][ORF1.first] = overlap_tuple;
                    } else {
                        std::pair<char, size_t> overlap_tuple(overlap_type, abs_overlap);
                        ORF_overlap_map[*colit][ORF1.first][ORF2.first] = overlap_tuple;
                    }
                }
            }
        }
    }
    auto ORF_overlap_pair = std::make_pair(ORF_overlap_map, full_ORF_map);
    return ORF_overlap_pair;
}
