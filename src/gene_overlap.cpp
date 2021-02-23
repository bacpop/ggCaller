#include "ggCaller_classes.h"

std::pair<ORFOverlapMap, FullORFMap> calculate_overlaps(const unitigMap& unitig_map,
                                                        const ORFNodeMap& ORF_node_paths,
                                                        const std::unordered_map<std::string, NodeStrandMap>& pos_strand_map,
                                                        const std::pair<ORFColoursMap, std::vector<std::string>>& ORF_colours_pair,
                                                        const int& DBG_overlap,
                                                        const size_t& max_overlap)
{
    // initialise full ORF map for easy iteration of gene scoring and gene removal
    FullORFMap full_ORF_map;

    // initialise overlap map for each ORF per colour (each first ORF is the first ORF on positive strand etc.)
    ORFOverlapMap ORF_overlap_map;


    // iterate over each colour combination in ORF_colours_map
    #pragma omp parallel
    {
        // initialise private thread items
        FullORFMap full_ORF_map_private;
        ORFOverlapMap ORF_overlap_map_private;

        // iterate over each colour
        #pragma omp for nowait
        for (auto colit = ORF_colours_pair.second.begin(); colit < ORF_colours_pair.second.end(); colit++)
        {
            // intialise Eigen Triplet
            std::vector<ET> tripletList;

            // initialise ORF ID count;
            size_t ORF_ID = 0;

            // map to determine ID of a gene
            robin_hood::unordered_map<size_t, std::string> ORF_ID_map;

            // iterate over each ORF sequence with specific colours combination
            for (const auto& ORF_seq : ORF_colours_pair.first.at(*colit))
            {
                // assign placeholder for ORF scores
                full_ORF_map_private[*colit][ORF_seq] = 'n';

                // assign ORF seq with unique id
                ORF_ID_map[ORF_ID] = ORF_seq;

                // iterate over nodes traversed by ORF
                for (const auto& node_traversed : (std::get<0>(ORF_node_paths.at(ORF_seq))))
                {
                    // add to triplet list, with ORRF_ID (row), node id (column) and set value as 1 (true)
                    tripletList.push_back(ET(ORF_ID, unitig_map.at(node_traversed).unitig_id, 1));
                }
                ORF_ID++;
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
                    auto ORF1_nodes = ((std::get<0>(ORF_node_paths.at(temp_ORF1)).size() >= std::get<0>(ORF_node_paths.at(temp_ORF2)).size()) ? ORF_node_paths.at(temp_ORF1) : ORF_node_paths.at(temp_ORF2));
                    auto ORF2_nodes = ((std::get<0>(ORF_node_paths.at(temp_ORF1)).size() >= std::get<0>(ORF_node_paths.at(temp_ORF2)).size()) ? ORF_node_paths.at(temp_ORF2) : ORF_node_paths.at(temp_ORF1));
                    auto ORF1 = ((std::get<0>(ORF_node_paths.at(temp_ORF1)).size() >= std::get<0>(ORF_node_paths.at(temp_ORF2)).size()) ? temp_ORF1 : temp_ORF2);
                    auto ORF2 = ((std::get<0>(ORF_node_paths.at(temp_ORF1)).size() >= std::get<0>(ORF_node_paths.at(temp_ORF2)).size()) ? temp_ORF2 : temp_ORF1);

                    // for debugging
                    //auto row = init.row();
                    //auto col = init.col();

                    // initialise overlap type
                    // n = no overlap
                    // u = unidirectional overlap (3' of first overlaps with 5' of second ->->)
                    // c = convergent overlap (3' of first overlaps with 3' of second -><-)
                    // d = divergent overlap (5' of first overlaps with 5' of second <-->)
                    // w = ORF lies completely within another
                    // i = incompatible overlap (i.e. if overlap is greater than maximum or same stop codon shared in same frame)
                    char overlap_type = 'n';

                    // set reversed variable to determine if ORFs are in same strand
                    bool reversed = false;

                    // set first ORF in relative ordering on strand
                    int first_ORF = 1;

                    // initialise absolute overlap
                    size_t abs_overlap = 0;

                    // initialise overlap_complete check
                    bool overlap_complete = false;

                    // get strand of 3p and 5p for reversal if necessary
                    bool ORF1_5p_strand = std::get<2>(ORF1_nodes)[0];
                    bool ORF1_3p_strand = std::get<2>(ORF1_nodes).back();
                    bool ORF2_5p_strand = std::get<2>(ORF2_nodes)[0];
                    bool ORF2_3p_strand = std::get<2>(ORF2_nodes).back();

                    // get ORF1 and ORF2 5' and 3' ends in pair <node_head_kmer, position>
                    std::pair<std::string, size_t> ORF1_5p(std::get<0>(ORF1_nodes)[0], std::get<0>(std::get<1>(ORF1_nodes)[0]));
                    std::pair<std::string, size_t> ORF1_3p(std::get<0>(ORF1_nodes).back(), std::get<1>(std::get<1>(ORF1_nodes).back()));
                    std::pair<std::string, size_t> ORF2_5p(std::get<0>(ORF2_nodes)[0], std::get<0>(std::get<1>(ORF2_nodes)[0]));
                    std::pair<std::string, size_t> ORF2_3p(std::get<0>(ORF2_nodes).back(), std::get<1>(std::get<1>(ORF2_nodes).back()));

                    // initialise values for start and end nodes
                    std::string ORF1_start_node = std::get<0>(ORF1_nodes)[0];
                    std::string ORF1_end_node = std::get<0>(ORF1_nodes).back();
                    std::string ORF2_start_node = std::get<0>(ORF2_nodes)[0];
                    std::string ORF2_end_node = std::get<0>(ORF2_nodes).back();

                    // work out if node 1 is negative by checking first node in pos_strand_map. If it doesn't match, it is negatively stranded.
                    bool negative = false;
                    if (ORF1_5p_strand != pos_strand_map.at(*colit).at(ORF1_start_node))
                    {
                        negative = true;
                    }

                    // work out if ORF2 is in same strand as ORF1, if so leave reversed as false. If not, set reversed as true
                    if ((ORF2_5p_strand != pos_strand_map.at(*colit).at(ORF2_start_node) && !negative) || (ORF2_5p_strand == pos_strand_map.at(*colit).at(ORF2_start_node) && negative))
                    {
                        reversed = true;
                    }

                    // if reversed is true, iterate through ORF2 coordinates and reverse
                    // check if both strands are the same. If not, reverse the nodes of ORF2 and their within-node coordinates
                    if (reversed)
                    {
                        // Reverse ORF2 node vector
                        std::reverse(std::get<0>(ORF2_nodes).begin(), std::get<0>(ORF2_nodes).end());
                        // Reverse ORF2 node coordinate vector
                        std::reverse(std::get<1>(ORF2_nodes).begin(), std::get<1>(ORF2_nodes).end());

                        // iterate over ORF2_nodes coordinate vector, reversing the coordinates relative to the end index of the node
                        for (auto & node_coord : std::get<1>(ORF2_nodes))
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

                        // correct ORF2 5p and 3p positions, reversing the coordinates
                        ORF2_5p.second = std::get<1>(std::get<1>(ORF2_nodes).back());
                        ORF2_3p.second = std::get<0>(std::get<1>(ORF2_nodes)[0]);
                        ORF2_start_node = std::get<0>(ORF2_nodes)[0];
                        ORF2_end_node = std::get<0>(ORF2_nodes).back();
                    }

                    // check if 3p matches between the two ORFs and node is in the same strand. If so, set as incompatible, and the overlap as the shorter of the two ORFS
                    if (ORF1_3p == ORF2_3p && ORF1_3p_strand == ORF2_3p_strand)
                    {
                        overlap_type = 'i';
                        overlap_complete = true;
                        abs_overlap = (ORF1.size() >= ORF2.size() ? ORF2.size() : ORF1.size()) - 16;

                        // set the first ORF to be the longest, regardless of strand
                        first_ORF = ((ORF1.size() >= ORF2.size()) ? 1 : 2);
                    }

                    // initialise vectors to capture overlapping node indices for ORFs
                    std::vector<size_t> ORF_1_overlap_node_index;
                    std::vector<size_t> ORF_2_overlap_node_index;

                    // initialise start and end index lists for ORF2 within ORF1
                    std::vector<size_t> start_index_list;
                    std::vector<size_t> end_index_list;

                    // if 3p ends don't match, look for indices

                    if (!overlap_complete)
                    {
                        // look for presence of ORF2 start node in ORF1 nodes
                        auto iter = find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_start_node);
                        while(iter != std::get<0>(ORF1_nodes).end())
                        {
                            size_t index = iter - std::get<0>(ORF1_nodes).begin();
                            start_index_list.push_back(index);
                            iter = find(iter + 1, std::get<0>(ORF1_nodes).end(), ORF2_start_node);
                        }

                        // look for presence of ORF2 end node in ORF1 nodes
                        iter = find(std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).end(), ORF2_end_node);
                        while (iter != std::get<0>(ORF1_nodes).end())
                        {
                            size_t index = iter - std::get<0>(ORF1_nodes).begin();
                            end_index_list.push_back(index);
                            iter = find(iter + 1, std::get<0>(ORF1_nodes).end(), ORF2_end_node);
                        }

                        // if both start and end indexes are found, check if ORF2 sits within ORF1
                        if (!start_index_list.empty() && !end_index_list.empty())
                        {
                            for (const auto& start_index : start_index_list)
                            {
                                for (const auto& end_index : end_index_list)
                                {
                                    // check moving in forward direction
                                    std::vector<std::string> ORF1_nodes_sliced;
                                    if (start_index <= end_index)
                                    {
                                        // slice ORF1 node vector from first entry to last
                                        ORF1_nodes_sliced = std::vector<std::string> (std::get<0>(ORF1_nodes).begin() + start_index, std::get<0>(ORF1_nodes).begin() + end_index + 1);
                                    }
                                    // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
                                    if (ORF1_nodes_sliced == std::get<0>(ORF2_nodes))
                                    {
                                        for (size_t i1 = start_index; i1 < end_index + 1; i1++)
                                        {
                                            ORF_1_overlap_node_index.push_back(i1);
                                        }
                                        for (size_t i2 = 0; i2 < std::get<0>(ORF2_nodes).size(); i2++)
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
                    }
                    // if overlap not found, check if beginning of ORF2 overlaps with end of ORF1. May be case that start and end nodes are found, but not complete match across ORF1
                    if (!overlap_complete && !start_index_list.empty())
                    {
                        for (const auto& start_index : start_index_list)
                        {
                            // check moving in forward direction
                            std::vector<std::string> ORF1_nodes_sliced;
                            std::vector<std::string> ORF2_nodes_sliced;

                            // slice ORF1 node vector from first entry to last
                            ORF1_nodes_sliced = std::vector<std::string> (std::get<0>(ORF1_nodes).begin() + start_index, std::get<0>(ORF1_nodes).end());
                            // check if ORF1_slice is too large to slice ORF1 (meaning likely ORF2 is reversed
                            if (ORF1_nodes_sliced.size() <= std::get<0>(ORF2_nodes).size())
                            {
                                // slice ORF1 node vector from first entry to the equivalent length of ORF2_nodes_sliced
                                ORF2_nodes_sliced = std::vector<std::string>(std::get<0>(ORF2_nodes).begin(), std::get<0>(ORF2_nodes).begin() + ORF1_nodes_sliced.size());
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
                            std::vector<std::string> ORF1_nodes_sliced;
                            std::vector<std::string> ORF2_nodes_sliced;

                            // slice ORF1 node vector from first entry to last
                            ORF1_nodes_sliced = std::vector<std::string> (std::get<0>(ORF1_nodes).begin(), std::get<0>(ORF1_nodes).begin() + end_index + 1);
                            // check that slice is small enough to correctly slice ORF2_nodes
                            if (ORF1_nodes_sliced.size() <= std::get<0>(ORF2_nodes).size())
                            {
                                // slice ORF2 node vector from last entry - end_index to the equivalent length of ORF2_nodes_sliced
                                ORF2_nodes_sliced = std::vector<std::string>(std::get<0>(ORF2_nodes).end() - end_index - 1, std::get<0>(ORF2_nodes).end());
                            }

                            // check if sliced vectors are equivalent. If they are, add indexes to overlap indexes
                            if (ORF1_nodes_sliced == ORF2_nodes_sliced)
                            {
                                for (size_t i1 = 0; i1 < ORF1_nodes_sliced.size(); i1++)
                                {
                                    ORF_1_overlap_node_index.push_back(i1);
                                }
                                for (size_t i2 = std::get<0>(ORF2_nodes).size() - end_index - 1; i2 < std::get<0>(ORF2_nodes).size(); i2++)
                                {
                                    ORF_2_overlap_node_index.push_back(i2);
                                }
                                overlap_complete = true;

                                break;
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
                            size_t overlap_start;
                            size_t overlap_end;

                            for (size_t i = 0; i < ORF_1_overlap_node_index.size(); i++)
                            {
                                // get the index of the corresponding overlapping nodes between ORF1 and ORF2
                                size_t ORF1_overlap_node = ORF_1_overlap_node_index[i];
                                size_t ORF2_overlap_node = ORF_2_overlap_node_index[i];

                                // get the node coordinates traversed by each ORF
                                size_t ORF1_start = std::get<0>(std::get<1>(ORF1_nodes)[ORF1_overlap_node]);
                                size_t ORF1_end = std::get<1>(std::get<1>(ORF1_nodes)[ORF1_overlap_node]);
                                size_t ORF2_start = std::get<0>(std::get<1>(ORF2_nodes)[ORF2_overlap_node]);
                                size_t ORF2_end = std::get<1>(std::get<1>(ORF2_nodes)[ORF2_overlap_node]);
                                size_t node_size = std::get<2>(std::get<1>(ORF2_nodes)[ORF2_overlap_node]);

                                // get the first node involved in the overlap
                                std::string overlap_node = std::get<0>(ORF1_nodes)[ORF1_overlap_node];

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

                                        // calculate overlap, negating node overlap. Add one as zero indexed.
                                        abs_overlap += ((high_index - low_index) - node_overlap) + 1;
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
                            ORF_overlap_map_private[*colit][ORF2][ORF1] = overlap_tuple;
                        } else {
                            std::pair<char, size_t> overlap_tuple(overlap_type, abs_overlap);
                            ORF_overlap_map_private[*colit][ORF1][ORF2] = overlap_tuple;
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            // merge all thread maps
            full_ORF_map.insert(full_ORF_map_private.begin(), full_ORF_map_private.end());
            ORF_overlap_map.insert(ORF_overlap_map_private.begin(), ORF_overlap_map_private.end());
        }
    }

    auto ORF_overlap_pair = std::make_pair(ORF_overlap_map, full_ORF_map);
    return ORF_overlap_pair;
}
