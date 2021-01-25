// ggCaller header
#include "ggCaller_classes.h"

// generate ORFs from paths
ORFNodeMap generate_ORFs(const ColoredCDBG<>& ccdbg,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const std::vector<std::pair<std::string, bool>>& unitig_path,
                         const int& overlap,
                         const size_t& min_len)
{
    // initialise path sequence and ORF list
    std::string path_sequence;
    std::vector<std::string> ORF_list;

    // here generate vector/map which contains the start/end of each node in basepairs. This will then be used to determine which nodes an ORF traverses
    std::vector<std::string> nodelist;
    std::vector<std::vector<size_t>> node_ranges;
    size_t node_start = 0;

    // generate path sequence by merging nodes in sequence
    for (const auto& node : unitig_path)
    {
        const Kmer kmer = Kmer(node.first.c_str());
        auto unitig = ccdbg.find(kmer, true);
        unitig.strand = node.second;

        // add node to node list for path coordinates
        nodelist.push_back(node.first);
        // 0th entry is start index of node within the unitig, 1st entry is end index of node within unitig, 2nd entry is node length.
        std::vector<size_t> node_range(3);

        // if strand is negative, calculate reverse complement
        std::string unitig_seq;
        if (unitig.strand)
        {
            unitig_seq = unitig.referenceUnitigToString();
        } else {
            unitig_seq = reverse_complement(unitig.referenceUnitigToString());
        }

        // calculate length of unitig and get end coordinates
        size_t node_end = unitig_seq.size();

        if (path_sequence.empty())
        {
            path_sequence = unitig_seq;
            // start index of node in path
            node_range[0] = node_start;
            // start index of next node (one past the end index)
            node_range[1] = node_end;
            // absolute end index of node
            node_range[2] = node_end - 1;
            node_start = node_end;
        } else {
            path_sequence.append(unitig_seq.begin()+overlap, unitig_seq.end());
            // start index of node in path
            node_range[0] = node_start - overlap;
            // absolute end index of node
            node_range[2] = node_end - 1;
            node_end += node_start - overlap;
            // start index of next node in path (one past the end index)
            node_range[1] = node_end;
            node_start = node_end;

        }
        node_ranges.push_back(node_range);
    }

    // generate codon indexes using FindIndex function for start and stop codons
    std::vector<size_t> start_codon_indices;
    std::vector<size_t> stop_codon_indices;
    std::vector<std::size_t> found_indices;

    for (const auto& codon : stop_codons)
    {
        found_indices = findIndex(path_sequence, codon, 0, 0, false);
        stop_codon_indices.insert(stop_codon_indices.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();
    }
    for (const auto& codon : start_codons)
    {
        found_indices = findIndex(path_sequence, codon, 0, 0, false);
        start_codon_indices.insert(start_codon_indices.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();
    }

    // create pair for each paired start and stop codon
    std::vector<std::pair<std::size_t, std::size_t>> ORF_index_pairs;

    // sort codon index arrays into ascending order
    std::sort(start_codon_indices.begin(), start_codon_indices.end());
    std::sort(stop_codon_indices.begin(), stop_codon_indices.end());

    // generate dictionaries for start and stop codon indices for each frame
    std::unordered_map<size_t, std::vector<size_t>> start_codon_dict;
    std::unordered_map<size_t, std::vector<size_t>> stop_codon_dict;

    std::vector<size_t> frame1;
    std::vector<size_t> frame2;
    std::vector<size_t> frame3;

    for (auto& index : start_codon_indices)
    {
        if (index % 3 == 0 || index == 0)
        {
            frame1.push_back(index);
        } else if (index % 3 == 1 || index == 1)
        {
            frame2.push_back(index);
        } else if (index % 3 == 2 || index == 2)
        {
            frame3.push_back(index);
        }
    }
    // update start codon indexes
    start_codon_dict[0] = std::move(frame1);
    start_codon_dict[1] = std::move(frame2);
    start_codon_dict[2] = std::move(frame3);
    // clear all index vectors
    frame1.clear();
    frame2.clear();
    frame3.clear();
    start_codon_indices.clear();

    for (auto& index : stop_codon_indices)
    {
        if (index % 3 == 0 || index == 0)
        {
            frame1.push_back(index);
        } else if (index % 3 == 1 || index == 1)
        {
            frame2.push_back(index);
        } else if (index % 3 == 2 || index == 2)
        {
            frame3.push_back(index);
        }
    }
    // update stop codon indexes
    stop_codon_dict[0] = std::move(frame1);
    stop_codon_dict[1] = std::move(frame2);
    stop_codon_dict[2] = std::move(frame3);
    // clear all index vectors
    frame1.clear();
    frame2.clear();
    frame3.clear();
    stop_codon_indices.clear();

    // iterate through frames, pair sequential start+stop codons after first stop codon
    for (int modulus = 0; modulus < 3; modulus++)
    {
        if (!stop_codon_dict[modulus].empty())
        {
            // initialise first stop codon
            size_t last_index = stop_codon_dict[modulus][0];

            // cycle through start codons
            for (const auto& start_index : start_codon_dict[modulus])
            {
                // check that start index is in correct frame and further in sequence than last stop codon index
                if (start_index > last_index)
                {
                    // cycle through stop codons
                    for (const auto& stop_index : stop_codon_dict[modulus])
                    {
                        // check that stop index is in correct frame and further in sequence than last start codon index
                        if (stop_index > start_index)
                        {
                            // stop index is indexed at first base, therefore end of ORF is two bases after
                            ORF_index_pairs.push_back(std::make_pair(start_index, stop_index + 2));
                            last_index = start_index;
                            break;
                        }
                    }
                }
            }
        }
    }
    // ORF dictionary for locating ORF position in graph. Nodes traversed by ORF is ordered map to ensure ordering is kept for overlap comparisons.
    ORFNodeMap ORF_map;

    // generate sequences for ORFs from codon pairs
    for (const auto& codon_pair : ORF_index_pairs)
    {
        // add one as codon_pair is zero-indexed
        size_t ORF_len = (codon_pair.second - codon_pair.first) + 1;
        if (ORF_len >= min_len)
        {
            // make a pair containing a vector of each node name, and corresponding vector of positions traversed in the node
            ORFNodeVector ORF_node_vector;

            // pull 16bp upstream of start codon for TIS model if possible
            if (codon_pair.first >= 16) {
                for (size_t i = 0; i < nodelist.size(); i++)
                {
                    size_t traversed_node_start;
                    size_t traversed_node_end;
                    bool start_assigned = false;
                    bool end_assigned = false;

                    if (codon_pair.first < node_ranges[i][0]){
                        traversed_node_start = 0;
                        start_assigned = true;
                    } else if (codon_pair.first >= node_ranges[i][0] && codon_pair.first < node_ranges[i][1]){
                        traversed_node_start = codon_pair.first - node_ranges[i][0];
                        // check that the start difference to end is greater than the overlap. If not, then sequence is covered in next node traversal
                        if ((node_ranges[i][2] - traversed_node_start) >= overlap)
                        {
                            start_assigned = true;
                        }
                    }

                    if (codon_pair.second >= node_ranges[i][1]) {
                        traversed_node_end = node_ranges[i][2];
                        end_assigned = true;
                    } else if (codon_pair.second >= node_ranges[i][0] && codon_pair.second < node_ranges[i][1]){
                        traversed_node_end = codon_pair.second - node_ranges[i][0];
                        // check that the end is greater than the overlap. If not, then sequence is already covered in prior node traversal
                        if (traversed_node_end >= overlap)
                        {
                            end_assigned = true;
                        }
                    }

                    if (start_assigned && end_assigned){
                        size_t node_end = node_ranges[i][2];
                        indexTriplet node_coords = std::make_tuple(traversed_node_start, traversed_node_end, node_end);
                        ORF_node_vector.first.push_back(nodelist[i]);
                        ORF_node_vector.second.push_back(std::move(node_coords));
                    }

                }
                std::string ORF = path_sequence.substr((codon_pair.first - 16), (ORF_len + 16));
                //std::string ORF = path_sequence.substr((codon_pair.first), (ORF_len));
                ORF_map.emplace(std::move(ORF), std::move(ORF_node_vector));
            }
        }
    }
    return ORF_map;
}

std::pair<ORFColoursMap, std::vector<std::string>> filter_artificial_ORFS(StrandSeqORFMap& all_ORFs,
                                                                          StrandORFNodeMap& ORF_node_paths,
                                                                          const std::vector<std::string>& fasta_files,
                                                                          const bool write_index)
{
    // call strings edits all_ORFs in place
    // run call strings for positive strand
    call_strings(all_ORFs["+"], "+", ORF_node_paths["+"], fasta_files, write_index);
    // run call strings for negative strand
    call_strings(all_ORFs["-"], "-", ORF_node_paths["-"], fasta_files, write_index);

    // generate a colours dictionary for gene overlap analysis
    auto return_tuple = sort_ORF_colours(all_ORFs);
    return return_tuple;
}

std::pair<ORFColoursMap, std::vector<std::string>> sort_ORF_colours(const StrandSeqORFMap& all_ORFs)
{
    ORFColoursMap ORF_colours_map;
    std::unordered_set<std::string> ORF_colours_set;
    for (const auto& ORF_strand : all_ORFs)
    {
        for (const auto& ORF_seq : ORF_strand.second)
        {
            std::string colour_key;
            for (const auto colour : ORF_seq.second)
            {
                colour_key += std::to_string(colour);
            }
            ORF_colours_map[colour_key].push_back(ORF_seq.first);
            ORF_colours_set.insert(colour_key);
        }
    }
    // convert set to vector to enable parallelisation
    std::vector<std::string> ORF_colours_vector(ORF_colours_set.begin(), ORF_colours_set.end());

    // generate return tuple
    std::pair<ORFColoursMap, std::vector<std::string>> return_tuple = std::make_pair(ORF_colours_map, ORF_colours_vector);
    return return_tuple;
}

std::tuple<StrandSeqORFMap, StrandORFNodeMap> call_ORFs(const ColoredCDBG<>& ccdbg,
                                                  const PathTuple& path_tuple,
                                                  const std::vector<std::string>& stop_codons_for,
                                                  const std::vector<std::string>& start_codons_for,
                                                  const int& overlap,
                                                  const size_t& min_ORF_length)
{
//initialise all_ORFs
    StrandSeqORFMap all_ORFs;
    StrandORFNodeMap ORF_node_paths;

    #pragma omp parallel
    {
        StrandSeqORFMap all_ORFs_private;
        StrandORFNodeMap ORF_node_paths_private;

        #pragma omp for nowait
        // iterate over head_kmer_strings
        for (auto it = std::get<1>(path_tuple).begin(); it < std::get<1>(path_tuple).end(); it++)
        {
            const auto unitig_paths = std::get<0>(path_tuple).at(*it);
            // iterate over paths following head_kmer
            for (const auto& path : unitig_paths)
            {
                int i = 0;
                const std::string strand = (path.first[0].second ? "+" : "-");
                const std::vector<bool> path_colour = path.second;
                // iterate over each start codon, generate all ORFs within the path
                for (const auto& start_codon : start_codons_for)
                {
                    std::vector<std::string> start_codon_vector{start_codon};
                    auto ORF_map = generate_ORFs(ccdbg, stop_codons_for, start_codon_vector, path.first, overlap, min_ORF_length);

                    // check if item in all_ORFs already. If not, add colours array. If yes, update the colours array.
                    for (const auto& ORF : ORF_map)
                    {
                        ORF_node_paths_private[strand][ORF.first] = ORF.second;

                        if (all_ORFs_private[strand].find(ORF.first) == all_ORFs_private[strand].end())
                        {
                            all_ORFs_private[strand][ORF.first] = path_colour;
                        } else {
                            std::vector<bool> updated_colours = add_colours_array(all_ORFs_private[strand][ORF.first], path_colour);
                            all_ORFs_private[strand][ORF.first] = updated_colours;
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            // go through all private ORFs, update colours as before in all_ORFs
            ORF_node_paths["+"].insert(ORF_node_paths_private["+"].begin(), ORF_node_paths_private["+"].end());
            ORF_node_paths["-"].insert(ORF_node_paths_private["-"].begin(), ORF_node_paths_private["-"].end());

            for (const auto& ORF_strand : all_ORFs_private)
            {
                for (const auto& ORF_seq : ORF_strand.second)
                {
                    if (all_ORFs[ORF_strand.first].find(ORF_seq.first) == all_ORFs[ORF_strand.first].end())
                    {
                        all_ORFs[ORF_strand.first][ORF_seq.first] = ORF_seq.second;
                    } else {
                        std::vector<bool> updated_colours = add_colours_array(all_ORFs[ORF_strand.first][ORF_seq.first], ORF_seq.second);
                        all_ORFs[ORF_strand.first][ORF_seq.first] = updated_colours;
                    }
                }
            }
        }
    }
    const auto ORF_tuple = std::make_tuple(all_ORFs, ORF_node_paths);
    return ORF_tuple;
}

void write_to_file (const std::string& outfile_name,
                    const StrandSeqORFMap& all_ORFs)
{
    int gene_id = 1;
    ofstream outfile;
    outfile.open(outfile_name);

    //iterate over strands in all_ORFs
    for (const auto& strand_dict : all_ORFs)
    {
        // iterate over ORFs in all_ORFs
        for (const auto& ORF_dict : strand_dict.second)
        {
            // generate string for colours
            std::string colours;
            for (const auto& i : ORF_dict.second)
            {
                colours += std::to_string(i);
            }
            // append to file
            outfile << ">" << std::to_string(gene_id) << "_" << strand_dict.first << "_" << colours << "\n" << ORF_dict.first << "\n";
            gene_id++;
        }
    }
    outfile.close();
}
