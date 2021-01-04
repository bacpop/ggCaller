// ggCaller header
#include "ggCaller_classes.h"

// generate ORFs from paths
std::vector<std::string> generate_ORFs(const ColoredCDBG<>& ccdbg,
                                       const std::vector<std::string>& stop_codons,
                                       const std::vector<std::string>& start_codons,
                                       const std::vector<std::pair<std::string, bool>>& unitig_path,
                                       const int& overlap,
                                       const size_t& min_len)
{
    // initialise path sequence and ORF list
    std::string path_sequence;
    std::vector<std::string> ORF_list;

    // generate path sequence by merging nodes in sequence
    for (const auto& node : unitig_path)
    {
        const Kmer kmer = Kmer(node.first.c_str());
        auto unitig = ccdbg.find(kmer, true);
        unitig.strand = node.second;

        // if strand is negative, calculate reverse complement
        std::string unitig_seq;
        if (unitig.strand)
        {
            unitig_seq = unitig.referenceUnitigToString();
        } else {
            unitig_seq = reverse_complement(unitig.referenceUnitigToString());
        }

        if (path_sequence.empty())
        {
            path_sequence = unitig_seq;
        } else {
            path_sequence.append(unitig_seq.begin()+overlap, unitig_seq.end());
        }
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
                            ORF_index_pairs.push_back(std::make_pair(start_index, stop_index));
                            last_index = start_index;
                            break;
                        }
                    }
                }
            }
        }
    }
    // generate sequences for ORFs from codon pairs
    for (const auto& codon_pair : ORF_index_pairs)
    {
        // generate ORF via substring from path_sequence, checking for ORF size
        size_t ORF_len = (codon_pair.second - codon_pair.first) + 3;
        if (ORF_len >= min_len)
        {
            std::string ORF = path_sequence.substr(codon_pair.first, ORF_len);
            ORF_list.push_back(std::move(ORF));

            // pull 16bp upstream of start codon for TIS model (not yet implemented)
            //if (codon_pair.first >= 16) {
                //std::string ORF = path_sequence.substr((codon_pair.first), (ORF_len));
                //std::string ORF = path_sequence.substr((codon_pair.first - 16), (ORF_len + 16));
                //ORF_list.push_back(std::move(ORF));
            //}
        }
    }
    return ORF_list;
}

void filter_artificial_ORFS(robin_hood::unordered_map<std::string, robin_hood::unordered_map<std::string, std::vector<bool>>>& all_ORFs,
                            const std::vector<std::string>& fasta_files,
                            const bool write_index = 1)
{
    // call strings edits all_ORFs in place
    // run call strings for positive strand
    call_strings(all_ORFs["+"], "+", fasta_files, write_index);
    // run call strings for negative strand
    call_strings(all_ORFs["-"], "-", fasta_files, write_index);
}

robin_hood::unordered_map<std::string, robin_hood::unordered_map<std::string, std::vector<bool>>> call_ORFs(const ColoredCDBG<>& ccdbg,
                                                                                                            const std::tuple<robin_hood::unordered_map<std::string, std::vector<std::pair<std::vector<std::pair<std::string, bool>>, std::vector<bool>>>>, std::vector<std::string>>& path_tuple,
                                                                                                            const std::vector<std::string>& stop_codons_for,
                                                                                                            const std::vector<std::string>& start_codons_for,
                                                                                                            const int& overlap,
                                                                                                            const size_t& min_ORF_length)
{
//initialise all_ORFs
    robin_hood::unordered_map<std::string, robin_hood::unordered_map<std::string, std::vector<bool>>> all_ORFs;

    #pragma omp parallel
    {
        robin_hood::unordered_map<std::string, robin_hood::unordered_map<std::string, std::vector<bool>>> all_ORFs_private;
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
                    std::vector<std::string> ORF_list = generate_ORFs(ccdbg, stop_codons_for, start_codon_vector, path.first, overlap, min_ORF_length);

                    // check if item in all_ORFs already. If not, add colours array. If yes, update the colours array.
                    for (const auto& ORF : ORF_list)
                    {
                        if (all_ORFs_private[strand].find(ORF) == all_ORFs_private[strand].end())
                        {
                            all_ORFs_private[strand][ORF] = path_colour;
                        } else {
                            std::vector<bool> updated_colours = add_colours_array(all_ORFs_private[strand][ORF], path_colour);
                            all_ORFs_private[strand][ORF] = updated_colours;
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            // go through all private ORFs, update colours as before in all_ORFs
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
    return all_ORFs;
}

void write_to_file (const std::string& outfile_name,
                    const robin_hood::unordered_map<std::string, robin_hood::unordered_map<std::string, std::vector<bool>>>& all_ORFs)
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
