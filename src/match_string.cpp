// ggCaller header
#include "ggCaller_classes.h"

// index fasta files
fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx)
{
    fm_index_coll ref_index;
    // create fm index file name
    std::string idx_file_name = fasta_file + ".fm";

    // Read index if it already exists
    if (std::filesystem::exists(idx_file_name))
    {
        std::ifstream is{idx_file_name, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(ref_index);
    }
    else
    {
        // Create index - should be viable for both single and collections fmindexes
        seqan3::sequence_file_input reference_in{fasta_file};
        std::vector<seqan3::dna5_vector> reference_seq;

        for (auto & [seq, id, qual] : reference_in) {
            reference_seq.push_back(std::move(seq));
        }

        // this is where the sequence is indexed
        ref_index = fm_index_coll{reference_seq};

        // write indexes to file if specified
        if (write_idx)
        {
            std::ofstream os{idx_file_name, std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(ref_index);
        }
    }
    return ref_index;
}

//search for a specific sequence within an fm index array
int seq_search(const seqan3::dna5_vector& query,
               const fm_index_coll& ref_idx)
{
    int present = 0;
    int query_count = 0;
    //count number of occurrences in positive strand
    auto results = search(query, ref_idx);
    query_count = (int)std::ranges::distance(results);

    // if not found, check reverse strand
    if (query_count == 0)
    {
        auto results = search(query | std::views::reverse | seqan3::views::complement, ref_idx);
        query_count = (int)std::ranges::distance(results);
    }

    if (query_count != 0)
    {
        // debug_stream << "found" << std::endl;
        present = 1;
    }
    return present;
}

////run fmindex workflow
//void call_strings(ORFNodeMap& ORF_node_paths,
//                  const std::vector<std::string>& assembly_list,
//                  const bool& write_idx)
//{
//    // compute number of colours
//    const size_t num_colours = assembly_list.size();
//
//    // generate list of genes for iteration
//    std::vector<std::string> ORF_list;
//    for (const auto& ORF : ORF_node_paths)
//    {
//        ORF_list.push_back(ORF.first);
//    }
//
//    // Read all sequences into memory as Fasta objects (threaded)
//    std::vector<fm_index_coll> seq_idx(num_colours);
//
//    // multithread with openMP
//    #pragma omp parallel for
//    for (size_t i = 0; i < num_colours; i++)
//    {
//        seq_idx[i] = std::move(index_fasta(assembly_list[i], write_idx));
//    }
//
//    // Run searches using openMP, looping over genes and multithreading searches of genes calls in specific fm-indexes
//    robin_hood::unordered_map<std::string, std::vector<bool>> to_remove;
//    #pragma omp parallel
//    {
//        robin_hood::unordered_map<std::string, std::vector<bool>> to_remove_private;
//        #pragma omp for nowait
//        for (auto itr = ORF_list.begin(); itr < ORF_list.end(); itr++)
//        {
//            // convert string to dn5 vector, get colours
//            seqan3::dna5_vector query = *itr | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;
//            auto& colours = ORF_node_paths.at(*itr).first;
//
//
//            //iterate over colours, if colour present then add to hits
//            for (size_t i = 0; i < num_colours; i++)
//            {
//                if (colours[i])
//                {
//                    int hits = seq_search(query, seq_idx[i]);
//                    if (!hits)
//                    {
//                        colours[i] = 0;
//                    }
//                }
//            }
//            // compute number of colours
//            int sum_colours = accumulate(colours.begin(), colours.end(), 0);
//
//            //if number of colours is 0, remove item from dictionary
//            if (sum_colours == 0)
//            {
//                ORF_node_paths.erase(*itr);
//            }
//        }
//    }
//}
//
////new workflow
//// generate list of fmindexes
//std::vector<fm_index_coll> generate_fmindex(const std::vector<std::string>& assembly_list,
//                                            const bool& write_idx)
//{
//    // compute number of colours
//    const size_t num_colours = assembly_list.size();
//
//    // Read all sequences into memory as Fasta objects (threaded)
//    std::vector<fm_index_coll> seq_idx(num_colours);
//
//    // multithread with openMP
//    #pragma omp parallel for
//    for (size_t i = 0; i < num_colours; i++)
//    {
//        seq_idx[i] = std::move(index_fasta(assembly_list[i], write_idx));
//    }
//
//    return seq_idx;
//}

// determine true colours of sequence
bool check_colours(const std::string& path_sequence,
                   const fm_index_coll& fm_idx)
{
    // initialise present
    bool present = true;

    // convert string to dn5 vector, get colours
    seqan3::dna5_vector query = path_sequence | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;

    //if sequence present then add to hits
    int hits = seq_search(query, fm_idx);
    if (!hits)
    {
        present = false;
    }
}