//
// Created by Sam Horsfield on 25/08/2020.
//

#include "match_string.h"

//run fmindex workflow
std::unordered_map<std::string, bool> call_strings(const std::vector<std::string>& assembly_list,
                                                   const std::unordered_map<std::string, std::string>& query_dict,
                                                   const bool write_idx,
                                                   const size_t num_threads)
{
    // Create threaded queue for computation
    assert(num_threads >= 1);
    const unsigned long int database_size = (unsigned long int)assembly_list.size();
    const unsigned long int calc_per_thread = database_size / num_threads;
    const unsigned int num_big_threads = database_size % num_threads;

    size_t start = 0;
    std::vector<size_t> start_points;
    for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) // Loop over threads
    {
        start_points.push_back(start);

        // First 'big' threads have an extra job
        if (thread_idx < num_big_threads)
        {
            start += calc_per_thread + 1;
        }
        else
        {
            start += calc_per_thread;
        }
    }
    start_points.push_back(start);

    // Read all sequences into memory as Fasta objects (threaded)
    std::cerr << "Constructing indexes for all input sequences..." << std::endl;
    std::vector<fm_index_coll> seq_idx;

    // multithread with openMP
    #pragma omp parallel num_threads(num_threads)
    {
        std::vector<fm_index_coll> seq_idx_private;
        #pragma omp for nowait
        for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx)
        {
            seq_idx_private = std::move(index_fasta(assembly_list, start_points[thread_idx], start_points[thread_idx + 1], write_idx));
        }
        #pragma omp critical
        seq_idx.insert(seq_idx.end(), seq_idx_private.begin(), seq_idx_private.end());
    }

    std::cerr << "Searching for gene sequences..." << std::endl;

    // Run searches using openMP, looping over genes and multithreading searches across all fmindexes
    // create unordered map object for output of query
    std::unordered_map<std::string, bool> results;

    for (auto& entry : query_dict)
    {
        // initialise hits variable for specific gene
        int hits = 0;
        // convert string to dn5 vector
        seqan3::dna5_vector query = entry.second | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;

        #pragma omp parallel for reduction(+:hits) num_threads(num_threads)
        for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx)
        {
            hits += seq_search(query, seq_idx, start_points[thread_idx], start_points[thread_idx + 1]);
        }

        //set present to true if matches are detected
        bool present = false;
        if (hits)
        {
            present = true;
        }

        // Print results for query to unordered map
        results[entry.first] = present;
    }

    std::cerr << "Done." << std::endl;
    return results;
}

// index fasta files
std::vector<fm_index_coll> index_fasta(const std::vector<std::string>& fasta_files,
                                       const size_t start,
                                       const size_t end,
                                       const bool& write_idx)
    {
    std::vector<fm_index_coll> seq_idx;
    for (auto file_it = fasta_files.begin() + start; file_it != fasta_files.begin() + end; ++file_it)
    {
        fm_index_coll ref_index;
        // create fm index file name
        std::string idx_file_name = *file_it + ".fm";

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
            seqan3::sequence_file_input reference_in{*file_it};
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
        seq_idx.push_back(ref_index);
    }
    return seq_idx;
}

//search for a specific sequence within an fm index array
int seq_search(const seqan3::dna5_vector& query,
               const std::vector<fm_index_coll>& seq_idx,
               const size_t start,
               const size_t end)
{
    int present = 0;
    int query_count = 0;
    for (auto ref_it = seq_idx.begin() + start; ref_it != seq_idx.begin() + end; ref_it++)
    {
        //count number of occurrences
        auto results = search(query, *ref_it);
        query_count = (int)std::ranges::distance(results);

        if (query_count == 0)
        {
            //find reverse complement, count number of occurrences
            auto results = search(query | std::views::reverse | seqan3::views::complement, *ref_it);
            query_count = (int)std::ranges::distance(results);
        }

        if (query_count != 0) {
            // debug_stream << "found" << std::endl;
            int present = 1;
            //std::cerr << present << std::endl;
            return present;
        }
        //std::cerr << "Not found!" << std::endl;
    }
    return present;
}