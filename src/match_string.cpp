// ggCaller header
#include "match_string.h"

// index fasta files
fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx)
{
    fm_index_coll ref_index;
    fm_index_coll test_index;
    // create fm index file name
    std::string idx_file_name = fasta_file + ".fm";

    // if fm_index not available, generate it
    if (!load_from_file(ref_index, idx_file_name))
    {
        std::string reference_seq;

        // open the file handler
        gzFile fp = gzopen(fasta_file.c_str(), "r");

        if(fp == 0) {
            perror("fopen");
            exit(1);
        }
        // initialize seq
        kseq_t *seq = kseq_init(fp);

        // read sequence
        int l;
        while ((l = kseq_read(seq)) >= 0)
        {
            reference_seq += seq->seq.s;
            reference_seq += ",";
        }

        // destroy seq and fp objects
        kseq_destroy(seq);
        gzclose(fp);

        sdsl::construct(ref_index, reference_seq, 1); // generate index
        if (write_idx)
        {
            store_to_file(ref_index, idx_file_name); // save it
        }
    }

    return ref_index;
}

//search for a specific sequence within an fm index array
int seq_search(const std::string& query,
               const fm_index_coll& ref_idx)
{
    int present = 0;
    //count number of occurrences in positive strand
    size_t query_count = sdsl::count(ref_idx, query);

    // if not found, check reverse strand
    if (query_count == 0)
    {
        const std::string rev_query = reverse_complement(query);
        query_count = sdsl::count(ref_idx, rev_query);
    }

    if (query_count != 0)
    {
        // debug_stream << "found" << std::endl;
        present = 1;
    }
    return present;
}

// determine true colours of sequence
bool check_colours(const std::string& query,
                   const fm_index_coll& fm_idx)
{
    // initialise present
    bool present = true;

    //if sequence present then add to hits
    int hits = seq_search(query, fm_idx);
    if (!hits)
    {
        present = false;
    }

    return present;
}