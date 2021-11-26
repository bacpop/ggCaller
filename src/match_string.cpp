// ggCaller header
#include "match_string.h"

// index fasta files
std::pair<fm_index_coll, std::vector<size_t>> index_fasta(const std::string& fasta_file,
                                                           const bool& write_idx)
{
    fm_index_coll ref_index;

    // create fm index file name
    std::string idx_file_name = fasta_file + ".fm";

    // create entry for start and end of contigs within fm_index
    std::vector<size_t> contig_locs;

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
            contig_locs.push_back(reference_seq.size());
            reference_seq += ",";
        }

        // destroy seq and fp objects
        kseq_destroy(seq);
        gzclose(fp);

        sdsl::construct_im(ref_index, reference_seq.c_str(), 1); // generate index
        if (write_idx)
        {
            store_to_file(ref_index, idx_file_name); // save it
        }
    } else
    {
        auto locations = sdsl::locate(ref_index, ",");
        sort(locations.begin(), locations.end());
        contig_locs.insert(contig_locs.end(), locations.begin(), locations.end());
    }

    return {ref_index, contig_locs};
}

//search for a specific sequence within an fm index array
std::pair<int, bool> seq_search(const std::string& query,
                                const fm_index_coll& ref_idx)
{
    int query_loc = -1;
    //count number of occurrences in positive strand
    auto locations = sdsl::locate(ref_idx, query);

    // determine if sequence reversed
    bool rev_comp = false;

    // if not found, check reverse strand
    if (locations.empty())
    {
        const std::string rev_query = reverse_complement(query);
        locations = sdsl::locate(ref_idx, rev_query);
        rev_comp = true;
    }

    // take first entry from locations
    if (!locations.empty())
    {
        // debug_stream << "found" << std::endl;
        sort(locations.begin(), locations.end());
        query_loc = locations[0];
    }
    return {query_loc, rev_comp};
}

// determine true colours of sequence
std::pair<ContigLoc, bool> check_colours(const std::string& query,
                                           const fm_index_coll& fm_idx,
                                           const std::vector<size_t>& contig_locs)
{
    // initialise location pair
    ContigLoc contig_loc;

    const auto query_pair = seq_search(query, fm_idx);
    const int& query_loc = std::get<0>(query_pair);
    const int& rev_comp = std::get<1>(query_pair);

    //if sequence present then determine contig coordinates
    if (query_loc > 0)
    {
        // go through contig_locs to determine in which contig sequence sits
        for (int i = 0; i < contig_locs.size(); i++)
        {
            if (query_loc < contig_locs.at(i))
            {
                if (i == 0)
                {
                    contig_loc = {1, {query_loc + 1, query_loc + query.size()}};
                } else
                {
                    size_t relative_loc = (query_loc - contig_locs.at(i - 1));
                    contig_loc = {i + 1, {relative_loc, relative_loc + query.size() - 1}};
                }
                break;
            }
        }
    }

    return {contig_loc, rev_comp};
}