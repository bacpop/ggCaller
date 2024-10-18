// ggCaller header
#include "match_string.h"

// code from
// https://stackoverflow.com/questions/735204/convert-a-string-in-c-to-upper-case
char ascii_toupper_char(char c) {
  return ('a' <= c && c <= 'z')
             ? c ^ 0x20
             : c; // ^ autovectorizes to PXOR: runs on more ports than paddb
}

// index fasta files
std::pair<fm_index_coll, std::vector<size_t>> index_fasta(const std::string& fasta_file,
                                                          const bool write_idx,
                                                          const std::string& path_dir)
{
    fm_index_coll ref_index;

    // create fm index file name
    const std::string base_filename = fasta_file.substr(fasta_file.find_last_of("/\\") + 1);

    std::string idx_file_name = path_dir + base_filename + ".fms";

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

        // convert string to uppercase to avoid indexing issues
        std::transform(reference_seq.begin(), reference_seq.end(), reference_seq.begin(), ::ascii_toupper_char);

        sdsl::construct_im(ref_index, reference_seq, 1); // generate index
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
    auto locations = sdsl::locate(ref_idx, query.begin(), query.end());

    // determine if sequence reversed
    bool strand = true;

    // if not found, check reverse strand
    if (locations.empty())
    {
        const std::string rev_query = reverse_complement(query);
        locations = sdsl::locate(ref_idx, rev_query.begin(), rev_query.end());
        strand = false;
    }

    // take first entry from locations
    if (!locations.empty())
    {
        // debug_stream << "found" << std::endl;
        sort(locations.begin(), locations.end());
        query_loc = locations[0];
    }
    return {query_loc, strand};
}

// determine true colours of sequence
std::pair<ContigLoc, bool> get_ORF_coords(const std::string& query,
                                          const fm_index_coll& fm_idx,
                                          const std::vector<size_t>& contig_locs)
{
    // initialise location pair
    ContigLoc contig_loc;

    const auto query_pair = seq_search(query, fm_idx);
    const int& query_loc = std::get<0>(query_pair);
    const bool strand = std::get<1>(query_pair);

    //if sequence present then determine contig coordinates
    if (query_loc >= 0)
    {
        // go through contig_locs to determine in which contig sequence sits
        for (int i = 0; i < contig_locs.size(); i++)
        {
            if (query_loc < contig_locs.at(i))
            {
                if (i == 0)
                {
                    contig_loc = {1, {query_loc, query_loc + query.size()}};
                } else
                {
                    size_t relative_loc = (query_loc - contig_locs.at(i - 1));
                    contig_loc = {i + 1, {relative_loc, relative_loc + (query.size() - 1)}};
                }
                break;
            }
        }
    }

    return {contig_loc, strand};
}

std::vector<int> reverse_unitig_path(const std::vector<int>& unitig_path)
{
    // fill reversed_path backwards, multiplying by -1
    std::vector<int> reversed_path(unitig_path.size());
    for (int i = 0; i < unitig_path.size(); i++)
    {
        reversed_path[(unitig_path.size() - 1) - i] = unitig_path[i] * - 1;
    }
    return reversed_path;
}

//search for a specific sequence within an fm index array
std::pair<bool, bool> path_search(const std::vector<int>& query_path,
                                 const fm_index_coll& ref_idx)
{
    bool present = false;

    // convert query into string
    std::ostringstream oss;
    std::copy(query_path.begin(), query_path.end(),
              std::ostream_iterator<int>(oss, ","));

    // add start delimeter
    std::string query = "," + oss.str();

    //count number of occurrences in positive strand
    auto count = sdsl::count(ref_idx, query);

    // determine if sequence reversed
    bool rev_comp = false;

    // if not found, check reverse strand
    if (count == 0)
    {
        const auto rev_query_path = reverse_unitig_path(query_path);
        oss.str(std::string());
        std::copy(rev_query_path.begin(), rev_query_path.end(),
                  std::ostream_iterator<int>(oss, ","));
        query = "," + oss.str();
        count = sdsl::count(ref_idx, query);
        rev_comp = true;
    }

    // take first entry from locations
    if (count > 0)
    {
        present = true;
    }
    return {present, rev_comp};
}