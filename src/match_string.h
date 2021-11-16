#ifndef MATCH_STRING_H
#define MATCH_STRING_H

#include "unitigDict.h"
#include "kseq.h"

// match_strings
std::pair<fm_index_coll, std::vector<size_t>> index_fasta(const std::string& fasta_file,
                                                          const bool& write_idx);

std::pair<int, bool> seq_search(const std::string& query,
                                const fm_index_coll& ref_idx);

std::pair<ContigLoc, bool> check_colours(const std::string& query,
                                           const fm_index_coll& fm_idx,
                                           const std::vector<size_t>& contig_locs);

#endif //MATCH_STRING_H
