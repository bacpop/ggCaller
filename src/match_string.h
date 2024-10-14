#ifndef MATCH_STRING_H
#define MATCH_STRING_H

#include "unitigDict.h"
#include "kseq.h"

// match_strings
// code from
// https://stackoverflow.com/questions/735204/convert-a-string-in-c-to-upper-case
char ascii_toupper_char(char c);

std::pair<fm_index_coll, std::vector<size_t>> index_fasta(const std::string& fasta_file,
                                                          const bool write_idx,
                                                          const std::string& path_dir);

std::pair<int, bool> seq_search(const std::string& query,
                                const fm_index_coll& ref_idx);

std::pair<ContigLoc, bool> get_ORF_coords(const std::string& query,
                                          const fm_index_coll& fm_idx,
                                          const std::vector<size_t>& contig_locs);

std::vector<int> reverse_unitig_path(const std::vector<int>& unitig_path);

std::pair<bool, bool> path_search(const std::vector<int>& query_path,
                                  const fm_index_coll& ref_idx);

#endif //MATCH_STRING_H
