#ifndef MATCH_STRING_H
#define MATCH_STRING_H

#include "unitigDict.h"
#include "kseq.h"

// match_strings
fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx);

int seq_search(const std::string& query,
               const fm_index_coll& ref_idx);

bool check_colours(const std::string& query,
                   const fm_index_coll& fm_idx);

#endif //MATCH_STRING_H
