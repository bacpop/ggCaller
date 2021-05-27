#ifndef MATCH_STRING_H
#define MATCH_STRING_H

#include "ggCaller_class.h"

// match_strings
fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx);

int seq_search(const seqan3::dna5_vector& query,
               const fm_index_coll& ref_idx);

bool check_colours(const std::string& path_sequence,
                   const fm_index_coll& fm_idx);

#endif //MATCH_STRING_H
