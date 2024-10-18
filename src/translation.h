#ifndef TRANSLATION_H
#define TRANSLATION_H

#include "definitions.h"

std::string translate (const std::string& dna_seq);

std::vector<int> tokenise (const std::string& aa_seq);

int nuc_encode (const char seq);

int start_encode (const std::string& seq);


#endif //TRANSLATION_H
