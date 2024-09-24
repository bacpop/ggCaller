#include "translation.h"

robin_hood::unordered_map<std::string, char> codonMap_ =
{
    {"TTT", 'F'},
    {"TTC", 'F'},
    {"TTA", 'L'},
    {"TTG", 'L'},
    {"TCT", 'S'},
    {"TCC", 'S'},
    {"TCA", 'S'},
    {"TCG", 'S'},
    {"TAT", 'Y'},
    {"TAC", 'Y'},
    {"TAA", '*'},
    {"TAG", '*'},
    {"TGT", 'C'},
    {"TGC", 'C'},
    {"TGA", '*'},
    {"TGG", 'W'},
    {"CTT", 'L'},
    {"CTC", 'L'},
    {"CTA", 'L'},
    {"CTG", 'L'},
    {"CCT", 'P'},
    {"CCC", 'P'},
    {"CCA", 'P'},
    {"CCG", 'P'},
    {"CAT", 'H'},
    {"CAC", 'H'},
    {"CAA", 'Q'},
    {"CAG", 'Q'},
    {"CGT", 'R'},
    {"CGC", 'R'},
    {"CGA", 'R'},
    {"CGG", 'R'},
    {"ATT", 'I'},
    {"ATC", 'I'},
    {"ATA", 'I'},
    {"ATG", 'M'},
    {"ACT", 'T'},
    {"ACC", 'T'},
    {"ACA", 'T'},
    {"ACG", 'T'},
    {"AAT", 'N'},
    {"AAC", 'N'},
    {"AAA", 'K'},
    {"AAG", 'K'},
    {"AGT", 'S'},
    {"AGC", 'S'},
    {"AGA", 'R'},
    {"AGG", 'R'},
    {"GTT", 'V'},
    {"GTC", 'V'},
    {"GTA", 'V'},
    {"GTG", 'V'},
    {"GCT", 'A'},
    {"GCC", 'A'},
    {"GCA", 'A'},
    {"GCG", 'A'},
    {"GAT", 'D'},
    {"GAC", 'D'},
    {"GAA", 'E'},
    {"GAG", 'E'},
    {"GGT", 'G'},
    {"GGC", 'G'},
    {"GGA", 'G'},
    {"GGG", 'G'}
};

// translate dna sequence
std::string translate (const std::string& dna_seq)
{
    std::string aa_seq;
    for (size_t i = 0; i < dna_seq.size(); i += 3)
    {
        const std::string codon = dna_seq.substr(i, 3);
        aa_seq += codonMap_[codon];
        if (codonMap_[codon] == '*') {
            break;
        }
    }
    return aa_seq;
};

// balrog tokenisation
robin_hood::unordered_map<char, int> aatokenMap_ =
{
        {'L', 1},
        {'V', 2},
        {'I', 3},
        {'M', 4},
        {'C', 5},
        {'A', 6},
        {'G', 7},
        {'S', 8},
        {'T', 9},
        {'P', 10},
        {'F', 11},
        {'Y', 12},
        {'W', 13},
        {'E', 14},
        {'D', 15},
        {'N', 16},
        {'Q', 17},
        {'K', 18},
        {'R', 19},
        {'H', 20},
        {'*', 0},
        {'X', 0}
};

// tokenise aa sequence for BALROG
std::vector<int> tokenise (const std::string& aa_seq)
{
    std::vector<int> tokens_;
    for (size_t i = 0; i < aa_seq.size(); i++)
    {
        tokens_.push_back(aatokenMap_[aa_seq.at(i)]);
    }
    return tokens_;
};

robin_hood::unordered_map<char, int> nucMap_ =
{
    {'A', 0},
    {'T', 1},
    {'G', 2},
    {'C', 3},
    {'N', 0},
    {'M', 0},
    {'R', 0},
    {'Y', 0},
    {'W', 0},
    {'K', 0}
};

// take a DNA char and encode
int nuc_encode (const char seq)
{
    return nucMap_[seq];
};

robin_hood::unordered_map<std::string, int> startMap_ =
{
    {"ATG", 0},
    {"GTG", 1},
    {"TTG", 2}
};

// take a start codon string and encode
int start_encode (const std::string& seq)
{
    return startMap_[seq];
};
