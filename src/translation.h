#ifndef TRANSLATION_H
#define TRANSLATION_H

#include "definitions.h"

class translate
{
    public:
    // take a DNA string and translate
    translate (const std::string dna_seq)
    {
        for (size_t i = 0; i < dna_seq.size(); i += 3)
        {
            const std::string codon = dna_seq.substr(i, 3);
            aa_seq_ += tMap_[codon];
            // break if stop codon present
            if (tMap_[codon] == '*') {
                break;
            }
        }
    };
    std::string aa () {return aa_seq_;};
    private:
        robin_hood::unordered_map<std::string, char> tMap_ =
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
    std::string aa_seq_;
};

// balrog tokenisation
class tokenise
{
public:
    // take a DNA string and translate
    tokenise (const std::string& aa_seq)
    {
        for (size_t i = 0; i < aa_seq.size(); i++)
        {
            tokens_.push_back(tMap_[aa_seq.at(i)]);
        }
    };
    std::vector<int> out () {return tokens_;};
private:
    robin_hood::unordered_map<char, int> tMap_ =
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
    std::vector<int> tokens_;
};

class nuc_encode
{
public:
    // take a DNA string and translate
    nuc_encode (const char seq)
    {
        token_ = tMap_[seq];
    };
    int out () {return token_;};
private:
    robin_hood::unordered_map<char, int> tMap_ =
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
    int token_;
};

class start_encode
{
public:
    // take a DNA string and translate
    start_encode (const std::string& seq)
    {
        token_ = tMap_[seq];
    };
    int out () {return token_;};
private:
    robin_hood::unordered_map<std::string, int> tMap_ =
            {
                {"ATG", 0},
                {"GTG", 1},
                {"TTG", 2}
            };
    int token_;
};



#endif //TRANSLATION_H
