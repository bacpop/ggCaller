#ifndef TRANSLATION_H
#define TRANSLATION_H

#include "definintions.h"

class translate
{
    public:
    // take a DNA string and translate
    std::string translate (const std::string dna_seq)
    {
        std::string aa_seq;
        for (size_t i = 0; i < dna_seq.size(); i += 3)
        {
            const std::string codon = dna_seq(i, 3);
            aa_seq += tMap_[codon];
        }
        return aa_seq;
    }
    private:
        robin_hood::unordered_map<std::string, char>> tMap_
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
            {"CGU", 'R'},
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
            {"GGG", 'G')
        };
}


#endif //TRANSLATION_H
