#ifndef UNITIG_DICT_H
#define UNITIG_DICT_H

#include "definitions.h"
#include <boost/serialization/access.hpp>

// unitigDict class declaration
class MyUnitigMap : public CCDBG_Data_t<MyUnitigMap>, CDBG_Data_t<MyUnitigMap> {
    public:

    // Clear method for CompactedDBG
    void clear(const UnitigColorMap<MyUnitigMap>& um_dest);

    // Concatenation method for ColoredCDBG
    void concat(const UnitigColorMap<MyUnitigMap>& um_dest, const UnitigColorMap<MyUnitigMap>& um_src);

    // Extraction method for ColoredCDBG
    void extract(const UnitigColors* uc_dest, const UnitigColorMap<MyUnitigMap>& um_src, const bool last_extraction);

    // serialisation method
    string serialize(const const_UnitigColorMap<MyUnitigMap>& um_src) const;

    // add and return untig kmers and ids
    void set_id (size_t new_id) {_id = new_id;};
    size_t get_id () const {return _id;};

    // add codon information
    void add_codon (const bool forward, const std::bitset<3>& array);
    void add_codon (const bool forward, const std::bitset<9>& array);
    void set_forward_stop (bool choice) {_forward_stop = choice;};
    void set_reverse_stop (bool choice) {_reverse_stop = choice;};
    bool forward_stop () const {return _forward_stop;};
    bool reverse_stop () const {return _reverse_stop;};

    // get codon information
    std::bitset<3> get_codon_arr (const bool full, const bool forward, const int frame) const;

    // add and return unitig colours
    void add_full_colour(boost::dynamic_bitset<> colours) {_unitig_full_colour = colours;};
    boost::dynamic_bitset<> full_colour() const {return _unitig_full_colour;}

    // access end_contig
    void set_end_contig (int colour, size_t nb_colours);
    bool end_contig (int colour) const;

    private:
    // Bifrost data accessor methods
    inline uint8_t get() const { return da_id; };
    inline void set(const uint8_t id) { da_id = id; };
    // Bifrost data accessor ID
    uint8_t da_id = 0;

    // custom unitig ID
    int _id;

    // codon arrays, initialise with two strands and 3 frames for each (6 reading frames total)
    // for full codon, only require 6 bits, 3 for forward, 3 for reverse in 0th frame
    // for part codons from 0th position - first 9 are forward strand, second 9 are reverse.
    // Of each 9, first 3 = 0th frame, second 3 = 1st frame, third 3 = 2nd frame
    std::bitset<6> _full_codon;
    std::bitset<18> _part_codon;

    // unitig colours
    boost::dynamic_bitset<> _unitig_full_colour;

    // bitset to determine if unitig is end of a contig/assembly
    boost::dynamic_bitset<> _end_contig;

    // forward_stop presence/absence
    bool _forward_stop = false;
    bool _reverse_stop = false;
};

// tuple of GraphVector, a mapping of colours to component nodes, the number of colours and the size of the overlap
typedef std::tuple<std::vector<std::string>, size_t, int, std::vector<bool>> GraphTuple;

// function to get um and pointer to unitig map data
std::pair<UnitigColorMap<MyUnitigMap>, MyUnitigMap*> get_um_data (ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                  const std::vector<Kmer>& head_kmer_arr,
                                                                  const int id);


// const qualified function to get pointer to unitig map data
std::pair<const_UnitigMap<DataAccessor<MyUnitigMap>, DataStorage<MyUnitigMap>>, MyUnitigMap*> get_um_data (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                                                           const std::vector<Kmer>& head_kmer_arr,
                                                                                                           const int id);

// get unitig data from head-kmer
std::pair<UnitigColorMap<MyUnitigMap>, MyUnitigMap*> get_um_data (ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                  const Kmer& head_kmer);

// const qualified function to get unitig data from head-kmer
std::pair<const_UnitigMap<DataAccessor<MyUnitigMap>, DataStorage<MyUnitigMap>>, MyUnitigMap*> get_um_data (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                                                           const Kmer& head_kmer);

// non-member functions for generating sequence from graph
std::string generate_sequence_nm(const std::vector<int>& nodelist,
                                 const std::vector<indexPair>& node_coords,
                                 const size_t& overlap,
                                 const ColoredCDBG<MyUnitigMap>& ccdbg,
                                 const std::vector<Kmer>& head_kmer_arr);

// non member function to remove overlap ends of ORFNodeVector
void simplify_ORFNodeVector (ORFNodeVector& ORF_info,
                             const int& overlap);

#endif //UNITIG_DICT_H