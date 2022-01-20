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
    void add_head(const std::string& head) {_head_kmer = head;};
    const std::string& head_kmer() const {return _head_kmer;};
    void set_id (size_t new_id) {_id = new_id;};
    size_t get_id () const {return _id;};

    // add codon information
    void add_codon (const bool& full, const bool& forward, const std::bitset<9>& array);
    void set_forward_stop (bool choice) {_forward_stop = choice;};
    void set_reverse_stop (bool choice) {_reverse_stop = choice;};
    bool forward_stop () const {return _forward_stop;};
    bool reverse_stop () const {return _reverse_stop;};

    // get codon information
    std::bitset<3> get_codon_arr (const bool full, const bool forward, const int frame) const;
//    std::vector<uint8_t> get_codon_dict (bool full, bool forward) const;

    // add contig mapping information
//    void add_contig_coords(const size_t& colour_ID, const std::tuple<size_t, size_t, size_t, size_t, bool>& contig_coords) {_contigcoords[colour_ID].push_back(contig_coords);};
//    const std::vector<std::tuple<size_t, size_t, size_t, size_t, bool>>& get_contig_coords(const size_t& colour_ID) const {return _contigcoords.at(colour_ID);};

    // add and return unitig colours
    void add_head_colour(boost::dynamic_bitset<> colours) {_unitig_head_colour = colours;};
    boost::dynamic_bitset<> head_colour() const {return _unitig_head_colour;};
    void add_full_colour(boost::dynamic_bitset<> colours) {_unitig_full_colour = colours;};
    boost::dynamic_bitset<> full_colour() const {return _unitig_full_colour;};
    void add_tail_colour(boost::dynamic_bitset<> colours) {_unitig_tail_colour = colours; _check_head_tail_equal();};
    boost::dynamic_bitset<> tail_colour() const {return _unitig_tail_colour;};
    bool head_tail_colours_equal() const {return _head_tail_colours_equal;};

    // access end_contig
    void set_end_contig (bool choice) {_end_contig = choice;};
    bool end_contig() const {return _end_contig;};

    //assign and return neighbours
//    void add_neighbour (bool strand, std::tuple<int, std::vector<uint8_t>, std::unordered_set<size_t>> neighbour) {_neighbours[strand].push_back(std::move(neighbour));};
//    void add_neighbour_colour (bool strand, int neighbour_ID, size_t colour_ID);

    // get neighbouring nodes
//    const NeighbourVector & get_neighbours (bool strand) const {return _neighbours[strand];};

    // assign traversing ORFs
    void set_ORFs (const size_t& colour_ID, const size_t& ORF_ID) {_traversing_ORFs[colour_ID].insert(ORF_ID);};
    bool ORFs_empty (const size_t& colour_ID) const {return (_traversing_ORFs.find(colour_ID) == _traversing_ORFs.end());};
    const std::unordered_set<size_t> & get_ORFs(const size_t& colour_ID) const {return _traversing_ORFs.at(colour_ID);};
    void clear_ORFs (const size_t& colour_ID) {_traversing_ORFs.erase(colour_ID);};

    private:
    // Bifrost data accessor methods
    inline uint8_t get() const { return da_id; };
    inline void set(const uint8_t id) { da_id = id; };
    // Bifrost data accessor ID
    uint8_t da_id = 0;

    // custom unitig IDs
    int _id;
    std::string _head_kmer;

    // codon arrays, initialise with two strands and 3 frames for each (6 reading frames total)
    // TODO remove 9-bit long full_codon as not used, place only first 3
    // for codons from 0th position - first 9 are forward strand, second 9 are reverse.
    // Of each 9, first 3 = 0th frame, second 3 = 1st frame, third 3 = 2nd frame
    std::bitset<18> _full_codon;
    std::bitset<18> _part_codon;
//    std::vector<std::vector<uint8_t>> _full_codon{std::vector<uint8_t>(3, 0), std::vector<uint8_t>(3, 0)};
//    std::vector<std::vector<uint8_t>> _part_codon{std::vector<uint8_t>(3, 0), std::vector<uint8_t>(3, 0)};

    // unitig colours
    boost::dynamic_bitset<> _unitig_full_colour;
    boost::dynamic_bitset<> _unitig_head_colour;
    boost::dynamic_bitset<> _unitig_tail_colour;
    bool _head_tail_colours_equal;

    // check head and tail colours equal
    void _check_head_tail_equal();

    // bool to determine if unitig is end of a contig/assembly
    bool _end_contig = false;

    // forward_stop presence/absence
    bool _forward_stop = false;
    bool _reverse_stop = false;

    // node neighbours. Neighbours map contains successors (true) and predecessors (false)
    // TODO remove this, fall back on getsuccessors for Bifrost
//    std::vector<NeighbourVector> _neighbours{NeighbourVector(), NeighbourVector()};

    // traversing ORFs (key is colour_ID, entry is traversing ORF)
    std::unordered_map<size_t, std::unordered_set<size_t>> _traversing_ORFs;

    // mapping of colour (key) to tuple of contig ID, mapping position within contig
    // and mapping position within node (start and length in k-mers) and boolean of strand (1 if same strand, 0 if reversed)
    // TODO remove, replace with mapping of all ORFs to FM-indexes if required for GFF
//    std::unordered_map<size_t, std::vector<std::tuple<size_t, size_t, size_t, size_t, bool>>> _contigcoords;
};

// tuple of GraphVector, a mapping of colours to component nodes, the number of colours and the size of the overlap
typedef std::tuple<NodeColourVector, std::vector<std::string>, size_t, int> GraphTuple;

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

//size_t get_id (const ColoredCDBG<MyUnitigMap>& ccdbg,
//               const Kmer& head_kmer);

// non-member functions for generating sequence from graph
std::string generate_sequence_nm(const std::vector<int>& nodelist,
                                 const std::vector<indexPair>& node_coords,
                                 const size_t& overlap,
                                 const ColoredCDBG<MyUnitigMap>& ccdbg,
                                 const std::vector<Kmer>& head_kmer_arr);

#endif //UNITIG_DICT_H