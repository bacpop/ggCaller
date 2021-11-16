#ifndef UNITIG_DICT_H
#define UNITIG_DICT_H

#include "definitions.h"

// UnitigDict typedefs and functions
// Vector of neighbouring nodes by ID, orientation and map of stop codon frames
typedef std::vector<std::pair<int, std::vector<uint8_t>>> NeighbourVector;

// unitigDict class declaration
class unitigDict {
    public:

    // add and return untig kmers and ids
    void add_head(const std::string& head) {_head_kmer = head;};
    const std::string& head_kmer() const {return _head_kmer;};
    size_t id;
//    void set_id (size_t id) {_unitig_id = id};
//    size_t id () {return _unitig_id};

    // add codon information
    void add_codon (const bool& full, const bool& forward, const int& frame, const uint8_t& array);
    //void add_codon (const bool& full, const bool& forward, const int& frame, uint8_t& array);
    void set_forward_stop (bool choice) {_forward_stop = choice;};
    void set_reverse_stop (bool choice) {_reverse_stop = choice;};
    bool forward_stop () const {return _forward_stop;};
    bool reverse_stop () const {return _reverse_stop;};

    // get codon information
    uint8_t get_codon_arr (bool full, bool forward, bool frame) const;
    std::vector<uint8_t> get_codon_dict (bool full, bool forward) const;

    // add size information
    void add_size(const size_t& full_len, const size_t& part_len);
    void add_size(size_t& full_len, size_t& part_len);
    std::pair<std::size_t, std::size_t> size() const {return _unitig_size;};

    // add and return unitig seqs
    void add_seq(const std::string& seq) {_unitig_seq = seq;};
    std::string seq() const {return _unitig_seq;};

    // add and return unitig colours
    void add_head_colour(sdsl::bit_vector colours) {_unitig_head_colour = colours;};
    sdsl::bit_vector head_colour() const {return _unitig_head_colour;};
    void add_full_colour(sdsl::bit_vector colours) {_unitig_full_colour = colours;};
    sdsl::bit_vector full_colour() const {return _unitig_full_colour;};
    void add_tail_colour(sdsl::bit_vector colours) {_unitig_tail_colour = colours; _check_head_tail_equal();};
    sdsl::bit_vector tail_colour() const {return _unitig_tail_colour;};
    bool head_tail_colours_equal() const {return _head_tail_colours_equal;};

    // access end_contig
    void set_end_contig (bool choice) {_end_contig = choice;};
    bool end_contig() const {return _end_contig;};

    //assign and return neighbours
    void set_succ (std::vector<std::pair<std::string, bool>> vect) {_succ_heads = vect;};
    void clear_succ() {_succ_heads.clear(); _succ_heads.shrink_to_fit();}
    const std::vector<std::pair<std::string, bool>> & get_succs () {return _succ_heads;};
    void set_pred (std::vector<std::pair<std::string, bool>> vect) {_pred_heads = vect;};
    void clear_pred() {_pred_heads.clear(); _pred_heads.shrink_to_fit();}
    const std::vector<std::pair<std::string, bool>> & get_preds () const {return _pred_heads;};
    void add_neighbour (bool strand, std::pair<int, std::vector<uint8_t>> neighbour) {_neighbours[strand].push_back(std::move(neighbour));};

    // get neighbouring nodes
    const NeighbourVector & get_neighbours (bool strand) const {return _neighbours[strand];};

    // assign traversing ORFs
    void set_ORFs (const size_t& colour_ID, const size_t& ORF_ID) {_traversing_ORFs[colour_ID].insert(ORF_ID);};
    bool ORFs_empty (const size_t& colour_ID) const {return (_traversing_ORFs.find(colour_ID) == _traversing_ORFs.end());};
    const std::unordered_set<size_t> & get_ORFs(const size_t& colour_ID) const {return _traversing_ORFs.at(colour_ID);};

    private:
    std::string _head_kmer;

    // codon arrays, initialise with two strands and 3 frames for each (6 reading frames total)
    std::vector<std::vector<uint8_t>> _full_codon{std::vector<uint8_t>(3, 0), std::vector<uint8_t>(3, 0)};
    std::vector<std::vector<uint8_t>> _part_codon{std::vector<uint8_t>(3, 0), std::vector<uint8_t>(3, 0)};

    // unitig properties
    std::pair<std::size_t, std::size_t> _unitig_size;

    // unitig colours
    sdsl::bit_vector _unitig_full_colour;
    sdsl::bit_vector _unitig_head_colour;
    sdsl::bit_vector _unitig_tail_colour;
    bool _head_tail_colours_equal;

    // check head and tail colours equal
    void _check_head_tail_equal();

    // unitig sequence
    std::string _unitig_seq;

    // bool to determine if unitig is end of a contig/assembly
    bool _end_contig = false;

    bool _forward_stop_defined = false;
    bool _reverse_stop_defined = false;

    // forward_stop presence/absence
    bool _forward_stop = false;
    bool _reverse_stop = false;

    // node neighbours. Neighbours map contains successors (true) and predecessors (false)
    std::vector<std::pair<std::string, bool>> _succ_heads;
    std::vector<std::pair<std::string, bool>> _pred_heads;
    std::vector<NeighbourVector> _neighbours{NeighbourVector(), NeighbourVector()};

    // traversing ORFs (key is colour_ID, entry is traversing ORF)
    std::unordered_map<size_t, std::unordered_set<size_t>> _traversing_ORFs;
};

// unitigDict typedefs
// mapping of unitig IDs (size_t) to unitigDict class for each unitig
typedef std::vector<unitigDict> GraphVector;
// a tuple of GraphVector, unitigs that contain stop codons in forward/reverse, and mappings of head-kmers to node IDs
typedef std::pair<GraphVector, NodeColourVector> GraphPair;
// tuple of GraphVector, a mapping of colours to component nodes, the number of colours and the size of the overlap
typedef std::tuple<NodeColourVector, std::vector<std::string>, size_t, int> GraphTuple;

#endif //UNITIG_DICT_H
