#ifndef BIFROST_API_GGCALLER_H
#define BIFROST_API_GGCALLER_H

// C/C++/C++11/C++17 headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include <algorithm>
#include <set>
#include <cinttypes>
#include <cstdlib>
#include <iterator>
#include <functional>
#include <cassert>
#include <experimental/filesystem>
#include <mutex>

// openMP headers
#include <omp.h>

// robinhood hashing header
#include "robin_hood.h"

// seqan3 headers
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Eigen header
#include "Eigen/Sparse"

// bifrost header
#include <bifrost/ColoredCDBG.hpp>

// global variable declaration
namespace py = pybind11;

// UnitigDict typedefs
// Vector of neighbouring nodes by ID, orientation and map of stop codon frames
typedef std::vector<std::pair<int, std::unordered_map<int, uint8_t>>> NeighbourVector;

// Eigen typedef
typedef Eigen::Triplet<double> ET;

// fmindex typedef
using cust_sdsl_wt_index_type = sdsl::csa_wt<sdsl::wt_blcd<sdsl::bit_vector,
        sdsl::rank_support_v<>,
        sdsl::select_support_scan<>,
        sdsl::select_support_scan<0>>,
        16,
        10000000,
        sdsl::sa_order_sa_sampling<>,
        sdsl::isa_sampling<>,
        sdsl::plain_byte_alphabet>;
typedef seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection, cust_sdsl_wt_index_type> fm_index_coll;
using seqan3::operator""_dna5;

// class declaration
class unitigDict {
    public:

    // add untig ids
    void add_head(const std::string& head) {head_kmer = head;};

    // add codon information
    void add_codon (const bool& full, const bool& forward, const int& frame, const uint8_t& array);
    void add_codon (const bool& full, const bool& forward, const int& frame, uint8_t& array);

    // add size information
    void add_size(const size_t& full_len, const size_t& part_len);
    void add_size(size_t& full_len, size_t& part_len);

    // Public class items
    size_t unitig_id;
    std::string head_kmer;

    // codon arrays
    std::unordered_map<bool, std::unordered_map<int, uint8_t>> full_codon;
    std::unordered_map<bool, std::unordered_map<int, uint8_t>> part_codon;

    // unitig properties
    std::pair<size_t, std::size_t> unitig_size;
    std::vector<bool> unitig_full_colour;
    std::vector<bool> unitig_head_colour;
    std::vector<bool> unitig_tail_colour;
    bool head_tail_colours_equal = true;

    // unitig sequence
    std::string unitig_seq;

    // bool to determine if unitig is end of a contig/assembly
    bool end_contig = false;

    // node neighbours. Neighbours map contains successors (true) and predecessors (false)
    std::vector<std::pair<std::string, bool>> succ_heads;
    std::vector<std::pair<std::string, bool>> pred_heads;
    std::unordered_map<bool, NeighbourVector> neighbours;

    // forward_stop presence/absence
    bool forward_stop = false;
    bool reverse_stop = false;

    private:
    bool forward_stop_defined = false;
    bool reverse_stop_defined = false;
};

// ggCaller typedefs
// mapping of unitig IDs (size_t) to unitigDict class for each unitig
typedef std::vector<unitigDict> UnitigVector;
// mapping of each colour to component nodes in graph
typedef std::vector<std::vector<size_t>> NodeColourVector;
// a tuple of unitigMap, unitigs that contain stop codons in forward/reverse, and mappings of head-kmers to node IDs
typedef std::pair<UnitigVector, NodeColourMap> GraphPair;
//a vector of start,stop and length coordinates and strand information for an ORF
typedef std::tuple<size_t, size_t, size_t> indexTriplet;
// tuple containing a vector of nodeIDs, a vector of start,stop and length coordinates, strand information, length of an ORF and TIS coordinate information
typedef std::tuple<std::vector<int>, std::vector<indexTriplet>, size_t, std::vector<int>, std::vector<indexTriplet>> ORFNodeVector;
// maps an ORF node sequence to its colours and path through graph
typedef robin_hood::unordered_map<std::string, std::pair<std::vector<bool>, ORFNodeVector>> ORFNodeMap;
// maps individual colour ids to ORF ids in ORFIDMap
typedef robin_hood::unordered_map<size_t, std::vector<size_t>> ORFColoursMap;
// maps a unique ID to a path through graph
typedef robin_hood::unordered_map<size_t, ORFNodeVector> ORFIDMap;
// A vector of paths following a head node, which containg complete stop-stop paths (pair of a vector of nodesID+orientation, and colours vector)
typedef std::vector<std::pair<std::vector<int>, std::vector<bool>>> PathVector;
// Mapping of header kmer ID to PathVector
typedef robin_hood::unordered_map<int, PathVector> PathMap;
// pairing of PathMap and a vector of head-kmer IDs (hashing for PathMap) for parrelisation
typedef std::pair<PathMap, std::vector<int>> PathPair;
// mapping of node ID to a orientation for a specific strand, used in overlap analysis
typedef robin_hood::unordered_map<size_t, bool> NodeStrandMap;
// mapping of colour ID to a NodeStrandMap
typedef robin_hood::unordered_map<size_t, NodeStrandMap> ColourNodeStrandMap;
// mapping colour to each colour to map of overlapping ORFs, detailed by ORFIDMap
typedef std::unordered_map<size_t, std::unordered_map<size_t, std::unordered_map<size_t, std::pair<char, size_t>>>> ORFOverlapMap;
// maps individual colour ids to ORF ids in ORFIDMap, output for pybind
typedef std::unordered_map<size_t, std::vector<size_t>> PyORFColoursMap;
// maps a unique ID to an ORF sequence and path through graph. output for pybind
typedef std::unordered_map<size_t, ORFNodeVector> PyORFIDMap;
// mapping of unitig IDs (size_t) to unitigDict class for each unitig
typedef std::unordered_map<size_t, unitigDict> PyUnitigMap;

// function headers
// indexing
std::vector<std::size_t> findIndex(const std::string& seq,
                                   const std::string& subseq,
                                   const int start_index,
                                   const int overlap,
                                   const bool reverse);

uint8_t calculateFrame_binary (const std::vector<std::size_t>& index_list);

uint8_t switchFrame_binary (const uint8_t& binary_array, const int& frame);

ColoredCDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool is_ref,
                          const int kmer,
                          const int threads,
                          const bool verb,
                          const bool write_graph,
                          const std::string& output_prefix);

template <class T, class U, bool is_const>
std::vector<bool> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                   const size_t& nb_colours,
                                   const size_t position);

template<class T>
std::vector<std::pair<std::string, bool>> get_neighbours (const T& neighbour_iterator);

template <class T, class U, bool is_const>
unitigDict analyse_unitigs_binary (const ColoredCDBG<>& ccdbg,
                                   UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
                                   const std::vector<std::string>& codon_for,
                                   const std::vector<std::string>& codon_rev,
                                   const int& kmer,
                                   const size_t& nb_colours);

GraphTuple index_graph(const ColoredCDBG<>& ccdbg,
                       const std::vector<std::string>& stop_codons_for,
                       const std::vector<std::string>& stop_codons_rev,
                       const int kmer,
                       const size_t nb_colours);

// traversal.cpp
std::vector<bool> compare_codon_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

std::vector<bool> negate_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

std::vector<bool> add_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

PathVector recur_nodes_binary (const unitigMap& graph_map,
                               const std::vector<int>& head_kmer_list,
                               const uint8_t& codon_arr,
                               std::vector<bool> colour_arr,
                               const std::unordered_set<int>& kmer_set,
                               const size_t& length,
                               const size_t& length_max,
                               const bool& repeat,
                               const vector<bool>& empty_colour_arr);

PathPair traverse_graph(const GraphTuple& graph_tuple,
                         const bool repeat,
                         const vector<bool>& empty_colour_arr,
                         const size_t max_path_length);

// match_strings

fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx);

int seq_search(const seqan3::dna5_vector& query,
               const fm_index_coll& ref_idx);

void call_strings(ORFNodeMap& ORF_node_paths,
                  const std::vector<std::string>& assembly_list,
                  const bool& write_idx);

std::vector<fm_index_coll> generate_fmindex(const std::vector<std::string>& assembly_list,
                                            const bool& write_idx);

void check_colours(std::vector<bool>& path_colours,
                   const std::string& path_sequence,
                   const std::vector<fm_index_coll>& seq_idx,
                   const size_t& nb_colours);

// call_ORFs
ORFNodeMap generate_ORFs(const unitigMap& graph_map,
                         const std::vector<std::string>& stop_codons,
                         const std::string& start_codon,
                         const std::vector<std::pair<size_t, bool>>& unitig_path,
                         const int overlap,
                         const size_t min_len,
                         const bool is_ref,
                         const std::vector<fm_index_coll>& seq_idx,
                         std::vector<bool>& path_colours,
                         const size_t& nb_colours);

std::tuple<std::string, std::vector<int>, std::vector<indexTriplet>> calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                                                                                      const std::vector<int>& nodelist,
                                                                                      const std::vector<std::vector<size_t>>& node_ranges,
                                                                                      const int& overlap);

std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>, ColourNodeStrandMap> call_ORFs(const PathPair& path_pair,
                                                                                             const unitigMap& graph_map,
                                                                                             const std::vector<std::string>& stop_codons_for,
                                                                                             const std::vector<std::string>& start_codons_for,
                                                                                             const int overlap,
                                                                                             const size_t min_ORF_length,
                                                                                             const bool is_ref,
                                                                                             const std::vector<fm_index_coll>& seq_idx,
                                                                                             const size_t& nb_colours);

std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>> sort_ORF_colours(ORFNodeMap& ORF_node_paths);

//std::tuple<ORFColoursMap, ORFIDMap, std::vector<std::size_t>> filter_artificial_ORFS(ORFNodeMap& ORF_node_paths,
//                                                                                     const std::vector<std::string>& fasta_files,
//                                                                                     const bool write_index);

ColourNodeStrandMap calculate_pos_strand(const ORFNodeMap& ORF_node_paths);

// gene_overlap.cpp
ORFOverlapMap calculate_overlaps(const unitigMap& unitig_map,
                                 const std::tuple<ORFColoursMap, ORFIDMap, std::vector<size_t>, ColourNodeStrandMap>& ORF_colours_tuple,
                                 const int DBG_overlap,
                                 const size_t max_overlap);

// ggCaller_bindings
std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyUnitigMap, size_t, size_t> py_ggCaller_graphexists (const std::string& graphfile,
                                                                                                            const std::string& coloursfile,
                                                                                                           const std::vector<std::string>& start_codons,
                                                                                                           const std::vector<std::string>& stop_codons_for,
                                                                                                           const std::vector<std::string>& stop_codons_rev,
                                                                                                           size_t num_threads,
                                                                                                           const bool is_ref,
                                                                                                           const bool write_idx,
                                                                                                           const bool repeat,
                                                                                                           const bool no_filter,
                                                                                                           const size_t max_path_length,
                                                                                                           const size_t min_ORF_length,
                                                                                                           const size_t max_ORF_overlap);

std::tuple<ORFOverlapMap, PyORFColoursMap, PyORFIDMap, PyUnitigMap, size_t, size_t> py_ggCaller_graphbuild (const std::string& infile1,
                                                                                                          const int kmer,
                                                                                                          const std::vector<std::string>& start_codons,
                                                                                                          const std::vector<std::string>& stop_codons_for,
                                                                                                          const std::vector<std::string>& stop_codons_rev,
                                                                                                          size_t num_threads,
                                                                                                          bool is_ref,
                                                                                                          const bool write_idx,
                                                                                                          const bool repeat,
                                                                                                          const bool write_graph,
                                                                                                          const bool no_filter,
                                                                                                          const size_t max_path_length,
                                                                                                          const size_t min_ORF_length,
                                                                                                          const size_t max_ORF_overlap,
                                                                                                          const std::string& infile2);

#endif //BIFROST_API_GGCALLER_H
