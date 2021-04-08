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
    std::pair<std::size_t, std::size_t> unitig_size;
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
typedef std::vector<std::unordered_set<size_t>> NodeColourVector;
// a tuple of UnitigVector, unitigs that contain stop codons in forward/reverse, and mappings of head-kmers to node IDs
typedef std::pair<UnitigVector, NodeColourVector> GraphPair;
// tuple of UnitigVector, a mapping of colours to component nodes, the number of colours and the size of the overlap
typedef std::tuple<UnitigVector, NodeColourVector, std::vector<std::string>, size_t, int> GraphTuple;
////a vector of start,stop and length coordinates and strand information for an ORF
//typedef std::tuple<size_t, size_t, size_t> indexTriplet;
//// tuple containing a vector of nodeIDs, a vector of start,stop and length coordinates, strand information, length of an ORF and TIS coordinate information
//typedef std::tuple<std::vector<int>, std::vector<indexTriplet>, size_t, std::vector<int>, std::vector<indexTriplet>> ORFNodeVector;
//a pair of start and end coordinates for an ORF across a node
typedef std::pair<size_t, size_t> indexPair;
// tuple containing a vector of nodeIDs, a vector of start,stop and length coordinates, strand information, length of an ORF and TIS coordinate information
typedef std::tuple<std::vector<int>, std::vector<indexPair>, size_t, std::vector<int>, std::vector<indexPair>> ORFNodeVector;
// maps an ORF node sequence to its path through graph
typedef robin_hood::unordered_map<std::string, ORFNodeVector> ORFNodeMap;
// vector of ORF paths throughg graphs
typedef std::vector<ORFNodeVector> ORFVector;
// A vector of paths following a head node, which containg complete stop-stop paths (a vector of nodesID+orientation)
typedef std::vector<std::vector<int>> PathVector;
// A vector of all paths generated from recursive traversal
typedef std::vector<PathVector> AllPaths;
// mapping of node ID to a orientation for a specific strand, used in overlap analysis
typedef robin_hood::unordered_map<size_t, bool> NodeStrandMap;
// mapping of overlapping ORFs, detailed by ORFIDMap
typedef std::unordered_map<size_t, std::unordered_map<size_t, std::pair<char, size_t>>> ORFOverlapMap;

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

std::vector<bool> negate_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

std::vector<bool> add_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

template<class T>
std::vector<std::pair<std::string, bool>> get_neighbours (const T& neighbour_iterator);

template <class T, class U, bool is_const>
unitigDict analyse_unitigs_binary (const ColoredCDBG<>& ccdbg,
                                   UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
                                   const std::vector<std::string>& codon_for,
                                   const std::vector<std::string>& codon_rev,
                                   const int& kmer,
                                   const size_t& nb_colours);

void update_neighbour_index(UnitigVector& graph_vector,
                            robin_hood::unordered_map<std::string, size_t> head_kmer_map);

GraphPair index_graph(const ColoredCDBG<>& ccdbg,
                       const std::vector<std::string>& stop_codons_for,
                       const std::vector<std::string>& stop_codons_rev,
                       const int kmer,
                       const size_t nb_colours);

// traversal.cpp
PathVector recur_nodes_binary (const UnitigVector& graph_vector,
                               const std::vector<int>& head_kmer_list,
                               const uint8_t& codon_arr,
                               const size_t& colour_ID,
                               const std::unordered_set<int>& kmer_set,
                               const size_t& length,
                               const size_t& length_max,
                               const bool& repeat);

AllPaths traverse_graph(const UnitigVector& graph_vector,
                         const size_t& colour_ID,
                         const std::unordered_set<size_t>& node_ids,
                         const bool repeat,
                         const size_t max_path_length);

// match_strings
fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx);

int seq_search(const seqan3::dna5_vector& query,
               const fm_index_coll& ref_idx);

bool check_colours(const std::string& path_sequence,
                   const fm_index_coll& fm_idx);

// call_ORFs
ORFNodeMap generate_ORFs(const UnitigVector& graph_vector,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const std::vector<int>& unitig_path,
                         const int& overlap,
                         const size_t min_len,
                         const bool is_ref,
                         const fm_index_coll& fm_idx);

std::tuple<std::string, std::vector<int>, std::vector<indexPair>> calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                                                                                    const std::vector<int>& nodelist,
                                                                                    const std::vector<std::vector<size_t>>& node_ranges,
                                                                                    const int& overlap);

std::pair<ORFVector, NodeStrandMap> call_ORFs(const AllPaths& all_paths,
                                             const UnitigVector& graph_vector,
                                             const std::vector<std::string>& stop_codons_for,
                                             const std::vector<std::string>& start_codons_for,
                                             const int overlap,
                                             const size_t min_ORF_length,
                                             const bool is_ref,
                                             const fm_index_coll& fm_idx);

ORFVector sort_ORF_indexes(ORFNodeMap& ORF_node_map);

NodeStrandMap calculate_pos_strand(const ORFNodeMap& ORF_node_map);

// gene_overlap.cpp
ORFOverlapMap calculate_overlaps(const UnitigVector& graph_vector,
                                 const std::pair<ORFVector, NodeStrandMap>& ORF_pair,
                                 const int DBG_overlap,
                                 const size_t max_overlap);

// ggCaller_bindings
GraphTuple py_index_graph_exists(const std::string& graphfile,
                                 const std::string& coloursfile,
                                 const std::vector<std::string>& stop_codons_for,
                                 const std::vector<std::string>& stop_codons_rev,
                                 size_t num_threads,
                                 const bool is_ref);

GraphTuple py_index_graph_build(const std::string& infile1,
                                const int kmer,
                                const std::vector<std::string>& stop_codons_for,
                                const std::vector<std::string>& stop_codons_rev,
                                size_t num_threads,
                                bool is_ref,
                                const bool write_graph,
                                const std::string& infile2);

std::pair<ORFOverlapMap, ORFVector> py_calculate_ORFs (const UnitigVector& graph_vector,
                                                     const size_t& colour_ID,
                                                     const std::unordered_set<size_t>& node_ids,
                                                     const bool& repeat,
                                                     const size_t& overlap,
                                                     const size_t& max_path_length,
                                                     bool& is_ref,
                                                     const bool& no_filter,
                                                     const std::vector<std::string>& stop_codons_for,
                                                     const std::vector<std::string>& start_codons_for,
                                                     const size_t min_ORF_length,
                                                     const size_t max_overlap,
                                                     const bool write_idx,
                                                     const std::string& FM_fasta_file);

#endif //BIFROST_API_GGCALLER_H
