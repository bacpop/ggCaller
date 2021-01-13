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
const vector<bool> empty_codon_arr(3, 0);
namespace py = pybind11;

// class declaration
class unitigDict {
    public:

    // add untig ids
    void add_head(const std::string& head) {head_kmer = head;};
    void add_id(int id) {unitig_id = id;};

    // add codon information
    void add_codon (const bool& full, const bool& forward, const int& frame, const uint8_t& array);
    void add_codon (const bool& full, const bool& forward, const int& frame, uint8_t& array);

    // add size information
    void add_size(const size_t& full_len, const size_t& part_len);
    void add_size(size_t& full_len, size_t& part_len);

    // add colour information
    void add_colour(const std::vector<bool>& array);
    void add_colour(std::vector<bool>& array);

//    // reset unitigDict
//    void reset();

    // access private members
//    const int& returnID() {return unitig_id;};
//    const std::string& returnHead() {return head_kmer;};
//    const robin_hood::unordered_map<bool, robin_hood::unordered_map<int, std::vector<bool>>>& returnFullcodon() {return full_codon;};
//    const robin_hood::unordered_map<bool, robin_hood::unordered_map<int, std::vector<bool>>>& returnPartcodon() {return part_codon;};
//    const std::pair<size_t, std::size_t>& returnSize() {return unitig_size;};
//    const std::vector<bool>& returnColours() {return unitig_colour;};
//    const bool& returnForstop() {return forward_stop;};
//    const bool& returnRevstop() {return reverse_stop;};

    // to add: add unitig sequence, predecessor/successor ID and orientation for fast access with robin-hood maps

    size_t unitig_id;
    std::string head_kmer;
    robin_hood::unordered_map<bool, robin_hood::unordered_map<int, uint8_t>> full_codon;
    robin_hood::unordered_map<bool, robin_hood::unordered_map<int, uint8_t>> part_codon;
    std::pair<size_t, std::size_t> unitig_size;
    std::vector<bool> unitig_colour;
    bool forward_stop = false;
    bool reverse_stop = false;

    private:
    bool forward_stop_defined = false;
    bool reverse_stop_defined = false;
};

// ggCaller typedefs
typedef std::tuple<size_t, size_t, size_t> indexTriplet;
typedef robin_hood::unordered_map<std::string, unitigDict> unitigMap;
typedef std::pair<std::vector<std::string>, std::vector<indexTriplet>> ORFNodeVector;
typedef robin_hood::unordered_map<std::string, ORFNodeVector> ORFNodeMap;
typedef robin_hood::unordered_map<std::string, ORFNodeMap> StrandORFNodeMap;
typedef robin_hood::unordered_map<std::string, std::vector<bool>> SeqORFMap;
typedef robin_hood::unordered_map<std::string, SeqORFMap> StrandSeqORFMap;


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
                          const bool& is_ref,
                          const int kmer,
                          const int threads,
                          const bool verb,
                          const bool& write_graph,
                          const std::string& output_prefix);

template <class T, class U, bool is_const>
std::vector<bool> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                   const size_t& nb_colours);

template <class T, class U, bool is_const>
unitigDict analyse_unitigs_binary (const ColoredCDBG<>& ccdbg,
                                   const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
                                   const std::vector<std::string>& codon_for,
                                   const std::vector<std::string>& codon_rev,
                                   const int& kmer,
                                   const size_t& nb_colours);

std::tuple<unitigMap, std::vector<std::string>, std::vector<std::string>, robin_hood::unordered_map<size_t, std::string>> index_graph(const ColoredCDBG<>& ccdbg,
                                                                                                                                      const std::vector<std::string>& stop_codons_for,
                                                                                                                                      const std::vector<std::string>& stop_codons_rev,
                                                                                                                                      const int& kmer,
                                                                                                                                      const size_t& nb_colours);
// traversal.cpp
std::vector<bool> compare_codon_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

std::vector<bool> negate_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

std::vector<bool> add_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

std::vector<std::pair<std::vector<std::pair<std::string, bool>>, std::vector<bool>>> recur_nodes_binary (const ColoredCDBG<>& ccdbg,
                                                                                                         const robin_hood::unordered_map<std::string, unitigDict>& graph_map,
                                                                                                         const robin_hood::unordered_map<std::string, std::vector<std::pair<std::vector<std::pair<std::string, bool>>, std::vector<bool>>>>& previous_paths,
                                                                                                         const std::vector<std::pair<std::string, bool>> head_kmer_list,
                                                                                                         const uint8_t& codon_arr,
                                                                                                         const std::vector<bool>& colour_arr,
                                                                                                         const std::set<std::pair<std::string, bool>> kmer_set,
                                                                                                         const size_t length,
                                                                                                         const bool& forward,
                                                                                                         const size_t& length_max,
                                                                                                         const bool repeat,
                                                                                                         const vector<bool>& empty_colour_arr);

std::tuple<robin_hood::unordered_map<std::string, std::vector<std::pair<std::vector<std::pair<std::string, bool>>, std::vector<bool>>>>, std::vector<std::string>> traverse_graph(const ColoredCDBG<>& ccdbg,
                                                                                                                                                                                  const std::tuple<robin_hood::unordered_map<std::string, unitigDict>, std::vector<std::string>, std::vector<std::string>, robin_hood::unordered_map<size_t, std::string>>& graph_tuple,
                                                                                                                                                                                  const bool& repeat,
                                                                                                                                                                                  const vector<bool>& empty_colour_arr,
                                                                                                                                                                                  const size_t& max_path_length);

// match_strings

fm_index_coll index_fasta(const std::string& fasta_file,
                          const bool& write_idx);

int seq_search(const seqan3::dna5_vector& query,
               const fm_index_coll& ref_idx,
               const std::string& strand);

void call_strings(SeqORFMap& query_list,
                  const std::string& strand,
                  ORFNodeMap& ORF_node_paths,
                  const std::vector<std::string>& assembly_list,
                  const bool& write_idx);

// call_ORFs
ORFNodeMap generate_ORFs(const ColoredCDBG<>& ccdbg,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const std::vector<std::pair<std::string, bool>>& unitig_path,
                         const int& overlap,
                         const size_t& min_len);

std::tuple<StrandSeqORFMap, StrandORFNodeMap> call_ORFs(const ColoredCDBG<>& ccdbg,
                                                        const std::tuple<robin_hood::unordered_map<std::string, std::vector<std::pair<std::vector<std::pair<std::string, bool>>, std::vector<bool>>>>, std::vector<std::string>>& path_tuple,
                                                        const std::vector<std::string>& stop_codons_for,
                                                        const std::vector<std::string>& start_codons_for,
                                                        const int& overlap,
                                                        const size_t& min_ORF_length);

std::unordered_map<std::string, std::vector<std::string>> filter_artificial_ORFS(StrandSeqORFMap& all_ORFs,
                                                                                 StrandORFNodeMap& ORF_node_paths,
                                                                                 const std::vector<std::string>& fasta_files,
                                                                                 const bool write_index);

void write_to_file (const std::string& outfile_name,
                    const StrandSeqORFMap& all_ORFs);

// gene_overlap.cpp
void calculate_overlaps(const unitigMap& unitig_map,
                        const StrandORFNodeMap& ORF_node_paths,
                        const std::unordered_map<std::string, std::vector<std::string>> ORF_colours_map,
                        const int& node_overlap,
                        const size_t& max_overlap);

// ggCaller_bindings
int py_ggCaller_graphexists (const std::string& graphfile,
                             const std::string& coloursfile,
                             const std::string& outfile,
                             const std::vector<std::string>& start_codons,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             size_t num_threads,
                             const bool is_ref,
                             const bool write_idx,
                             const bool repeat,
                             const size_t& max_path_length,
                             const size_t& min_ORF_length);

int py_ggCaller_graphbuild (const std::string& infile1,
                            const int& kmer,
                            const std::string& outfile,
                            const std::vector<std::string>& start_codons,
                            const std::vector<std::string>& stop_codons_for,
                            const std::vector<std::string>& stop_codons_rev,
                            size_t num_threads,
                            bool is_ref,
                            const bool write_idx,
                            const bool repeat,
                            const bool write_graph,
                            const size_t& max_path_length,
                            const size_t& min_ORF_length,
                            const std::string& infile2);

#endif //BIFROST_API_GGCALLER_H
