#ifndef INDEXING_H
#define INDEXING_H

#include "unitigDict.h"

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
boost::dynamic_bitset<> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
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

void update_neighbour_index(GraphVector& graph_vector,
                            robin_hood::unordered_map<std::string, size_t> head_kmer_map);

GraphPair index_graph(const ColoredCDBG<>& ccdbg,
                       const std::vector<std::string>& stop_codons_for,
                       const std::vector<std::string>& stop_codons_rev,
                       const int kmer,
                       const size_t nb_colours);

#endif //INDEXING_H
