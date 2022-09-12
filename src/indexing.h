#ifndef INDEXING_H
#define INDEXING_H

#include "unitigDict.h"
#include "kseq.h"
#include "translation.h"

ColoredCDBG<MyUnitigMap> buildGraph (const std::string& infile_1,
                                     const std::string& infile_2,
                                     const bool is_ref,
                                     const int kmer,
                                     const int threads,
                                     const bool verb,
                                     const bool write_graph,
                                     const std::string& output_prefix);

std::vector<std::size_t> findIndex(const std::string& seq,
                                   const std::string& subseq,
                                   const int start_index,
                                   const bool reverse);

std::bitset<3> calculateFrame_binary_full (const std::vector<std::size_t>& index_list);

std::bitset<9> calculateFrame_binary_part (const std::vector<std::size_t>& index_list);

template <class T, class U, bool is_const>
boost::dynamic_bitset<> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                         const size_t nb_colours,
                                         const size_t position);

template<class T>
std::vector<std::pair<Kmer, bool>> get_neighbours (const T& neighbour_iterator);

template <class T, class U, bool is_const>
void analyse_unitigs_binary (ColoredCDBG<MyUnitigMap>& ccdbg,
                             UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
                             size_t& num_stops,
                             size_t& num_codons,
                             const std::vector<std::string>& stop_codon_for,
                             const std::vector<std::string>& stop_codon_rev,
                             const std::vector<std::string>& start_codon_for,
                             const std::vector<std::string>& start_codon_rev,
                             const int& kmer,
                             const size_t& nb_colours,
                             tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_set<int>>& start_freq_set);

void calculate_genome_paths(const std::vector<Kmer>& head_kmer_arr,
                            ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::string& fasta_file,
                            const int kmer,
                            const int colour_ID,
                            const size_t nb_colours);

NodeColourVector index_graph(std::vector<Kmer>& head_kmer_arr,
                             ColoredCDBG<MyUnitigMap>& ccdbg,
                             float& stop_codon_freq,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             const std::vector<std::string>& start_codons_for,
                             const std::vector<std::string>& start_codons_rev,
                             const int kmer,
                             const size_t nb_colours,
                             const std::vector<std::string>& input_colours,
                             const boost::dynamic_bitset<>& ref_set,
                             tbb::concurrent_unordered_map<size_t, size_t>& start_freq);

#endif //INDEXING_H
