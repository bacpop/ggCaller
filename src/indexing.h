#ifndef INDEXING_H
#define INDEXING_H

#include "unitigDict.h"
#include "kseq.h"

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
                                   const int overlap,
                                   const bool reverse);

std::bitset<9> calculateFrame_binary (const std::vector<std::size_t>& index_list);

//std::bitset<3> switchFrame_binary (const std::bitset<3> binary_array, const int frame);

template <class T, class U, bool is_const>
boost::dynamic_bitset<> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                         const size_t nb_colours,
                                         const size_t position);

template<class T>
std::vector<std::pair<Kmer, bool>> get_neighbours (const T& neighbour_iterator);

template <class T, class U, bool is_const>
void analyse_unitigs_binary (ColoredCDBG<MyUnitigMap>& ccdbg,
                             UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
                             const std::vector<std::string>& codon_for,
                             const std::vector<std::string>& codon_rev,
                             const int& kmer,
                             const size_t& nb_colours);

void update_neighbour_index(ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::vector<Kmer>& head_kmer_arr);

void calculate_genome_paths(const std::vector<Kmer>& head_kmer_arr,
                            ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::string& fasta_file,
                            const int kmer,
                            const int colour_ID);

NodeColourVector index_graph(std::vector<Kmer>& head_kmer_arr,
                             ColoredCDBG<MyUnitigMap>& ccdbg,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             const int kmer,
                             const size_t nb_colours,
                             const bool is_ref,
                             const std::vector<std::string>& input_colours);

#endif //INDEXING_H
