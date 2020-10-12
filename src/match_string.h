//
// Created by sth19 on 25/08/2020.
//

// C/C++/C++11/C++17 headers
#include <cinttypes>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iterator>
#include <vector>
#include <functional>
#include <cassert>
#include <experimental/filesystem>
#include <unordered_map>

// openMP headers
#include <omp.h>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// seqan3 headers
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

namespace py = pybind11;

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



// Function headers
// match_strings.cpp
std::unordered_map<std::string, std::vector<std::string>> call_strings(const std::vector<std::string>& assembly_list,
                                                                       std::unordered_map<std::string, std::vector<std::string>>& query_list,
                                                                       const bool& write_idx,
                                                                       const size_t& num_threads);

std::vector<fm_index_coll> index_fasta(const std::vector<std::string>& fasta_files,
                                       const size_t start,
                                       const size_t end,
                                       const bool& write_idx);

int seq_search(const seqan3::dna5_vector& query,
               const std::vector<fm_index_coll>& seq_idx,
               const size_t start,
               const size_t end);

// match_bindings.cpp
std::unordered_map<std::string, std::vector<std::string>> py_call_strings(const std::vector<std::string>& assembly_list,
                                                                          std::unordered_map<std::string, std::vector<std::string>>& query_list,
                                                                          const bool& write_idx,
                                                                          size_t& num_threads);