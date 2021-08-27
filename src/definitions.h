#ifndef DEFINITIONS_H
#define DEFINITIONS_H

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
typedef std::vector<std::pair<int, std::vector<uint8_t>>> NeighbourVector;

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

// general typedefs
// mapping of each colour to component nodes in graph
typedef std::vector<std::vector<size_t>> NodeColourVector;
//a pair of start and end coordinates for an ORF across a node
typedef std::pair<size_t, size_t> indexPair;
// tuple holding ORF path ID, nodes traversed, node coordinates, coordinates in path, 5p and 3p coordinates
typedef std::tuple<std::vector<int>, std::vector<indexPair>, indexPair> ORFCoords;
// tuple containing a vector of nodeIDs, a vector of start,stop and length coordinates, strand information, length of an ORF, TIS coordinate information, relative strand, vector of paths originated from and coordinates, and 5p and 3p coordinates
typedef std::tuple<std::vector<int>, std::vector<indexPair>, size_t, std::vector<int>, std::vector<indexPair>, bool, size_t, indexPair> ORFNodeVector;
// maps an ORF node sequence to its path through graph
typedef robin_hood::unordered_map<size_t, ORFNodeVector> ORFNodeMap;
// vector of ORF paths through graphs
typedef std::vector<ORFNodeVector> ORFVector;
// tuple for holding node information during traversal (1st = path index, 2nd = node id, 3rd = codon array, 4th = colour array, 5th = path length)
typedef std::tuple<size_t, int, uint8_t, sdsl::bit_vector, size_t> NodeTuple;
// stack for holding nodes during DFS traversal for ORF identification
typedef std::stack<NodeTuple> NodeStack;
// stack for holding node ids, previous ORF ID, the next expected node, colour array and path length during DFS traversal for ORF ordering
typedef std::stack<std::tuple<int, size_t, std::tuple<int, size_t, bool>, sdsl::bit_vector, size_t>> ORFStack;
// A vector of paths following a head node, which contain complete stop-stop paths (a vector of nodesID+orientation)
typedef std::vector<std::vector<int>> PathVector;
// A vector of all paths generated from recursive traversal
typedef std::vector<PathVector> AllPaths;
// mapping of node ID to a orientation for a specific strand, used in overlap analysis
typedef robin_hood::unordered_map<size_t, bool> NodeStrandMap;
// mapping of overlapping ORFs, detailed by ORFIDMap
typedef std::unordered_map<size_t, std::unordered_map<size_t, std::pair<char, size_t>>> ORFOverlapMap;
// mapping for each path with overlapping path id and how overlapping path is orientated relative to current path
typedef std::unordered_map<size_t, std::vector<std::vector<std::tuple<int, size_t, size_t>>>> PathOverlapMap;
// stack for holding next paths to be traversed


#endif //DEFINITIONS_H
