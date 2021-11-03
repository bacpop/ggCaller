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
#include <mutex>
#include <math.h>

// openMP headers
#include <omp.h>

// robinhood hashing header
#include "robin_hood.h"

// sdsl header
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Eigen header
#include <Eigen/Sparse>

// edlib header
#include "edlib/edlib.h"

// bifrost header
#include <bifrost/ColoredCDBG.hpp>

// global variable declaration
namespace py = pybind11;

// UnitigDict typedefs
// Vector of neighbouring nodes by ID, orientation and map of stop codon frames
typedef std::vector<std::pair<int, std::vector<uint8_t>>> NeighbourVector;

// Eigen typedef
typedef Eigen::Triplet<double> ET;

//typedef sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<256>>,
//                    512,
//                    1024> fm_index_coll;
//csa_wt_int<sdsl::wt_huff<sdsl::rrr_vector<256>>, 512, 1024> fm_index_coll;
typedef sdsl::csa_wt<> fm_index_coll;

// general typedefs
// hasher for strings
typedef std::hash<std::string> hasher;
// mapping of each colour to component nodes in graph
typedef std::vector<std::vector<size_t>> NodeColourVector;
//a pair of start and end coordinates for an ORF across a node
typedef std::pair<size_t, size_t> indexPair;
// tuple holding ORF path ID, nodes traversed, node coordinates, coordinates in path, 5p and 3p coordinates
typedef std::tuple<std::vector<int>, std::vector<indexPair>> ORFCoords;
// tuple containing a vector of nodeIDs, a vector of start,stop and length coordinates, strand information, length of an ORF, TIS coordinate information, relative strand and population ID
typedef std::tuple<std::vector<int>, std::vector<indexPair>, size_t, std::vector<int>, std::vector<indexPair>, bool> ORFNodeVector;
// maps an ORFNodeVector sequence to its ID in order
typedef std::map<size_t, ORFNodeVector> ORFNodeMap;
// maps an map of ORFNodeVector sequence to its colour
typedef std::unordered_map<size_t, std::unordered_map<size_t, ORFNodeVector>> ColourORFMap;
// vector of ORF paths through graphs
typedef std::vector<ORFNodeVector> ORFVector;
// tuple for holding node information during traversal (1st = path index, 2nd = node id, 3rd = codon array, 4th = colour array, 5th = path length)
typedef std::tuple<size_t, int, uint8_t, sdsl::bit_vector, size_t> NodeTuple;
// stack for holding nodes during DFS traversal for ORF identification
typedef std::stack<NodeTuple> NodeStack;
// stack for holding node ids, previous ORF ID, the next expected node, colour array and path length during DFS traversal for ORF ordering
typedef std::stack<std::tuple<int, sdsl::bit_vector, size_t>> ORFStack;
// A vector of paths following a head node, which contain complete stop-stop paths (a vector of nodesID+orientation)
typedef std::vector<std::vector<int>> PathVector;
// mapping of node ID to a orientation for a specific strand, used in overlap analysis
typedef std::map<size_t, bool> NodeStrandMap;
// mapping of overlapping ORFs, detailed by ORFIDMap
typedef std::unordered_map<size_t, std::unordered_map<size_t, std::pair<char, size_t>>> ORFOverlapMap;
// vector that maps colour/ORF_ID to a new 1D index for fast searching, and maps homologous IDs in same vector
typedef std::vector<std::pair<size_t, size_t>> ORFMatrixVector;
// tuple containing grouping information for ORFs filtered by Balrog
typedef std::tuple<ORFMatrixVector, std::vector<std::unordered_set<size_t>>, std::vector<std::vector<std::pair<size_t, size_t>>>> ORFGroupTuple;
// map of ORFs to clusters, with centroid as key
typedef std::unordered_map<size_t, std::vector<size_t>> ORFClusterMap;
// tuple of ORF sequence, node list and node coordinates for orientation
typedef std::tuple<std::string, std::vector<int>, std::vector<std::vector<size_t>>> RefindTuple;

#endif //DEFINITIONS_H
