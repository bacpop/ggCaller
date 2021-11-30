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

// boost header
#include <boost/dynamic_bitset.hpp>

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

typedef sdsl::csa_wt<> fm_index_coll;

// general typedefs
// hasher for strings
typedef std::hash<std::string> hasher;
// mapping of each colour to component nodes in graph
typedef std::vector<std::vector<size_t>> NodeColourVector;
// vector of contig mappings for a node in the graph
typedef std::vector<std::pair<size_t, std::tuple<size_t, size_t, size_t, bool>>> NodeContigMapping;
//a pair of start and end coordinates for an ORF across a node
typedef std::pair<size_t, size_t> indexPair;
// pair that describes the contig locations of an ORF, 1-indexed for contig id (first) and locations within contig (second)
typedef std::pair<size_t, std::pair<size_t, size_t>> ContigLoc;
// tuple holding ORF path ID, nodes traversed, node coordinates, coordinates in path, 5p and 3p coordinates
typedef std::tuple<std::vector<int>, std::vector<indexPair>> ORFCoords;
// tuple containing a vector of nodeIDs, a vector of start,stop and length coordinates, length of an ORF, TIS coordinate information, relative strand and location in contigs
typedef std::tuple<std::vector<int>, std::vector<indexPair>, size_t, std::vector<int>, std::vector<indexPair>, bool, std::pair<ContigLoc, bool>> ORFNodeVector;
// maps an ORFNodeVector sequence to its ID in order
typedef std::map<size_t, ORFNodeVector> ORFNodeMap;
// maps an map of ORFNodeVector sequence to its colour
typedef std::map<size_t, std::map<size_t, ORFNodeVector>> ColourORFMap;
// vector of ORF paths through graphs
typedef std::vector<ORFNodeVector> ORFVector;
// tuple for holding node information during traversal (1st = path index, 2nd = node id, 3rd = codon array, 4th = colour array, 5th = path length)
typedef std::tuple<size_t, int, uint8_t, boost::dynamic_bitset<>, size_t> NodeTuple;
// stack for holding nodes during DFS traversal for ORF identification
typedef std::stack<NodeTuple> NodeStack;
// stack for holding node ids, previous ORF ID, the next expected node, colour array and path length during DFS traversal for ORF ordering
typedef std::stack<std::tuple<int, boost::dynamic_bitset<>, size_t>> ORFStack;
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
// tuple of ORF sequence, node list, node coordinates for orientation and the contig locations if using FM index
typedef std::tuple<std::string, std::vector<int>, std::vector<std::vector<size_t>>, std::pair<ContigLoc, bool>> RefindTuple;
// map containing nodeID, search sequence and refind tuple
typedef std::map<int, std::map<std::string, RefindTuple>> RefindMap;

#endif //DEFINITIONS_H