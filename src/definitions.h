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
#include <memory>
#include <memory>

// openMP headers
#include <omp.h>

// robinhood hashing header
#include "robin_hood.h"

// sdsl header
#include <sdsl/suffix_arrays.hpp>

// boost algorithm headers
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>

// boost serialising headers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/dynamic_bitset/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/map.hpp>
#include <serialize_tuple.h>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Eigen header
#include <Eigen/Sparse>

// edlib header
#include "edlib/edlib.h"

// Intel tbb headers
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>

// bifrost header
#include <bifrost/ColoredCDBG.hpp>

// global variable declaration
namespace py = pybind11;

// UnitigDict typedefs
// Vector of neighbouring nodes by ID, map of stop codon frames and set of colours in which edge is found in
typedef std::vector<std::tuple<int, std::vector<uint8_t>, std::unordered_set<size_t>>> NeighbourVector;

// Eigen typedef
typedef Eigen::Triplet<double> ET;

typedef sdsl::csa_wt<> fm_index_coll;

// general typedefs
// hasher using robin_hood hash
typedef robin_hood::hash<std::string> hasher;
// mapping of each colour to component nodes in graph
typedef std::vector<std::vector<size_t>> NodeColourVector;
//a pair of start and end coordinates for an ORF across a node
typedef std::pair<unsigned int, unsigned int> indexPair;
// pair that describes the contig locations of an ORF, 1-indexed for contig id (first) and locations within contig (second)
typedef std::pair<size_t, std::pair<int, int>> ContigLoc;
// tuple holding ORF path ID, nodes traversed, node coordinates, coordinates in path, 5p and 3p coordinates
typedef std::pair<std::vector<int>, std::vector<indexPair>> ORFCoords;
// tuple containing a vector of nodeIDs, a vector of start,stop and length coordinates, length of an ORF, relative strand and score and strings for protein/DNA
typedef std::tuple<std::vector<int>, std::vector<indexPair>, size_t, bool, float, std::string, std::string> ORFNodeVector;
// maps an ORFNodeVector sequence to its ID in order
typedef std::map<int, ORFNodeVector> ORFNodeMap;
// maps an map of ORFNodeMap to its colour
typedef std::map<size_t, ORFNodeMap> ColourORFMap;
// maps colours to the edges for every ORF in that colour
typedef std::map<size_t, std::unordered_map<size_t, std::unordered_set<size_t>>> ColourEdgeMap;
// map of ORF paths through graphs
typedef robin_hood::unordered_map<size_t, ORFNodeVector> ORFNodeRobMap;
// map of ORFNodeRobMap, one for each colour
typedef std::unordered_map<size_t, ORFNodeRobMap> ColourORFVectorMap;
// tuple for holding node information during traversal (1st = path index, 2nd = node id, 3rd = codon array, 4th = colour array, 5th = path length)
typedef std::tuple<size_t, int, std::bitset<3>, boost::dynamic_bitset<>, size_t> NodeTuple;
// stack for holding nodes during DFS traversal for ORF identification
typedef std::stack<NodeTuple> NodeStack;
// stack for holding pos_idx, node ids, previous ORF ID, the next expected node, colour array and path length during DFS traversal for ORF ordering
typedef std::stack<std::tuple<size_t, int, boost::dynamic_bitset<>, size_t>> ORFStack;
// A vector of paths
typedef std::vector<std::vector<int>> PathVector;
// A vector of paths, with associated lengths
typedef std::vector<std::pair<size_t, std::vector<int>>> RefindPathVector;
// mapping of node ID to a orientation for a specific strand, used in overlap analysis
typedef std::map<size_t, bool> NodeStrandMap;
// mapping of overlapping ORFs, detailed by ORFIDMap
typedef std::unordered_map<size_t, std::unordered_map<size_t, std::pair<char, size_t>>> ORFOverlapMap;
// tuple containing grouping information for ORFs filtered by Balrog
typedef std::pair<std::vector<std::tuple<int, int, size_t, size_t, std::shared_ptr<std::string>>>, std::vector<std::pair<size_t, std::pair<size_t, size_t>>>> ORFGroupPair;
// map of ORFs to clusters, with centroid as first entry
typedef std::unordered_map<std::string, std::vector<std::pair<size_t, size_t>>> ORFClusterMap;
// tuple of ORF sequence, node list, node coordinates for orientation and the contig locations if using FM index
typedef std::tuple<std::string, std::vector<int>, std::vector<std::vector<size_t>>, bool> RefindTuple;
// map containing nodeID, search sequence and refind tuple
typedef std::map<int, std::vector<RefindTuple>> RefindMap;
// dictionary passed by python for refinding
typedef std::unordered_map<int, std::pair<std::pair<std::vector<int>, std::vector<indexPair>>, std::vector<size_t>>> NodeSearchDict;

#endif //DEFINITIONS_H