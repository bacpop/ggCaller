#ifndef CONCURRENT_HASH_MAP_H
#define CONCURRENT_HASH_MAP_H

#include "definitions.h"

// Intel tbb headers
#include <tbb/tbb.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/blocked_range.h>

// global variable declaration
using namespace tbb;


// A concurrent hash table that maps strings to ORFNodeVector.
typedef concurrent_hash_map<std::string, ORFNodeVector> ORFNodeMap;
// A concurrent hash table that node ids to bool (dummy variable).
typedef concurrent_hash_map<std::string, bool> NodeMap;
// A vector of concurrent hash table that node ids to bool (dummy variable).
typedef std::vector<NodeMap> NodeColourMap;

#endif //CONCURRENT_HASH_MAP
