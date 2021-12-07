//
// Created by sth19 on 07/12/2021.
//

#ifndef GGCALLER_SEARCH_DBG_H
#define GGCALLER_SEARCH_DBG_H

#include "definitions.h"

MappingCoords query_DBG(const ColoredCDBG<>& ccdbg,
                        const std::string& query,
                        const int& kmer,
                        const std::unordered_map<std::string, size_t>& kmer_map);

#endif //GGCALLER_SEARCH_DBG_H
