//
// Created by sth19 on 07/12/2021.
//

#ifndef GGCALLER_SEARCH_DBG_H
#define GGCALLER_SEARCH_DBG_H

#include "definitions.h"

std::vector<std::tuple<std::string, bool, std::pair<size_t, size_t>>> query_DBG(const ColoredCDBG<>& ccdbg,
                                                                                const std::string& query,
                                                                                const int& kmer);

#endif //GGCALLER_SEARCH_DBG_H
