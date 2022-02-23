//
// Created by sth19 on 07/12/2021.
//

#ifndef GGCALLER_SEARCH_DBG_H
#define GGCALLER_SEARCH_DBG_H

#include "definitions.h"
#include "unitigDict.h"

std::unordered_set<int> query_DBG(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                    const std::string& query,
                                    const int& kmer,
                                    const double& id_cutoff);

#endif //GGCALLER_SEARCH_DBG_H
