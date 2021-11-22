#ifndef GGCALLER_DISTANCES_H
#define GGCALLER_DISTANCES_H

#include "unitigDict.h"
#include <boost/dynamic_bitset.hpp>

template <typename T>
std::vector<T> combine_vectors(const std::vector<std::vector<T>> &vec,
                               const size_t len);

std::vector<double> get_distances_align(const std::vector<std::string>& matrix_in,
                                        const int n_threads);

std::vector<double> get_distances_pa(const std::vector<std::vector<bool>>& matrix_in,
                                    const int n_threads);

#endif //GGCALLER_DISTANCES_H

