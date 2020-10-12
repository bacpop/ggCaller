//
// Created by sth19 on 25/08/2020.
//

#include "match_string.h"

std::unordered_map<std::string, std::vector<std::string>> py_call_strings(const std::vector<std::string>& assembly_list,
                                                                          std::unordered_map<std::string, std::vector<std::string>>& query_list,
                                                                          const bool& write_idx,
                                                                          size_t& num_threads)
{
    // Check input

    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // Convert python objs to C++
    // Here done automatically with pybind11/stl.h

    // call pure C++ function
    query_list = call_strings(assembly_list, query_list, write_idx, num_threads);

    // return C++ map (pybind11 converts to python dictionary)
    return query_list;
}

PYBIND11_MODULE(match_string, m)
{
    m.doc() = "Finds presence/absence of gene sequences in source genomes";

    m.def("query_gene", &py_call_strings, "Removes values from a dictionary based on whether they are artifically generated sequences.",
    py::arg("assembly_files"),
    py::arg("queries"),
    py::arg("write_idx") = 1,
    py::arg("threads") = 1);
}