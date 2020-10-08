//
// Created by sth19 on 25/08/2020.
//

#include "match_string.h"

std::unordered_map<std::string, bool> py_call_strings(const std::vector<std::string>& assembly_list,
                                                      const std::unordered_map<std::string, std::string>& query_dict,
                                                      const bool write_idx,
                                                      size_t num_threads)
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
    std::unordered_map<std::string, bool> result = call_strings(assembly_list, query_dict, write_idx, num_threads);

    // return C++ map (pybind11 converts to python dictionary)
    return result;
}

PYBIND11_MODULE(fmindex, m)
{
    m.doc() = "Finds presence/absence of gene sequences in source genomes";

    m.def("query_gene", &py_call_strings, "Generate dictionary for dictionary or queries, whether they represent true sequences.",
    py::arg("assembly_files"),
    py::arg("queries"),
    py::arg("write_idx") = 1,
    py::arg("threads") = 1);
}