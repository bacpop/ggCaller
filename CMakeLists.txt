cmake_minimum_required (VERSION 3.4)
project(match_string)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)

find_package(pybind11 REQUIRED)
pybind11_add_module(match_string src/match_string.cpp src/match_bindings.cpp)

find_package(SeqAn3 3.0.0 REQUIRED)
target_link_libraries(match_string PRIVATE seqan3::seqan3 pthread)

if(OpenMP_CXX_FOUND)
    target_link_libraries(match_string PRIVATE OpenMP::OpenMP_CXX)
endif()