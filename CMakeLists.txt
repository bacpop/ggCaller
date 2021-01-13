cmake_minimum_required(VERSION 3.4)
project(ggCaller_cpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -march=native")

set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(OpenMP)
find_library(pthread REQUIRED)
find_library(z REQUIRED)
find_library(bifrost REQUIRED)

find_package(pybind11 REQUIRED)
pybind11_add_module(ggCaller_cpp src/call_ORFs.cpp src/ggCaller_bindings.cpp src/indexing.cpp src/match_string.cpp src/traversal.cpp src/unitigDict.cpp src/gene_overlap.cpp)

IF( DEFINED ENV{CONDA_PREFIX} )
  SET(CONDA_DIR "$ENV{CONDA_PREFIX}")
  find_path(SEQAN3_CLONE_DIR name bin PATHS ${CONDA_DIR} NO_DEFAULT_PATH)
  find_path(SEQAN3_INCLUDE_DIR name seqan PATHS "${SEQAN3_CLONE_DIR}/include" NO_DEFAULT_PATH)
  find_path(SEQAN3_SUBMODULES_DIR name submodules PATHS "${SEQAN3_INCLUDE_DIR}/seqan3" NO_DEFAULT_PATH)
  find_package(SeqAn3 REQUIRED HINT "${CONDA_PREFIX}/share/cmake/seqan3")
ELSE()
  find_package(SeqAn3 REQUIRED)
ENDIF()

target_link_libraries(ggCaller_cpp PRIVATE seqan3::seqan3 pthread z bifrost)

if(OpenMP_CXX_FOUND)
    target_link_libraries(ggCaller_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()