cmake_minimum_required(VERSION 3.4)
project(ggCaller_cpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -march=native")

#set(CMAKE_VERBOSE_MAKEFILE ON)

# find pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# find pybind and create python module
find_package(pybind11 REQUIRED)
pybind11_add_module(ggCaller_cpp src/bindings.cpp src/call_ORFs.cpp src/graph.cpp src/indexing.cpp src/match_string.cpp src/traversal.cpp src/unitigDict.cpp src/gene_overlap.cpp src/ORF_connection.cpp src/ORF_clustering.cpp src/edlib/edlib.cpp)

# check for conda environment
IF( DEFINED ENV{CONDA_PREFIX} )
  # set conda directory
  SET(CONDA_DIR "$ENV{CONDA_PREFIX}")
  # find Seqan3 packages
  find_path(SEQAN3_CLONE_DIR name bin PATHS ${CONDA_DIR} NO_DEFAULT_PATH)
  find_path(SEQAN3_INCLUDE_DIR name seqan PATHS ${SEQAN3_CLONE_DIR}/include NO_DEFAULT_PATH)
  find_path(SEQAN3_SUBMODULES_DIR name submodules PATHS ${SEQAN3_INCLUDE_DIR}/seqan3 NO_DEFAULT_PATH)
ENDIF()

#Find include dirs for z and bifrost
FIND_PATH(z_INCLUDE_DIR NAMES zlib.h HINTS ${CONDA_DIR}/include)
FIND_PATH(bifrost_INCLUDE_DIR bifrost/CompactedDBG.hpp HINTS ${CONDA_DIR}/include)

# Find libraries for z and bifrost
FIND_LIBRARY(z_LIBRARY NAMES z HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_LIBRARY(bifrost_LIBRARY NAMES bifrost HINTS ${CONDA_DIR}/lib REQUIRED)

# Set package found as true
FIND_PACKAGE_HANDLE_STANDARD_ARGS(z DEFAULT_MSG z_INCLUDE_DIR z_LIBRARY)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(bifrost DEFAULT_MSG bifrost_INCLUDE_DIR bifrost_LIBRARY)

# Set dirs if multiple include/libraries found
SET( z_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( z_LIBRARIES ${z_LIBRARY} )
SET( bifrost_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( bifrost_LIBRARIES ${bifrost_LIBRARY} )

# find Seqan3 package
find_package(SeqAn3 REQUIRED HINTS ${CONDA_DIR}/share/cmake/seqan3)

# link libraries
target_link_libraries(ggCaller_cpp PRIVATE seqan3::seqan3 ${bifrost_LIBRARIES} ${z_LIBRARIES} Threads::Threads)

#find OpenMP library
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
  target_link_libraries(ggCaller_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()
