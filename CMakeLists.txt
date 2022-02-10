cmake_minimum_required(VERSION 3.10)
project(ggCaller_cpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -march=native")

#set(CMAKE_VERBOSE_MAKEFILE ON)

# find packages with cmake methods
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED system serialization graph)
include_directories(${Boost_INCLUDE_DIRS})

# find pybind and create python module
find_package(pybind11 REQUIRED)
pybind11_add_module(ggCaller_cpp src/bindings.cpp src/call_ORFs.cpp src/graph.cpp src/indexing.cpp src/match_string.cpp src/traversal.cpp src/unitigDict.cpp src/gene_overlap.cpp src/ORF_connection.cpp src/ORF_clustering.cpp src/gene_refinding.cpp src/distances.cpp src/search_DBG.cpp src/edlib/edlib.cpp src/ORF_scoring.cpp src/gene_graph.cpp)

# check for conda environment
IF( DEFINED ENV{CONDA_PREFIX} )
  # set conda directory
  SET(CONDA_DIR "$ENV{CONDA_PREFIX}")
  SET(Protobuf_PROTOC_EXECUTABLE ${CONDA_DIR}/bin/protoc)
ENDIF()

# find correct protobuf library
FIND_LIBRARY(Protobuf_LIBRARY NAMES protobuf HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_LIBRARY(Protobuf_PROTOC_LIBRARY NAMES protoc HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_PATH(Protobuf_INCLUDE_DIR zlib.h HINTS ${CONDA_DIR}/include)

#Find include dirs
FIND_PATH(z_INCLUDE_DIR NAMES zlib.h HINTS ${CONDA_DIR}/include)
FIND_PATH(bifrost_INCLUDE_DIR bifrost/CompactedDBG.hpp HINTS ${CONDA_DIR}/include)
FIND_PATH(tbb_INCLUDE_DIR tbb/concurrent_hash_map.h HINTS ${CONDA_DIR}/include)

# Find libraries
FIND_LIBRARY(z_LIBRARY NAMES z HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_LIBRARY(bifrost_LIBRARY NAMES bifrost HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_LIBRARY(tbb_LIBRARY NAMES tbb libtbb.so libtbb.so.0 libtbb.so.1 libtbb.so.2 HINTS ${CONDA_DIR}/lib REQUIRED)

# Set package found as true
FIND_PACKAGE_HANDLE_STANDARD_ARGS(z DEFAULT_MSG z_INCLUDE_DIR z_LIBRARY)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(bifrost DEFAULT_MSG bifrost_INCLUDE_DIR bifrost_LIBRARY)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(tbb DEFAULT_MSG tbb_INCLUDE_DIR tbb_LIBRARY)

# Set dirs if multiple include/libraries found
SET( z_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( z_LIBRARIES ${z_LIBRARY} )
SET( tbb_INCLUDE_DIRS ${tbb_INCLUDE_DIR} )
SET( tbb_LIBRARIES ${tbb_LIBRARY} )
SET( bifrost_INCLUDE_DIRS ${z_INCLUDE_DIR} )
# hack to avoid local MKL being linked
SET( MKL_LIBRARIES ${bifrost_LIBRARY} )
# avoid local MKL
SET(MKL_FOUND TRUE)

# add torch libraries
FIND_PATH(PYTHON_SITE_PACKAGES LICENSE.txt HINTS ${CONDA_DIR}/lib NO_DEFAULT_PATH)
#MESSAGE(${PYTHON_SITE_PACKAGES})
list(APPEND CMAKE_PREFIX_PATH ${PYTHON_SITE_PACKAGES}/torch/share/cmake)
list(APPEND CMAKE_PREFIX_PATH ${MKL_ROOT})
find_package(Torch REQUIRED)

# add sdsl include directory
target_include_directories(ggCaller_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${tbb_INCLUDE_DIR})

# link libraries
target_link_libraries(ggCaller_cpp PRIVATE ${MKL_LIBRARIES} ${z_LIBRARIES} ${tbb_LIBRARIES} Threads::Threads Eigen3::Eigen ${Boost_LIBRARIES} ${TORCH_LIBRARIES})

#find OpenMP library
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
  target_link_libraries(ggCaller_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()
