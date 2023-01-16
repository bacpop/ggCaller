cmake_minimum_required(VERSION 3.10)
project(ggCaller_cpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

set(CMAKE_VERBOSE_MAKEFILE ON)

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
  # find torch cmake files
  SET(Torch_SUFF "site-packages/torch/share/cmake/Torch")
  FIND_PATH(Torch_DIR NAMES TorchConfig.cmake HINTS ${CONDA_DIR}/lib/python3.8/${Torch_SUFF} ${CONDA_DIR}/lib/python3.9/${Torch_SUFF} ${CONDA_DIR}/lib/python3.10/${Torch_SUFF} ${CONDA_DIR}/lib/python3.11/${Torch_SUFF} ${CONDA_DIR}/lib/python3.12/${Torch_SUFF} ${CONDA_DIR}/lib/python3.13/${Torch_SUFF} ${CONDA_DIR}/lib/python3.14/${Torch_SUFF} ${CONDA_DIR}/lib/python3.15/${Torch_SUFF} NO_DEFAULT_PATH)
  #MESSAGE(${Torch_DIR})
  #MESSAGE(${CONDA_DIR})
ENDIF()

# find correct protobuf library
FIND_LIBRARY(Protobuf_LIBRARY NAMES protobuf HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_LIBRARY(Protobuf_PROTOC_LIBRARY NAMES protoc HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_PATH(Protobuf_INCLUDE_DIR zlib.h HINTS ${CONDA_DIR}/include)

include(FindPackageHandleStandardArgs)
#Find zlib
FIND_PATH(z_INCLUDE_DIR NAMES zlib.h HINTS ${CONDA_DIR}/include)
FIND_LIBRARY(z_LIBRARY NAMES z HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(z DEFAULT_MSG z_INCLUDE_DIR z_LIBRARY)
SET( z_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( z_LIBRARIES ${z_LIBRARY} )

#Find bifrost
FIND_PATH(bifrost_INCLUDE_DIR bifrost/CompactedDBG.hpp HINTS ${CONDA_DIR}/include)
FIND_LIBRARY(bifrost_LIBRARY NAMES bifrost HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(bifrost DEFAULT_MSG bifrost_INCLUDE_DIR bifrost_LIBRARY)
SET( bifrost_INCLUDE_DIRS ${bifrost_INCLUDE_DIR} )
SET( bifrost_LIBRARIES ${bifrost_LIBRARY} )

# find TBB if not found already
FIND_PATH(tbb_INCLUDE_DIR concurrent_unordered_map.h HINTS ${CONDA_DIR}/include/tbb NO_DEFAULT_PATH)
FIND_LIBRARY(tbb_LIBRARY NAMES tbb libtbb.so libtbb.so.0 libtbb.so.1 libtbb.so.2 HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(tbb DEFAULT_MSG tbb_INCLUDE_DIR tbb_LIBRARY)
SET( tbb_INCLUDE_DIRS ${tbb_INCLUDE_DIR} )
SET( tbb_LIBRARIES ${tbb_LIBRARY} )

# find MKL
FIND_PATH(MKL_INCLUDE_DIR NAMES mkl.h HINTS ${MKL_ROOT}/include ${CONDA_DIR}/include)
FIND_LIBRARY(MKL_LIBRARY NAMES libmkl_core.so HINTS ${MKL_ROOT}/lib/intel64 ${CONDA_DIR}/lib REQUIRED)
SET(MKL_LIBRARIES)
list(APPEND MKL_LIBRARIES ${MKL_LIBRARY})
MESSAGE(${MKL_LIBRARIES})
MESSAGE(${MKL_INCLUDE_DIR})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARIES)
# avoid local MKL when looking for pytorch
SET(MKL_FOUND TRUE)

# add torch libraries
find_package(Torch REQUIRED)

# add sdsl include directory
target_include_directories(ggCaller_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${tbb_INCLUDE_DIR} ${MKL_INCLUDE_DIR})

# link libraries
target_link_libraries(ggCaller_cpp PRIVATE ${bifrost_LIBRARIES} ${MKL_LIBRARIES} ${z_LIBRARIES} ${tbb_LIBRARIES} Threads::Threads Eigen3::Eigen ${Boost_LIBRARIES} ${TORCH_LIBRARIES} )

# find openmp
if(DEFINED ENV{CONDA_PREFIX} AND (NOT APPLE OR CMAKE_COMPILER_IS_GNUCC))
    target_link_libraries(ggCaller_cpp PRIVATE gomp)
else()
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(ggCaller_cpp PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()