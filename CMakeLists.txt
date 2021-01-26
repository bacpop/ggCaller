cmake_minimum_required(VERSION 3.4)
project(ggCaller_cpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -march=native")

set(CMAKE_VERBOSE_MAKEFILE ON)

# find pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# find pybind and create python module
find_package(pybind11 REQUIRED)
pybind11_add_module(ggCaller_cpp src/call_ORFs.cpp src/ggCaller_bindings.cpp src/indexing.cpp src/match_string.cpp src/traversal.cpp src/unitigDict.cpp)

# check for conda environment
IF( DEFINED ENV{CONDA_PREFIX} )
  # set conda directory
  SET(CONDA_DIR "$ENV{CONDA_PREFIX}")
  # find Seqan3 packages
  find_path(SEQAN3_CLONE_DIR name bin PATHS ${CONDA_DIR} NO_DEFAULT_PATH)
  find_path(SEQAN3_INCLUDE_DIR name seqan PATHS "${SEQAN3_CLONE_DIR}/include" NO_DEFAULT_PATH)
  find_path(SEQAN3_SUBMODULES_DIR name submodules PATHS "${SEQAN3_INCLUDE_DIR}/seqan3" NO_DEFAULT_PATH)
  find_package(SeqAn3 REQUIRED HINT "${CONDA_PREFIX}/share/cmake/seqan3")
  # link bifrost, zlib, seqan and pthreads
  target_link_libraries(ggCaller_cpp PRIVATE seqan3::seqan3 ${CONDA_DIR}/lib/libbifrost.so ${CONDA_DIR}/lib/libz.so ${CONDA_DIR}/lib/libomp.so Threads::Threads)
ELSE()
  # find packages if conda not used
  find_package(SeqAn3 REQUIRED)
  find_library(z REQUIRED)
  find_library(bifrost REQUIRED)
  find_package(OpenMP)
  target_link_libraries(ggCaller_cpp PRIVATE seqan3::seqan3 bifrost z pthread)
  # link OpenMP if found
  if(OpenMP_CXX_FOUND)
    target_link_libraries(ggCaller_cpp PRIVATE OpenMP::OpenMP_CXX)
  endif()
ENDIF()