cmake_minimum_required (VERSION 3.4)
project(match_string)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)

find_package(pybind11 REQUIRED)
pybind11_add_module(match_string src/match_string.cpp src/match_bindings.cpp)

IF( DEFINED ENV{CONDA_PREFIX} )
  SET(CONDA_DIR "$ENV{CONDA_PREFIX}")
  find_path(SEQAN3_CLONE_DIR name bin PATHS ${CONDA_DIR} NO_DEFAULT_PATH)
  find_path(SEQAN3_INCLUDE_DIR name seqan PATHS "${SEQAN3_CLONE_DIR}/include" NO_DEFAULT_PATH)
  find_path(SEQAN3_SUBMODULES_DIR name submodules PATHS "${SEQAN3_INCLUDE_DIR}/seqan3" NO_DEFAULT_PATH)
  find_package(SeqAn3 REQUIRED HINT "${CONDA_PREFIX}/share/cmake/seqan3")
ELSE()
  find_package(SeqAn3 REQUIRED)
ENDIF()

target_link_libraries(match_string PRIVATE seqan3::seqan3 pthread)

if(OpenMP_CXX_FOUND)
    target_link_libraries(match_string PRIVATE OpenMP::OpenMP_CXX)
endif()