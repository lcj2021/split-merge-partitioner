cmake_minimum_required (VERSION 2.8)
project (FSM_partitioner)

set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_CXX_STANDARD 17)

find_package (glog REQUIRED)
find_package (gflags REQUIRED)

# SET a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # SET the possible values of build type for cmake-gui
endif()

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_EXE_LINK_FLAGS_RELEASE} -Ofast -DNDEBUG")
# 

FIND_PACKAGE(OpenMP REQUIRED)
INCLUDE_DIRECTORIES(${OpenMP_CXX_INCLUDE_DIRS})

# CMake will automatically detect the depending .hpp of the corresponding .cpp 
# and add them to the building phase 
add_executable (
    main
    src/main.cpp
    src/graph.cpp
    src/hep_graph.cpp
    src/fsm_partitioner.cpp
    src/ne_partitioner.cpp
    src/hep_partitioner.cpp
    src/ebv_partitioner.cpp
    src/dbh_partitioner.cpp
    src/hdrf_partitioner.cpp
    src/fennel_partitioner.cpp
    src/bpart_partitioner.cpp
    src/hybridbl_partitioner.cpp
    src/hybrid_partitioner.cpp
    src/vertex2edgepart.cpp
    # src/edgelist2adjlist.hpp
    # src/test.cpp
    src/conversions.cpp
)


add_executable (
    get_stat
    src/get_stat.cpp
)

target_link_libraries(
    main 
    glog 
    gflags
    OpenMP::OpenMP_CXX
)
