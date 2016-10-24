cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)

find_package(libflatarray)
include_directories(${libflatarray_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -march=native -DHAVE__M128")
set(CMAKE_BUILD_TYPE "Release")

add_executable(soa_block.test soa_block.test.cpp)
add_executable(soa_compare.benchmark soa_compare.benchmark.cpp)
add_executable(udt.test udt.test.cpp)
add_executable(vec_offsets.test vec_offsets.test.cpp)
