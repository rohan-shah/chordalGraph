project(chordalSubgraph)
cmake_minimum_required(VERSION 3.1)
set(Boost_USE_STATIC_LIBS	OFF)
set(Boost_USE_MULTITHREADED	ON)
set(Boost_USE_STATIC_RUNTIME	OFF)

add_definitions("-DBOOST_NO_AUTO_PTR")
if(CMAKE_COMPILER_IS_GNUCC)
	add_definitions("-Wall -fPIC")
endif()

#Find boost packages
find_package(Boost 1.53.0 COMPONENTS graph program_options system serialization filesystem regex REQUIRED)

add_subdirectory(chordalSubgraphCommon)
add_subdirectory(berry2006Test)