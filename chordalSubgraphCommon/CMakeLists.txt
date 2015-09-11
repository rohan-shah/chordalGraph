cmake_minimum_required(VERSION 3.1)
set(CppFiles cliqueTree.cpp)
set(HeaderFiles cliqueTree.h)

add_library(chordalSubgraphCommon STATIC ${CppFiles} ${HeaderFiles})
target_include_directories(chordalSubgraphCommon PUBLIC ${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/chordalSubgraphCommon ${IGRAPH_INCLUDES})
target_link_libraries(chordalSubgraphCommon PUBLIC  ${Boost_LIBRARIES})

target_compile_definitions(chordalSubgraphCommon PUBLIC BOOST_ALL_NO_LIB BOOST_ALL_DYN_LINK)
if(WIN32)
	set_target_properties(chordalSubgraphCommon PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else()
	target_compile_features(chordalSubgraphCommon PUBLIC cxx_rvalue_references)
endif()
