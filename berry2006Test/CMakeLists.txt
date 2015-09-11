cmake_minimum_required(VERSION 3.1)
add_executable(berry2006Test test.cpp)
target_link_libraries(berry2006Test PUBLIC chordalSubgraphCommon)

if(WIN32)
	set_target_properties(berry2006Test PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
