cmake_minimum_required(VERSION 3.1)
add_executable(testMCMC2 test.cpp)
target_link_libraries(testMCMC2 PUBLIC chordalGraphCommon)

if(WIN32)
	set_target_properties(testMCMC2 PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
