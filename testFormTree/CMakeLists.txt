cmake_minimum_required(VERSION 3.1)
add_executable(testFormTree testFormTree.cpp)
target_link_libraries(testFormTree PUBLIC chordalGraphCommon)

if(WIN32)
	set_target_properties(testFormTree PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
