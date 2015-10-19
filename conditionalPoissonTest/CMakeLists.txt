cmake_minimum_required(VERSION 3.1)
add_executable(conditionalPoissonTest conditionalPoissonTest.cpp)
target_link_libraries(conditionalPoissonTest PUBLIC chordalGraphCommon)

if(WIN32)
	set_target_properties(conditionalPoissonTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
