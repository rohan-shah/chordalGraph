cmake_minimum_required(VERSION 3.1)
add_executable(stochasticEnumeration stochasticEnumeration.cpp)
target_link_libraries(stochasticEnumeration PUBLIC chordalGraphCommon)

if(WIN32)
	set_target_properties(stochasticEnumeration PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
