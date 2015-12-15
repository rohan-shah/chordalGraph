if(USE_NAUTY)
	cmake_minimum_required(VERSION 3.1)
	add_executable(testComputeAutomorphismGroupSize test.cpp)
	target_link_libraries(testComputeAutomorphismGroupSize PUBLIC chordalGraphCommon)
	
	if(WIN32)
		set_target_properties(testComputeAutomorphismGroupSize PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
	endif()
endif()
