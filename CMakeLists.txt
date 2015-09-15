project(chordalSubgraph)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(Boost_USE_STATIC_LIBS	OFF)
set(Boost_USE_MULTITHREADED	ON)
set(Boost_USE_STATIC_RUNTIME	OFF)

option(BUILD_R_PACKAGE "Build R package" ON)

add_definitions("-DBOOST_NO_AUTO_PTR")
if(CMAKE_COMPILER_IS_GNUCC)
	add_definitions("-Wall -fPIC")
endif()

#Find boost packages
find_package(Boost 1.53.0 COMPONENTS graph program_options system serialization filesystem regex REQUIRED)

add_subdirectory(chordalSubgraphCommon)
add_subdirectory(berry2006Test)
if(BUILD_R_PACKAGE)
	add_subdirectory(RPackage)
	find_package(R REQUIRED)
	if(WIN32)
		install(CODE "execute_process(COMMAND \"${R_COMMAND}\" CMD INSTALL --no-multiarch RPackage/\$\{CMAKE_INSTALL_CONFIG_NAME\} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})")
	else()
		install(CODE "execute_process(COMMAND \"${R_COMMAND}\" CMD INSTALL RPackage/ WORKING_DIRECTORY ${CMAKE_BINARY_DIR})")
	endif()
endif()