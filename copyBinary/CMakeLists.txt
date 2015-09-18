cmake_minimum_required(VERSION 3.1)
function(copyBinary buildType destination resultName)
	set(COPIED_BINARY_FILES)
	#copy boost serialization
	get_filename_component(NAME_PART ${Boost_SERIALIZATION_LIBRARY_${buildType}} NAME_WE)
	get_filename_component(PATH_PART ${Boost_SERIALIZATION_LIBRARY_${buildType}} PATH)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${PATH_PART}/${NAME_PART}.dll" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#copy boost program_options
	get_filename_component(NAME_PART ${Boost_PROGRAM_OPTIONS_LIBRARY_${buildType}} NAME_WE)
	get_filename_component(PATH_PART ${Boost_PROGRAM_OPTIONS_LIBRARY_${buildType}} PATH)
	add_custom_command(OUTPUT ${destination}/${NAME_PART}.dll COMMAND ${CMAKE_COMMAND} ARGS -E copy "${PATH_PART}/${NAME_PART}.dll" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#copy boost graph
	get_filename_component(NAME_PART ${Boost_GRAPH_LIBRARY_${buildType}} NAME_WE)
	get_filename_component(PATH_PART ${Boost_GRAPH_LIBRARY_${buildType}} PATH)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${PATH_PART}/${NAME_PART}.dll" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#copy boost regex
	get_filename_component(NAME_PART ${Boost_REGEX_LIBRARY_${buildType}} NAME_WE)
	get_filename_component(PATH_PART ${Boost_REGEX_LIBRARY_${buildType}} PATH)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${PATH_PART}/${NAME_PART}.dll" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#Now for MPFR
	get_target_property(MPFR_LOCATION mpfr "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${MPFR_LOCATION}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${MPFR_LOCATION}" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#Now for MPIR
	get_target_property(MPIR_LOCATION mpir "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${MPIR_LOCATION}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${MPIR_LOCATION}" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	set(${resultName} ${COPIED_BINARY_FILES} PARENT_SCOPE)
endfunction()
function(copyQt buildType destination resultName)
	set(COPIED_BINARY_FILES)

	get_target_property(QtCore_location Qt5::Core "LOCATION_${buildType}")
	get_filename_component(NAME_PART ${QtCore_location} NAME_WE)
	get_filename_component(QtBinaryPath ${QtCore_location} PATH)
	add_custom_command(OUTPUT ${destination}/${NAME_PART}.dll COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtCore_location} ${destination}/${NAME_PART}.dll)
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	if(EXISTS ${QtBinaryPath}/icuin51.dll)
		add_custom_command(OUTPUT ${destination}/icuin51.dll COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtBinaryPath}/icuin51.dll "${destination}/icuin51.dll")
		add_custom_command(OUTPUT ${destination}/icuuc51.dll COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtBinaryPath}/icuuc51.dll "${destination}/icuuc51.dll")
		add_custom_command(OUTPUT ${destination}/icudt51.dll COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtBinaryPath}/icudt51.dll "${destination}/icudt51.dll")
		list(APPEND COPIED_BINARY_FILES "${destination}/icuin51.dll" "${destination}/icuuc51.dll" "${destination}/icudt51.dll")
	endif()

	get_target_property(QtWidgets_location Qt5::Widgets "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${QtWidgets_location}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtWidgets_location} ${destination}/${NAME_PART}.dll)
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	get_target_property(QtGuiEgl_location Qt5::Gui_EGL "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${QtGuiEgl_location}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtGuiEgl_location} ${destination}/${NAME_PART}.dll)
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	get_target_property(QtGuiGlesV2_location Qt5::Gui_GLESv2 "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${QtGuiGlesV2_location}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtGuiGlesV2_location} ${destination}/${NAME_PART}.dll)
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	get_target_property(QtGui_location Qt5::Gui "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${QtGui_location}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy ${QtGui_location} ${destination}/${NAME_PART}.dll)
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#Now for platform stuff
	get_filename_component(QT_LIB_DIR "${QtGui_location}" DIRECTORY)
	if(EXISTS ${QT_LIB_DIR}/../plugins/platforms/qwindowsd.dll)
		add_custom_command(OUTPUT "${destination}/platforms/qwindowsd.dll" COMMAND ${CMAKE_COMMAND} -E copy "${QT_LIB_DIR}/../plugins/platforms/qwindowsd.dll" "${destination}/platforms/qwindowsd.dll")
		list(APPEND COPIED_BINARY_FILES "${destination}/platforms/qwindowsd.dll")
	endif()
	if(EXISTS "${QT_LIB_DIR}/../plugins/platforms/qwindows.dll")
		add_custom_command(OUTPUT "${destination}/platforms/qwindows.dll" COMMAND ${CMAKE_COMMAND} -E copy "${QT_LIB_DIR}/../plugins/platforms/qwindows.dll" "${destination}/platforms/qwindows.dll")
		list(APPEND COPIED_BINARY_FILES "${destination}/platforms/qwindows.dll")
	endif()

	set(${resultName} ${COPIED_BINARY_FILES} PARENT_SCOPE)
endfunction()
if(WIN32)
	#Copy files for the command line executables
	set(COPIED_BINARY_FILES)
	#If we're on windows, Start copying dlls
	set(BUILD_TYPES Release Debug RelWithDebInfo MinSizeRel)
	set(BOOST_BUILD_TYPES RELEASE DEBUG RELEASE RELEASE)
	foreach(val RANGE 3)
		list(GET BUILD_TYPES ${val} BUILD_TYPE)
		list(GET BOOST_BUILD_TYPES ${val} BOOST_BUILD_TYPE)
		
		copyBinary("${BOOST_BUILD_TYPE}" "${CMAKE_BINARY_DIR}/${BUILD_TYPE}" additionalBinaryFiles)
		list(APPEND COPIED_BINARY_FILES ${additionalBinaryFiles})
	endforeach()

	if(BUILD_R_PACKAGE)
		foreach(val RANGE 3)
			list(GET BUILD_TYPES ${val} BUILD_TYPE)
			list(GET BOOST_BUILD_TYPES ${val} BOOST_BUILD_TYPE)
		
			copyBinary("${BOOST_BUILD_TYPE}" "${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/inst/libs/x64/" additionalBinaryFiles)
			list(APPEND COPIED_BINARY_FILES ${additionalBinaryFiles})
		endforeach()
	endif()
	add_custom_target(copyBinary ALL DEPENDS ${COPIED_BINARY_FILES})
endif()