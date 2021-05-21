#
#
#

## https://stackoverflow.com/questions/1208681/with-cmake-how-would-you-disable-in-source-builds

#
# This function will prevent in-source builds
function(PreventInTreeBuilds)
	get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
	get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

	find_path(CMAKE_LISTS_IN_BUILD "CMakeLists.txt" PATHS "${CMAKE_BINARY_DIR}"
			NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
		NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
	)
	if(CMAKE_LISTS_IN_BUILD)
		set(ATTEMPTED_IN_TREE_BUILD TRUE)
	endif()

	# disallow in-source builds
	if("${srcdir}" STREQUAL "${bindir}")
		set(ATTEMPTED_IN_TREE_BUILD TRUE)
	endif()

	if(ATTEMPTED_IN_TREE_BUILD)
		message("##########################")
		message("\nFATAL ERROR\n")
		message("	You are attempting an in-tree build.")
		message("	(Or to build in some other source tree...)")
		message("	Don't do that.

	ALWAYS build cmake projects into a new build directory,
	either a subdirectory of the project, or elsewhere.

### Now you need to remove CMakeCache.txt and CMakeFiles from the source tree. ###
### You can probably do this with \`git clean -dfX\` in the project root.      ###
	")
	message("##########################")
	#doesn't work. Can't remove these. (Probably actually written after the fatal error...)
	#file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeFiles" "${CMAKE_BINARY_DIR}/CMakeCache.txt")
	message(FATAL_ERROR "Attempted in-tree build.")
	endif()


endfunction()

PreventInTreeBuilds()
