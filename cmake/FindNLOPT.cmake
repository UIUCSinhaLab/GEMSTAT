#stuff to find and download/check out nlopt or find one installed on the system.


#include(ExternalProject)
#ExternalProject_Add(nlopt_project
#  GIT_REPOSITORY https://github.com/stevengj/nlopt.git
#  STEP_TARGETS   build
#  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/nlopt
#  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#)
#ExternalProject_Get_Property(nlopt_project install_dir)
#include_directories(${install_dir}/include)

#message("NLOPT LIBS ")

#include(ExternalProject)
#ExternalProject_Add(foobar
#  GIT_REPOSITORY https://github.com/stevengj/nlopt.git
#  GIT_TAG        origin/release/1.2.3
#  STEP_TARGETS   build
#)

#TODO: detect from system, automatically download, etc. etc.
#
#https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
#
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/nlopt/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()


#TODO: Make sure to statically link

option(NLOPT_CXX ON)
option (NLOPT_FORTRAN OFF)
option (BUILD_SHARED_LIBS OFF) #We want to force static linkage
option (NLOPT_PYTHON OFF)
option (NLOPT_OCTAVE OFF)
option (NLOPT_MATLAB OFF)
option (NLOPT_GUILE OFF)
option (NLOPT_SWIG OFF)
option (NLOPT_TESTS OFF)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/nlopt ${CMAKE_BINARY_DIR}/nlopt)
#what about the .hpp?
