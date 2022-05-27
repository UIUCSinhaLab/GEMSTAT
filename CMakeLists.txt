cmake_minimum_required(VERSION 3.10)


option(DOWNLOAD_BOOST "Use this to download boost into this directory, instead of searching the system." OFF)


########## CMAKE SETUP ###########
#For building on much newer CMAKE
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif(POLICY CMP0079)

if(POLICY CMP0022)
    cmake_policy(SET CMP0022 NEW)
endif(POLICY CMP0022)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )
include("PreventInTreeBuilds")



#
######### PROJECT ###########
# set the project name
project(GEMSTAT VERSION 1.9.0 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


add_subdirectory("src")
