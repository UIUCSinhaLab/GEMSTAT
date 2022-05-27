
find_package(GSL REQUIRED)

#
############## LIBRARIES ############
#
add_library(GS_tools STATIC Tools.cpp Tools.h)
target_link_libraries(GS_tools ${GSL_LIBRARIES})
target_include_directories(GS_tools PUBLIC ".")


#TODO: Don't we need to setup include directories for the .h files???
#TODO: Make this an object once we upgrade to a newer CMake
