cmake_minimum_required(VERSION 2.8)


option(DOWNLOAD_BOOST "Use this to download boost into this directory, instead of searching the system." OFF)


########## CMAKE SETUP ###########
#For building on much newer CMAKE
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

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


#The configure file from the old autotools build...
configure_file(config.h.in config.h)

#required packages
find_package(GSL REQUIRED)
find_package(NLOPT REQUIRED)

#
############## LIBRARIES ############
#
add_library(GS_tools STATIC Tools.cpp Tools.h)
target_link_libraries(GS_tools ${GSL_LIBRARIES})
#TODO: Don't we need to setup include directories for the .h files???
#TODO: Make this an object once we upgrade to a newer CMake

add_library(GS_annot STATIC SeqAnnotator.cpp SeqAnnotator.h)
target_link_libraries(GS_tools) #TODO: Figure out how to break this dependency


add_library(gem_everything DataSet.cpp ExprFunc.cpp ExprModel.cpp ExprPar.cpp ExprPredictor.cpp FactorIntFunc.cpp IO.cpp ObjFunc.cpp PredictorTrainer.cpp)
target_link_libraries(gem_everything GS_tools)
#target_link_libraries(gem_everything GS_annot)
target_link_libraries(gem_everything "nlopt")#should get includes
target_include_directories(gem_everything PUBLIC ${NLOPT_INCLUDE_DIRS})


#
############# BINARIES ####################
#
# add the executables
add_executable(seq2expr seq2expr.cpp)
target_link_libraries(seq2expr gem_everything GS_annot GS_tools)
target_link_libraries(seq2expr ${NLOPT_LIBRARIES} "nlopt")
target_link_libraries(seq2expr ${GSL_LIBRARIES})



add_executable(seqannot seqannot.cpp utils/gs_errors.h utils/gs_parsing.h)
target_link_libraries(seqannot GS_tools GS_annot)
