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


add_library(gem_everything DataSet.cpp ExprFunc.cpp ExprModel.cpp ExprPar.cpp ExprPredictor.cpp FactorIntFunc.cpp IO.cpp ObjFunc.cpp PredictorTrainer.cpp SeqAnnotator.cpp Tools.cpp)
target_link_libraries(gem_everything "nlopt")#should get includes
target_include_directories(gem_everything PUBLIC ${NLOPT_INCLUDE_DIRS})
# add the executables
add_executable(seq2expr seq2expr.cpp)
add_executable(seqannot seqannot.cpp)

target_link_libraries(seq2expr ${NLOPT_LIBRARIES} "nlopt")
target_link_libraries(seq2expr ${GSL_LIBRARIES})
target_link_libraries(seq2expr gem_everything)


target_link_libraries(seqannot ${NLOPT_LIBRARIES})
target_link_libraries(seqannot ${GSL_LIBRARIES})
target_link_libraries(seqannot gem_everything)
