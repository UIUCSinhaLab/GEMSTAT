cmake_minimum_required(VERSION 3.10)



#The configure file from the old autotools build...
configure_file(config.h.in config.h)
add_compile_definitions("BETAOPTTOGETHER=true")

#required packages
find_package(GSL REQUIRED)
find_package(NLOPT REQUIRED)



#
############## LIBRARIES ############
#

#imported package in our own source-tree
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/snot")

add_subdirectory("utils")

add_subdirectory("tools")

add_subdirectory("annotation")

add_subdirectory("aspects")


#Does this default to being an object on older CMake?
add_library(gem_everything DataSet.cpp ExprFunc.cpp ExprModel.cpp ExprPar.cpp ExprPredictor.cpp FactorIntFunc.cpp IO.cpp ObjFunc.cpp PredictorTrainer.cpp)
target_link_libraries(gem_everything GS_tools GS_annot)
target_link_libraries(gem_everything snot)
target_link_libraries(gem_everything "nlopt")#should get includes
target_include_directories(gem_everything PUBLIC ${NLOPT_INCLUDE_DIRS})

#Signalling extensions
add_library(GS_signalling DataSet_signal.cpp DataSet_signal.h)
target_link_libraries(GS_signalling GS_tools GS_annot)

#
############# BINARIES ####################
#
# add the executables
add_executable(seq2expr "seq2expr.cpp")
target_link_libraries(seq2expr main_aspects)
target_link_libraries(seq2expr gem_everything GS_annot GS_tools)
target_link_libraries(seq2expr GS_signalling)
target_link_libraries(seq2expr ${NLOPT_LIBRARIES} "nlopt")
target_link_libraries(seq2expr ${GSL_LIBRARIES})
