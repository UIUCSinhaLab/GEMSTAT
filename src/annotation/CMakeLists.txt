cmake_minimum_required(VERSION 3.10)

add_library(GS_annot STATIC SeqAnnotator.cpp SeqAnnotator.h)
target_link_libraries(GS_annot PRIVATE GS_tools) #TODO: Figure out how to break this dependency
target_link_libraries(GS_annot PUBLIC GS_utils )
target_include_directories(GS_annot PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(seqannot seqannot.cpp )
target_link_libraries(seqannot GS_utils GS_tools GS_annot)
