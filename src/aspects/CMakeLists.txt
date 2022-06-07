cmake_minimum_required(VERSION 3.10)

add_subdirectory("regularization")

add_library(main_aspects STATIC "aspects.cpp")
target_link_libraries(main_aspects PUBLIC regularization_aspect)
