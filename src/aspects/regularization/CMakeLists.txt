cmake_minimum_required(VERSION 3.10)

add_library(regularization_aspect OBJECT "regularization.cpp" "regularization.hpp")
target_include_directories(regularization_aspect PUBLIC "../.." "../../annotation" "../../tools") #TODO: refactor everything
target_include_directories(regularization_aspect PUBLIC ".")
