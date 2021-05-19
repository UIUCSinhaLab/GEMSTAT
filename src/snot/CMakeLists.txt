cmake_minimum_required(VERSION 2.8)


########## CMAKE SETUP ###########
#For building on much newer CMAKE
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(SNOT-gemstat VERSION 0.1.0 LANGUAGES CXX)

add_library(snot STATIC snot.tab.cpp snot.tab.hpp snot.lex.cpp snot.yy.hpp param_storage.h)
