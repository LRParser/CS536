cmake_minimum_required(VERSION 3.5)
project(CG_hw2)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(SOURCE_FILES CG_hw2.cpp Node.cpp Node.h)
configure_file(cpts_test0.txt cpts_test0.txt COPYONLY)

add_executable(CG_hw2 ${SOURCE_FILES})
