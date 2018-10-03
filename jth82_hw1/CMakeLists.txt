cmake_minimum_required(VERSION 3.5)
project(CG_hw1)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(SOURCE_FILES CG_hw1.cpp Node.cpp Node.h)
configure_file(cpts_in.txt cpts_in.txt COPYONLY)

add_executable(CG_hw1 ${SOURCE_FILES})