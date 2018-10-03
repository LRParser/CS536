cmake_minimum_required(VERSION 3.5)
project(CG_hw3)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(SOURCE_FILES CG_hw3.cpp Node.cpp Node.h)
configure_file(patchPoints.txt patchPoints.txt COPYONLY)

add_executable(CG_hw3 ${SOURCE_FILES})
