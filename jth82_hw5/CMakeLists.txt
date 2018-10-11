cmake_minimum_required(VERSION 3.5)
project(CG_hw5)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(SOURCE_FILES CG_hw5.cpp Node.cpp Node.h)

add_executable(CG_hw5 ${SOURCE_FILES})
