cmake_minimum_required(VERSION 3.5)
project(CG_hw4)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(SOURCE_FILES CG_hw4.cpp Node.cpp Node.h)

add_executable(CG_hw4 ${SOURCE_FILES})
