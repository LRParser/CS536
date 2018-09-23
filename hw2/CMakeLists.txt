cmake_minimum_required(VERSION 3.5)
project(CG_hw2)

set(CMAKE_CXX_FLAGS "-I /usr/include/Inventor -lInventor -std=c++11")
add_library(Inventor SHARED IMPORTED)
set(SOURCE_FILES CG_hw2.cpp)
configure_file(cpts_test0.txt cpts_test0.txt COPYONLY)

add_executable(CG_hw2 ${SOURCE_FILES})
target_link_libraries(CG_hw2 /home/joe/open-inventor/usr/lib/libInventor.so)
