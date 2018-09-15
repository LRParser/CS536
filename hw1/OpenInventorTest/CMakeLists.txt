cmake_minimum_required(VERSION 3.5)
project(dummy)

set(CMAKE_CXX_FLAGS "-I /usr/include/Inventor -lInventor -std=c++11")
add_library(Inventor SHARED IMPORTED)
set(SOURCE_FILES OITest.cpp)

add_executable(dummy ${SOURCE_FILES})
target_link_libraries(dummy /home/joe/open-inventor/usr/lib/libInventor.so)
