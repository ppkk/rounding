project(rounding)

set (CMAKE_CXX_FLAGS " -g -Wall -std=c++11 ${CMAKE_CXX_FLAGS}")
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

