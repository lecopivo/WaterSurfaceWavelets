cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")

set(CMAKE_CXX_STANDARD 17)
#SET(CMAKE_CXX_COMPILER clang++-5.0)
#add_compile_options(-std=c++17)

find_package(Eigen REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(~/local/include/)

add_subdirectory(src)
