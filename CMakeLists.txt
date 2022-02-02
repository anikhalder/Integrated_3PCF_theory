cmake_minimum_required (VERSION 2.8.11)

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
project(Integrated_3PCF_theory)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -Wall -O3 -pthread -fopenmp -march=native -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O3 -pthread -fopenmp -march=native -fPIC")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -O4 -ffast-math -march=native -fopenmp -g -fpic")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Ofast -march=native -fopenmp -g -fPIC")

find_package(GSL REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include/ ${CMAKE_SOURCE_DIR}/external/class_files/class_include/ ${CMAKE_SOURCE_DIR}/external/cubature_files/cubature_include/)

add_subdirectory(src)