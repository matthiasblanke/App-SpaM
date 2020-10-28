cmake_minimum_required (VERSION 2.6)
project (APPSPAM)

SET(GCC_COMPILE_FLAGS "-O3 -std=c++11 -g -Wno-write-strings -fopenmp")
add_definitions(${GCC_COMPILE_FLAGS})
set(CMAKE_CXX_STANDARD 11)  # enable C++11 standard

# include
include_directories(./include ./src)

file(GLOB appspam_SOURCES "./src/*.cpp")

# target
add_executable(appspam ${appspam_SOURCES})

# OpenMP
FIND_PACKAGE(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})
