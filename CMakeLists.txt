cmake_minimum_required (VERSION 2.6)
project (APPSPAM)

SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -std=c++11 -g -Wno-write-strings -fopenmp")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
set(CMAKE_CXX_STANDARD 11)  # enable C++11 standard


# include
include_directories(./include ./src)

file(GLOB appspam_SOURCES "./src/*.cpp")

# target
add_executable(appspam ./src/Seed.cpp ./src/APPSPAM.cpp ./src/Pattern.cpp ./src/BucketManager.cpp ./src/ReadManager.cpp ./src/Tree.cpp ./src/Node.cpp ./src/GenomeManager.cpp ./src/Algorithms.cpp ./src/Scoring.cpp ./src/GlobalParameters.cpp ./src/SeqIO.cpp ./src/Placement.cpp ./src/Crc32.cpp)

FIND_PACKAGE( OpenMP)
      if(OPENMP_FOUND)
         message("OPENMP FOUND")
         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
         set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
      endif()
     include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})
     