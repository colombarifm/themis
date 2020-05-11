# CMake project file for THEMIS


cmake_minimum_required(VERSION 2.6)
project (themis)
enable_language (Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -g -O3 -fcheck=bounds -Wno-compare-reals -Wno-conversion -fbacktrace -fcheck=all -Wextra -no-pie")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -g -O3")
endif()

#
# Compile
#

set(CMAKE_VERBOSE_MAKEFILE on)

file(GLOB_RECURSE sources src/*.f90)

set(EXECUTABLE_OUTPUT_PATH "bin")
add_executable(themis ${sources})
