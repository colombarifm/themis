# CMake project file for THEMIS

cmake_minimum_required(VERSION 2.6)
project (themis)
enable_language (Fortran CXX)
cmake_policy(SET CMP0015 NEW)


if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -g -O3 -ffpe-summary=none -fcheck=bounds -Wno-compare-reals -Wno-conversion -fbacktrace -fcheck=all -Wextra -no-pie")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -g -O3")
endif()


#
# Compile
#

link_directories(lib)

set(CMAKE_VERBOSE_MAKEFILE on)
set(EXECUTABLE_OUTPUT_PATH "bin")

file(GLOB_RECURSE special_src_files "xdr.f90")
set_property(SOURCE ${special_src_files} PROPERTY COMPILE_FLAGS -cpp)

file(GLOB_RECURSE sources src/mod*.f90 src/themis.f90)

add_executable(themis ${sources} ${special_src_files})

target_link_libraries(themis xdrfile) 
