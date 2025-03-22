# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2021 Vissarion Fisikopoulos
# Copyright (c) 2018-2021 Apostolos Chalkis
# Contributed and/or modified by Vaibhav Thakkar
# Licensed under GNU LGPL.3, see LICENCE file



CMAKE_MINIMUM_REQUIRED(VERSION 3.11)
project( VolEsti )

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from ../../external" OFF)

  find_library(GMP_LIB gmp)
  if (NOT GMP_LIB)
    message(FATAL_ERROR "This program requires the gmp library, and will not be compiled.")
  else()
    message(STATUS "Library gmp found: ${GMP_LIB}")
  endif(NOT GMP_LIB)
if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
else()
  find_library(IFOPT NAMES libifopt_core.so PATHS /usr/local/lib)
  find_library(IFOPT_IPOPT NAMES libifopt_ipopt.so PATHS /usr/local/lib)
  find_library(GMP NAMES libgmp.a PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(MPSOLVE NAMES libmps.so PATHS /usr/local/lib)
  find_library(PTHREAD NAMES libpthread.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(FFTW3 NAMES libfftw3.so.3 PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)



  if (NOT IFOPT)

    message(FATAL_ERROR "This program requires the ifopt library, and will not be compiled.")

  elseif (NOT GMP)

    message(FATAL_ERROR "This program requires the gmp library, and will not be compiled.")

  elseif (NOT MPSOLVE)

    message(FATAL_ERROR "This program requires the mpsolve library, and will not be compiled.")

  elseif (NOT FFTW3)

    message(FATAL_ERROR "This program requires the fftw3 library, and will not be compiled.")

  else()
    message(STATUS "Library ifopt found: ${IFOPT}")
    message(STATUS "Library gmp found: ${GMP}")
    message(STATUS "Library mpsolve found: ${MPSOLVE}")
    message(STATUS "Library fftw3 found:" ${FFTW3})

  endif(NOT IFOPT)

endif(DISABLE_NLP_ORACLES)

find_library(CDDLIB NAMES libcddgmp.a PATHS /usr/local/lib /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)

if (NOT CDDLIB)
  message(FATAL_ERROR "This program requires the cddlib library, and will not be compiled.")
else()
  message(STATUS "Library cddlib found: ${CDDLIB}")
  include_directories(/usr/include/cdd)
endif()

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()


set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

include_directories (BEFORE ../../external)
include_directories (BEFORE ../../external/minimum_ellipsoid)
include_directories (BEFORE ../../include/generators)
include_directories (BEFORE ../../include/volume)
include_directories (BEFORE ../../include)
include_directories (BEFORE ../../include/convex_bodies)
include_directories (BEFORE ../../include/random_walks)
include_directories (BEFORE ../../include/annealing)
include_directories (BEFORE ../../include/ode_solvers)
include_directories (BEFORE ../../include/root_finders)
include_directories (BEFORE ../../include/samplers)
include_directories (BEFORE ../../include/lp_oracles)
include_directories (BEFORE ../../include/nlp_oracles)
include_directories (BEFORE ../../include/misc)
include_directories (BEFORE ../../include/optimization)

# for Eigen
if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
  add_compile_options(-D "EIGEN_NO_DEBUG")
else ()
  add_compile_definitions("EIGEN_NO_DEBUG")
endif ()

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++17")
add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
#add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")


add_executable (volume volume.cpp)
TARGET_LINK_LIBRARIES(volume lp_solve ${CDDLIB} cdd ${GMP_LIB})

add_executable (sample sampler.cpp)
TARGET_LINK_LIBRARIES(sample lp_solve ${CDDLIB} cdd ${GMP_LIB})
# target_link_libraries(${PROJECT_NAME} PUBLIC cdd::cdd)
