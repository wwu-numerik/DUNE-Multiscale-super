#set minimum cmake version
cmake_minimum_required(VERSION 3.1)
set( CMAKE_COLOR_MAKEFILE ON)
project(dune-random-diff CXX)
# ********************************************************************************************
#set(CMAKE_PREFIX_PATH /home/mercury/dune-prj/libdune2.4.1/lib/cmake)
#set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} /home/mercury/dune-prj/libdune2.4.1/lib/pkgconfig/)
# ********************************************************************************************

set(FFTW_ROOT $ENV{DUNE_FFTW_ROOT})

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("src")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules )

dune_enable_all_packages()


# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

# Set headers prop
# ********************************************************************************************
#FILE( GLOB_RECURSE dune_header "/home/mercury/dune-prj/libdune2.4.1/include/*.hh")
#FILE( GLOB_RECURSE local_dune "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )
#set( COMMON_HEADER ${local_header} ${local_dune} ${dune_header} )
#set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )
# ********************************************************************************************

#disable warnings from dependent modules, except ours
# ********************************************************************************************
#list(APPEND module_exceptions "")
#foreach(_mod ${ALL_DEPENDENCIES})
#  list (FIND module_exceptions "${_mod}" _index)
#  if (_index LESS 0)
#    dune_module_to_uppercase(_upper_case "${_mod}")
#    if(${_mod}_INCLUDE_DIRS)
#      foreach( _idir ${${_mod}_INCLUDE_DIRS} )
#        add_definitions("-isystem ${_idir}")
#      endforeach( _idir )
#    endif(${_mod}_INCLUDE_DIRS)
#  endif(_index GREATER -1)
#endforeach(_mod DEPENDENCIES)
# ********************************************************************************************
