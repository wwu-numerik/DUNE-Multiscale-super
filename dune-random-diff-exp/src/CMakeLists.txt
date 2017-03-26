file(GLOB SRC "*.h" "*.hh" "*.c" "*.cc")

add_subdirectory("Mlmc")
add_subdirectory("Permeability")
add_subdirectory("Utils")
add_subdirectory("randomDiff")

#set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -Wall -pedantic -Wshorten-64-to-32 -Wundef -Wnon-virtual-dtor -Wno-tautological-compare")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wfloat-equal -Wdisabled-macro-expansion  -Wused-but-marked-unused -Wnon-literal-null-conversion")

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-unused-local-typedef")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-order")

add_executable("dune-random-diff" dune-random-diff.cc ${SRC})

#DUNE default Library
# ********************************************************************************************
target_link_dune_default_libraries("dune-random-diff")
# ********************************************************************************************

#Additional Library
# ********************************************************************************************
#LIB FFT library
include_directories(${FFTW_INCLUDE_DIRS})
target_link_libraries("dune-random-diff" ${FFTW_LIBRARIES})
# ********************************************************************************************


#Includes
# ********************************************************************************************
include_directories(${CMAKE_SOURCE_DIR}/Mlmc)
include_directories(${CMAKE_SOURCE_DIR}/Permeability)
include_directories(${CMAKE_SOURCE_DIR}/randomDiff)
include_directories(${CMAKE_SOURCE_DIR}/Utils)
# ********************************************************************************************


# Set headers prop
# ********************************************************************************************
FILE( GLOB_RECURSE dune_header "/home/mercury/dune-prj/libdune2.5.0/include/*.hh")
FILE( GLOB_RECURSE local_dune "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )
set( COMMON_HEADER ${local_header} ${local_dune} ${dune_header} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )
# ********************************************************************************************

#disable warnings from dependent modules, except ours
# ********************************************************************************************
list(APPEND module_exceptions "")
foreach(_mod ${ALL_DEPENDENCIES})
  list (FIND module_exceptions "${_mod}" _index)
  if (_index LESS 0)
    dune_module_to_uppercase(_upper_case "${_mod}")
    if(${_mod}_INCLUDE_DIRS)
      foreach( _idir ${${_mod}_INCLUDE_DIRS} )
        add_definitions("-isystem ${_idir}")
      endforeach( _idir )
    endif(${_mod}_INCLUDE_DIRS)
  endif(_index GREATER -1)
endforeach(_mod DEPENDENCIES)
# ********************************************************************************************
