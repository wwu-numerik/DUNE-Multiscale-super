if(NOT dune-random-diff_FOUND)

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was dune-random-diff-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#report other information
set_and_check(dune-random-diff_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-random-diff_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set(dune-random-diff_CXX_FLAGS " -std=c++14 ")
set(dune-random-diff_CXX_FLAGS_DEBUG "-g3 -gdwarf-2")
set(dune-random-diff_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-random-diff_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-random-diff_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-random-diff_DEPENDS "dune-common;dune-geometry;dune-localfunctions;dune-grid;dune-typetree;dune-istl;dune-pdelab")
set(dune-random-diff_SUGGESTS "")
set(dune-random-diff_MODULE_PATH "${PACKAGE_PREFIX_DIR}/share/dune/cmake/modules")
set(dune-random-diff_LIBRARIES "")

# Lines that are set by the CMake buildsystem via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-random-diff_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-random-diff-targets.cmake")
endif(dune-random-diff_LIBRARIES)
endif(NOT dune-random-diff_FOUND)
