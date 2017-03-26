file(GLOB SRC "*.h" "*.hh" "*.c" "*.cc")

add_subdirectory("Mlmc")
add_subdirectory("Permeability")
add_subdirectory("Utils")
add_subdirectory("randomDiff")

set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


#set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -Wall -pedantic -Wundef -Wnon-virtual-dtor -Wno-tautological-compare")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wfloat-equal")

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
