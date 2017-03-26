#Try to find LibFFTW
#Once done this will define
# LIBFFTW_FOUND
# LIBFFTW_INCLUDE_DIR
# LIBFFTW_LIBRARIES
#
# Variables used by this module which you may want to set:
# FFTW_ROOT Path list to search for GMP

#User set Paths
find_path(FFTW_INCLUDE_DIR
  NAMES "fftw3-mpi.h" PATHS ${FFTW_PREFIX} ${FFTW_ROOT}
  PATH_SUFFIXES include NO_DEFAULT_PATH)

# try default paths now
find_path(FFTW_INCLUDE_DIRNAMES "fftw3-mpi.h")

# look for library FFTW3, only at positions given by the user
find_library(FFTW_LIB fftw3
  PATHS ${FFTW_PREFIX} ${FFTW_ROOT}
  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH
  DOC "FFTW library")

# try default paths now
find_library(FFTW_LIB fftw3)


# look for library FFTW3_MPI, only at positions given by the user
find_library(FFTW_MPI_LIB fftw3_mpi
  PATHS ${FFTW_PREFIX} ${FFTW_ROOT}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "FFTW3_mpi library")

# try default paths now
find_library(FFTW_MPI_LIB fftw3_mpi)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIR FFTW_LIB FFTW_MPI_LIB)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIB FFTW_MPI_LIB)


# if both headers and library are found, store results
if(FFTW_FOUND)
  set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
  set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTW_MPI_LIB})
  set(FFTW_COMPILE_FLAGS "-DENABLE_FFTW=1")
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of FFTW succeded:\n"
    "Include directory: ${FFTW_INCLUDE_DIRS}\n"
    "Library directory: ${FFTW_LIBRARIES}\n\n")

  # register in dune
#  dune_register_package_flags(COMPILE_DEFINITIONS ""
#      INCLUDE_DIRS "${FFTW_INCLUDE_DIRS}"
#      LIBRARIES "${FFTW_LIBRARIES}")


else(FFTW_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of FFTW failed:\n"
    "Include directory: ${FFTW_INCLUDE_DIR}\n"
    "fftw library directory: ${FFTW_LIB}\n"
    "fftw_mpi library directory: ${FFTW_MPI_LIB}\n\n")
endif(FFTW_FOUND)
