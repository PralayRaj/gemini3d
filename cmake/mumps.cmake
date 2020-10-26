# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
#
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
#
# easily obtain MUMPS without compiling:
# CentOS 6/7 EPEL: yum install mumps-devel
# Ubuntu / Debian: apt install libmumps-dev

# --- prereqs
include(${CMAKE_CURRENT_LIST_DIR}/lapack.cmake)

if(mpi)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack.cmake)
endif()
# --- MUMPS

if(mumps_external)
  include(${CMAKE_CURRENT_LIST_DIR}/mumps_external.cmake)
  return()
endif()

unset(_mumps_extra)

if(MUMPS_ROOT OR (DEFINED ENV{MUMPS_ROOT}) OR (CMAKE_Fortran_COMPILER_ID STREQUAL GNU))
  set(_comp ${arith})
  if(NOT mpi)
    list(APPEND _comp mpiseq)
  endif()

  if(autobuild)
    find_package(MUMPS COMPONENTS ${_comp})
  else()
    find_package(MUMPS COMPONENTS ${_comp} REQUIRED)
  endif()
else()
  message(VERBOSE "Skipping find_package(MUMPS)")
endif()

if(MUMPS_FOUND)
  set(mumps_external false CACHE BOOL "autobuild Mumps")
else()
  include(${CMAKE_CURRENT_LIST_DIR}/mumps_external.cmake)
  return()
endif()

if(metis)
  find_package(METIS REQUIRED)
  target_link_libraries(MUMPS::MUMPS INTERFACE METIS::METIS)
endif()

if(scotch)
  find_package(Scotch REQUIRED COMPONENTS ESMUMPS)
  target_link_libraries(MUMPS::MUMPS INTERFACE Scotch::Scotch)
endif()

# rather than appending libraries everywhere, just put them together here.
target_link_libraries(MUMPS::MUMPS INTERFACE SCALAPACK::SCALAPACK LAPACK::LAPACK)
if(OpenMP_FOUND)
  target_link_libraries(MUMPS::MUMPS INTERFACE OpenMP::OpenMP_Fortran OpenMP::OpenMP_C)
endif()

if(mumps_external OR scalapack_external OR lapack_external OR NOT mpi)
# pre-build checks can't be used when external library isn't built yet.
  return()
endif()

# -- minimal check that MUMPS is linkable
set(CMAKE_REQUIRED_LIBRARIES MUMPS::MUMPS MPI::MPI_Fortran)

check_fortran_source_compiles("
implicit none (type, external)
include '${arith}mumps_struc.h'
external :: ${arith}mumps
type(${arith}mumps_struc) :: mumps_par
end"
  MUMPS_link SRC_EXT f90)

if(NOT MUMPS_link)
  message(STATUS "MUMPS ${MUMPS_LIBRARIES} not working with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  if(NOT autobuild)
    message(FATAL_ERROR "autobuild=off, so cannot proceed")
  endif()
  include(${CMAKE_CURRENT_LIST_DIR}/mumps_external.cmake)
endif()
