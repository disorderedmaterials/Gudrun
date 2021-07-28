# CMake project file for Gudrun

cmake_minimum_required(VERSION 3.1)
project(Gudrun)
enable_language(Fortran)

# Build release version
set(CMAKE_BUILD_TYPE "Release")

include(ExternalProject)

# ZLIB
# Locate installed library, or build as external version (if BUILD_ZLIB)
option(BUILD_ZLIB "Build ZLib as an external project" OFF)
add_library(zliblib STATIC IMPORTED GLOBAL)
if(BUILD_ZLIB)
ExternalProject_Add(
  zlib
  PREFIX 3rdparty
  URL https://github.com/madler/zlib/archive/refs/tags/v1.2.11.zip
  BUILD_ALWAYS ON
  INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty/"
  BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/3rdparty/lib/libzlibstatic.a"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)

  set_property(
    TARGET zliblib PROPERTY IMPORTED_LOCATION
                             "${CMAKE_BINARY_DIR}/3rdparty/lib/libzlibstatic.a")
else(BUILD_ZLIB)
  find_package(ZLIB REQUIRED)
  set_property(
    TARGET zliblib PROPERTY IMPORTED_LOCATION
                             "${ZLIB_LIBRARIES}")
endif(BUILD_ZLIB)


# Build SZIP as an external project
ExternalProject_Add(
  szip
  PREFIX 3rdparty
  URL https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
  BUILD_ALWAYS ON
  INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty/"
  BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/3rdparty/lib/libszip-static.a"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  INSTALL_COMMAND  ${CMAKE_COMMAND} --install .
  COMMAND   ${CMAKE_COMMAND} -E remove -f ${CMAKE_BINARY_DIR}/3rdparty/src/szip-build/szip-config.cmake
  COMMAND   ${CMAKE_COMMAND} -E remove -f ${CMAKE_BINARY_DIR}/3rdparty/src/szip-build/szip-targets.cmake)
add_library(sziplib STATIC IMPORTED GLOBAL)
set_property(
  TARGET sziplib PROPERTY IMPORTED_LOCATION
                           "${CMAKE_BINARY_DIR}/3rdparty/lib/libszip-static.a")

# Build HDF5 as an external project
if(WIN32)
  set(HDF_FORTRAN_LINK_FLAGS "-Wl,--allow-multiple-definition")
endif(WIN32)
ExternalProject_Add(
  HDF5
  PREFIX 3rdparty
  URL https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_10_7.zip
  BUILD_ALWAYS ON
  INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty/"
  BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5.a" "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5_fortran.a" "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5_f90cstub.a"
  CMAKE_ARGS -DHDF5_BUILD_FORTRAN:bool=True
             -DHDF5_INSTALL_MOD_FORTRAN:string=STATIC
             -DHDF5_ENABLE_SZIP_SUPPORT:bool=True
             -DHDF5_ENABLE_Z_LIB_SUPPORT:bool=True
	     -DBUILD_TESTING:bool=False
             -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
	     -DCMAKE_Fortran_FLAGS:string=${HDF_FORTRAN_LINK_FLAGS}
             -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
add_library(hdf5main STATIC IMPORTED GLOBAL)
set_property(
  TARGET hdf5main
  PROPERTY IMPORTED_LOCATION
           "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5.a")
add_dependencies(hdf5main HDF5 szip)
add_library(hdf5fortran STATIC IMPORTED GLOBAL)
set_property(
  TARGET hdf5fortran
  PROPERTY IMPORTED_LOCATION
           "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5_fortran.a")
add_dependencies(hdf5fortran HDF5 szip)
add_library(hdf5f90cstub STATIC IMPORTED GLOBAL)
set_property(
  TARGET hdf5f90cstub
  PROPERTY IMPORTED_LOCATION
           "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5_f90cstub.a")
add_dependencies(hdf5f90cstub HDF5 szip)

# Set necessary include directories
include_directories(${CMAKE_BINARY_DIR}/3rdparty/include
                    ${CMAKE_BINARY_DIR}/src/neutron/libGudrunN)

# Define installation directory
get_filename_component(default_prefix "./" ABSOLUTE)
set(CMAKE_INSTALL_PREFIX
    ${default_prefix}
    CACHE
      STRING
      "Choose the installation directory (a bin/ directory will be created there)."
      FORCE)

# Define compile flags and linker lines for builds
set(CMAKE_Fortran_FLAGS_RELEASE " -funroll-all-loops -fno-f2c -O2")

add_subdirectory(src)
