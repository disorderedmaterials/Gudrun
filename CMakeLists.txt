# CMake project file for gudrun

cmake_minimum_required(VERSION 3.1)
project(gudrun)
enable_language(Fortran)

# Build release version
set(CMAKE_BUILD_TYPE "Release")

include(ExternalProject)

# Build SZIP as an external project
ExternalProject_Add(
  szip
  PREFIX 3rdparty
  URL https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
  BUILD_ALWAYS ON
  INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty/"
  BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/3rdparty/lib/libszip-shared.so"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
add_library(szip_lib SHARED IMPORTED GLOBAL)
set_property(
  TARGET szip_lib PROPERTY IMPORTED_LOCATION
                           "${CMAKE_BINARY_DIR}/3rdparty/lib/libszip-shared.so")

# Build HDF5 as an external project
ExternalProject_Add(
  HDF5
  PREFIX 3rdparty
  URL https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_10_7.zip
  BUILD_ALWAYS ON
  INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty/"
  BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5_fortran.so"
  CMAKE_ARGS -DHDF5_BUILD_FORTRAN:bool=True
             -DHDF5_INSTALL_MOD_FORTRAN:string=SHARED
             -DHDF5_ENABLE_SZIP_SUPPORT:bool=True
             -DHDF5_ENABLE_Z_LIB_SUPPORT:bool=True
             -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
add_library(hdf5fortran SHARED IMPORTED GLOBAL)
set_property(
  TARGET hdf5fortran
  PROPERTY IMPORTED_LOCATION
           "${CMAKE_BINARY_DIR}/3rdparty/lib/libhdf5_fortran.so")
add_dependencies(hdf5fortran HDF5 szip)

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
