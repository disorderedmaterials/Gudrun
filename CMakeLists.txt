# CMake project file for gudrun

cmake_minimum_required (VERSION 3.0)
project (gudrun)
enable_language (Fortran)

# Build release version
set(CMAKE_BUILD_TYPE "Release")

# Options to allow setting of paths to HDF and SZIP
SET(HDF5_DIR "C:/Users/axr61803/src/hdf5_build" CACHE PATH "Directory containing HDF5")

# Setup HDF (required)
find_package(HDF5 REQUIRED)

# Set search locations for header files (for both project and external dependencies)
# Two HDF5 directories required, as one contains the Fortran module file (HDF5_INCLUDE_DIR_FORTRAN)
if (WIN32 OR APPLE)
  include_directories(
	${CMAKE_BINARY_DIR}
	${CMAKE_BINARY_DIR}/src/neutron/libGudrunN
	${HDF5_FORTRAN_INCLUDE_DIRS}
	${HDF5_INCLUDE_DIR_FORTRAN}
  )
else(WIN32 OR APPLE)
  include_directories(
	${CMAKE_BINARY_DIR}
	${CMAKE_BINARY_DIR}/src/neutron/libGudrunN
	${HDF5_FORTRAN_INCLUDE_DIRS}
	${HDF5_INCLUDE_DIR_FORTRAN}
	/usr/include
  )
endif(WIN32 OR APPLE)

# Define installation directory
get_filename_component (default_prefix "./" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
	"Choose the installation directory (a bin/ directory will be created there)."
      FORCE)

# Define compile flags and linker lines for builds
set (CMAKE_Fortran_FLAGS_RELEASE " -funroll-all-loops -fno-f2c -O2")
SET(LINK_LIBS "${HDF5_FORTRAN_LIBRARIES}" CACHE STRING "Link-time libraries")

add_subdirectory(src)
