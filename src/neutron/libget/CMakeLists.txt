# Additional defines for building libget under Windows
if(WIN32)
	add_definitions(-D_WIN32 -DGCC_MINGW)
endif(WIN32)

add_library(get
  error_add_dummy.c
  float_convert.c
  getsect.c
  int_convert.c

  float_convert.h
  geniedefs.h
  int_convert.h
  libget.h

  ferror_add_dummy.f
  io.f
)

include_directories(
../
${CMAKE_SOURCE_DIR}
${CMAKE_BINARY_DIR}
${CMAKE_BINARY_DIR}/src
)
