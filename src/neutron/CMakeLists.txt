# Specify subdirectories to descend into and build
add_subdirectory(libget)
add_subdirectory(libGudrunN)

# Build GudrunN
add_executable(gudrun_dcs gudrun_dcs.f90)
target_link_libraries(
  gudrun_dcs
  PUBLIC gudrun get
  PRIVATE hdf5libs)

# Build purge_det
add_executable(purge_det purge_det.f90)
target_link_libraries(
  purge_det
  PUBLIC gudrun get
  PRIVATE hdf5libs)

# Install targets
install(TARGETS gudrun_dcs purge_det RUNTIME DESTINATION bin)
