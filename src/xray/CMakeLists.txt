# Build calc_corrsx
add_executable(
  calc_corrsx_in_out calc_corrsx_in_out.for cylmulx_azi.for fltmulx_azi.for
                     abs_mul_routines.for auxilliary_routines.for)

# Build tophatsub
add_executable(tophatsub tophatsub.for auxilliary_routines.for)

# Install targets
install(TARGETS calc_corrsx_in_out tophatsub RUNTIME DESTINATION bin)
