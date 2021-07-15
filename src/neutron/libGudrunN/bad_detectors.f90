!     
! File:   bad_detector.f90
! Author: aks45
!
! Created on 01 May 2012, 16:09
!

MODULE bad_detectors

    integer                                 :: nbread        !number of spectra read in spec.bad
    integer, dimension(:), allocatable      :: ibad        !bad spectrum index
    integer, dimension(:), allocatable      :: spike        !number of spikes in a spectrum
    integer                                 :: nchfir,nchlas        !range of channels for spike search
    integer                                 :: flagsp        !non-zero to treat spikes and steps
    real                                    :: spikedev            !square of deviation for spike checking

END MODULE bad_detectors
