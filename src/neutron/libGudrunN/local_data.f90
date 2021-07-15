!     
! File:   local_data.f90
! Author: aks45
!
! Created on 30 April 2012, 13:56
!

MODULE local_data
    
    use reallocation_routines

    implicit none
    
! temporary values common to routines and modules
    integer                                 :: nspecf                !first spectrum to read
    integer                                 :: nspecl                !last spectrum to read
    integer                                 :: nspecrq                !no. of spectra requested
    integer                                 :: nchanrq                !no. of channels requested
    integer                                 :: nf,nl                        !segment of c_DATA to retrieve
    integer                                 :: mcount           !size of array to be returned by getdat
    integer, dimension(:,:,:), allocatable  :: counts        !temporary array of counts from a spectrum
    integer, dimension(2)                   :: oldtcountsdims !Dimensions of the counts array
    integer, dimension(2)                   :: tcountsdims !Dimensions of the counts array
    integer, dimension(3)                   :: oldcountsdims !Dimensions of the counts array
    integer, dimension(3)                   :: countsdims !Dimensions of the counts array
    integer, dimension(3)                   :: countsoffsets  !Offsets for counts array - used for nxs files.
    integer                                 :: nspecproc                !no. of spectra to process
    integer, dimension(:), allocatable      :: specproc        !spectrum numbers to process
    integer                                 :: nspecprocgood                !no. of spectra to process after vanadium smooth
    integer, dimension(:), allocatable      :: specprocgood        !spectrum numbers to process after vanadium smooth
    real, dimension(:,:), allocatable       :: tcounts                !temporary store of sums of counts
    real, dimension(:), allocatable         :: detcount                !normalised data from a detector
    real, dimension(:), allocatable         :: errcount                !normalised errors from a detector
    real, dimension(:), allocatable         :: wavebin                !wavelength bin values
    real, dimension(:), allocatable         :: wavebound                !wavelentth boundary values
    real, dimension(:), allocatable         :: wavecorr          !wavelength values of corrections
    real, dimension(:), allocatable         :: corrbin                !interpolated correction
    real, dimension(:), allocatable         :: corrtemp        !correction at calculated wavelength scale
    real, dimension(:), allocatable         :: wavetemp,dettemp,errtemp !temporary arrays
    
    contains
    
    subroutine check_and_allocate_local_data(newsize)
        !Checks dimensions and allocates space for temporary arrays, but only if newsize > currentsize
        integer newsize
        
        call reallocate1d_r(detcount,newsize)
        call reallocate1d_r(errcount,newsize)
        call reallocate1d_r(wavebin,newsize)
        call reallocate1d_r(wavebound,newsize)
        call reallocate1d_r(wavecorr,newsize)
        call reallocate1d_r(corrbin,newsize)
        call reallocate1d_r(corrtemp,newsize)
        call reallocate1d_r(wavetemp,newsize)
        call reallocate1d_r(dettemp,newsize)
        call reallocate1d_r(errtemp,newsize)
        return
    end subroutine check_and_allocate_local_data

END MODULE local_data
