!     
! File:   spec_sam.f90
! Author: aks45
!
! Created on 01 November 2013, 14:03
!

MODULE spec_sam
    
    	integer, dimension(:), allocatable              :: nfiles!mcont)		!number of sample runs 
	character(len=256), dimension(:,:), allocatable :: runs!(mruns,mcont)	!run numbers of sample files
	integer, dimension(:), allocatable              :: npers!(mcont)		!period number for sample files
	real, dimension(:,:), allocatable               :: smosammon!(mchan,mcont)		!smoothed sample monitor
	real, dimension(:,:), allocatable               :: smosamtrans!(mchan,mcont)		!sample transmission monitor
	real, dimension(:,:,:), allocatable             :: normsamdet!(mchan,mproc,mcont)	!normalised sample detectors
	real, dimension(:,:), allocatable               :: normsamdetnosub!(mchan,mproc)	!normalised sample detectors but no self subtraction
	real, dimension(:,:,:), allocatable             :: errsamdet!(mchan,mproc,mcont) 	!error values for sample results
	real, dimension(:), allocatable                 :: speccheck!(mspec)		!spectrum numbers
	real, dimension(:), allocatable                 :: checksum,errchecksum!(mspec) !Stores check sum and error
	real, dimension(:,:), allocatable               :: samcor,samcore!(mspec,mcont)	!Normalised, corrected sample data summed over channels
        real, dimension(:,:), allocatable               :: samtransmission!(mcorrwav,mcont)     !Transmission of sample and containers

END MODULE
     
