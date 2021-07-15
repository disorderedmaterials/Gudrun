!     
! File:   spec_baks.f90
! Author: aks45
!
! Created on 01 November 2013, 13:55
!

MODULE spec_baks
    
	integer nfilebs		!number of vanadium runs 
	character(len=256), dimension(:), allocatable   :: runbs!(mruns)		!run numbers of vanadium files
	integer nperbs
	real, dimension(:), allocatable              :: smobsammon!(mchan)		!smoothed vanadium monitor
	real, dimension(:), allocatable              :: smobsamtrans!(mchan)	!vanadium transmission monitor
	real, dimension(:,:), allocatable              :: normbakdet!(mchan,mproc)	!normalised background detectors
	real, dimension(:,:), allocatable              :: errbakdet!(mchan,mproc) 	!error values for background
    
    
END MODULE
