!     
! File:   spec_van.f90
! Author: aks45
!
! Created on 01 November 2013, 13:47
!

MODULE spec_van
    
	integer nfilev		!number of vanadium runs 
	character(len=256), dimension(:), allocatable   :: runv!(mruns)		!run numbers of vanadium files
	integer nperv
        integer, dimension(:), allocatable              :: nsmoov!(mspec)			!period number of vanadium
	real, dimension(:), allocatable              :: smovanmon!(mchan)		!smoothed vanadium monitor
	real, dimension(:), allocatable              :: smovantrans!(mchan)	!vanadium transmission monitor
	real, dimension(:,:), allocatable              :: smovandet!(mchan,mproc)	!smoothed vanadium detectors
	real, dimension(:,:), allocatable              :: errvandet!(mchan,mproc)
	real, dimension(:), allocatable              :: vancor,vancore!(mspec)		!normalised, corrected vanadium data summed over channels
        real, dimension(:), allocatable              :: vantransmission!(mcorrwav)     !transmission of vanadium
    
    
END MODULE spec_van
   
