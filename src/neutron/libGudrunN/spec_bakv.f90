!     
! File:   spec_vanb.f90
! Author: aks45
!
! Created on 01 November 2013, 13:51
!

MODULE spec_bakv
    
	integer nfilebv		!number of vanadium runs 
	character(len=256), dimension(:), allocatable   :: runbv!(mruns)		!run numbers of vanadium files
	integer nperbv
	real, dimension(:), allocatable              :: smobvanmon!(mchan)		!smoothed vanadium monitor
	real, dimension(:), allocatable              :: smobvantrans!(mchan)	!vanadium transmission monitor
    
    
END MODULE

