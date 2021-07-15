!     
! File:   histogram.f90
! Author: aks45
!
! Created on 27 April 2012, 13:09
!

MODULE histogram

	integer                         :: nhist,mhist			!dimension of histogram arrays
	real, dimension(:), allocatable :: rathist		!histogram of ratios
	real, dimension(:), allocatable :: errhist		!histogram of ratio errors

END MODULE histogram
