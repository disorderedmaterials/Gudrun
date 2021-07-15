!     
! File:   spec_ratio.f90
! Author: aks45
!
! Created on 27 April 2012, 13:18
!

MODULE spec_ratio

	real, dimension(:,:), allocatable   :: specrat		!integrated spectrum ratios
	real, dimension(:,:), allocatable   :: errrat		!rms deviation on ratio

END MODULE spec_ratio
