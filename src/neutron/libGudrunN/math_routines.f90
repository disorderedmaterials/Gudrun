!     
! File:   trig_functions.f90
! Author: aks45
!
! Created on 05 November 2013, 21:25
!

MODULE math_routines
      
    implicit none
      
    real, parameter       :: degradconv=4.0*atan(1.0)/180.0
    
    CONTAINS
      
    function cosd(deg)
        real :: cosd,deg
	cosd=cos(deg*degradconv)
        return
    end

    function sind(deg)
        real :: sind,deg
   	sind=sin(deg*degradconv)
	return
    end

    function acosd(cosphi)
        real :: acosd,cosphi
	acosd=acos(cosphi)
	acosd=acosd/degradconv
        return
    end
    
    function asind(sinphi)
        real :: asind,sinphi
	asind=asin(sinphi)
	asind=asind/degradconv
	return
    end
    
    function atan2d(sinphi,cosphi)
        real :: atan2d,sinphi,cosphi
	atan2d=atan2(sinphi,cosphi)
	atan2d=atan2d/degradconv
	return
    end	

END MODULE math_routines