!     
! File:   runfactor_list.f90
! Author: aks45
!
! Created on 23 October 2013, 17:47
!

MODULE runfactor_list
    
   implicit none

!c
!c Common blocks containing the factors for each run. These are applied BEFORE ANY processing is 
!c performed.
!c
    integer no_of_runno,no_of_runno_dim		!number of runs requiring non-unity factors
    character(len=256), dimension(:), allocatable ::     runno_list  !List of run numbers requiring non-unity factors
    real, dimension(:), allocatable ::              runno_factor	!List of factors
    
    CONTAINS
    
    subroutine get_runfactor_list
            
        use reallocation_routines

        integer:: ierr
        character(len=256) list
        real::  factor
!
! Gets a list of factors for specified runs. Run numbers which are not specified are automatically
! given a factor of 1.0. If the file runfactor_list.dat does not exist, then all the runs 
! automatically will have a factor of 1.0.
!
        no_of_runno=0
        open(10,file='runfactor_list.dat',status='unknown')
        no_of_runno_dim=0
        read(10,*,iostat=ierr) list,factor
        do while (ierr.eq.0)
            no_of_runno=no_of_runno+1
            do while (no_of_runno.gt.no_of_runno_dim) 
                no_of_runno_dim=no_of_runno_dim+5
                call reallocate1d_c(runno_list,len(runno_list),no_of_runno_dim)
                call reallocate1d_r(runno_factor,no_of_runno_dim)
            end do
            runno_list(no_of_runno) = list
            runno_factor(no_of_runno) = factor
            read(10,*,iostat=ierr) list,factor
        end do
        close(10)
        return
    end subroutine get_runfactor_list

    function get_runfactor(runno)
        
        character(len=256) runno
        integer icheck,i
!c
!c Finds the factor associated with a particular filename. If not on the list the factor
!c for this run number is assumed to be unity.
!c
	real runfac,get_runfactor
	runfac=1.0
	if(no_of_runno.gt.0) then
            icheck=0
            i=1
            do while (icheck.eq.0.and.i.le.no_of_runno)
            	if(runno.eq.runno_list(i)) then
                    runfac=runno_factor(i)
                    icheck=1
		endif
		i=i+1
            end do
	endif	
	get_runfactor=runfac
	return
        
    end function get_runfactor
    
END MODULE
