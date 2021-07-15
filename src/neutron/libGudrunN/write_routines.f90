!     
! File:   write_routines.f90
! Author: aks45
!
! Created on 02 November 2013, 14:17
!

MODULE write_routines
    
    implicit none
    
    CONTAINS
    
    subroutine w_diag_file(fname,nx,x,y,e)
	
        character(len=256) fname
	integer ic,nx,ierr
	real x(*),y(*),e(*)
	real rat,err
        
!        write(6,*) 'w_diag_file> ',fname(1:len_trim(fname)),nx,x(nx),y(nx),e(nx)
	open(10,file=fname,status='unknown',iostat=ierr)
	if(ierr.eq.0) then
            write(10,100) '#',fname(1:len_trim(fname))
100         format(a1,1x,a)
            write(10,100) '#'
            write(10,100) '#'
            write(10,100) '#'
            do ic=1,nx
		rat=y(ic)
		if(e(ic).gt.0.0) then
		   err=sqrt(e(ic))
		else
		   err=0.0
		endif
		write(10,101) x(ic),rat,err
!		write(6,101) x(ic),rat,err
            end do
            close (10)
        else
            write(6,*) 'w_diag_file> Open error for: ',fname(1:len_trim(fname)), nx
        end if
101	format(1x,3(1x,e13.6))
	return
    end
        
END MODULE write_routines