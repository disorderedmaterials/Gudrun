!     
! File:   amutshield_routines.f90
! Author: aks45
!
! Created on 04 November 2013, 12:51
!

MODULE amutshield_routines
    
    implicit none
    
    CONTAINS
    
	subroutine get_amutshield(amut)
        use run_par
        use reallocation_routines
        use groups_routines
    
!c
!c Gets the list of shielding total attenuation coefficients for each group
!c
      integer ierr,ig
      real amut,value

!c Initialise the array in case no file exists

      call reallocate1d_r(amutshield,ngroup)
      do ig=1,ngroup
         amutshield(ig)=abs(amut)
      end do

!c Only read the values if amut > 0

      write(6,*) 'get_amutshield> ',amut
      if(amut.gt.0.0) then

         write(6,*) 'get_amutshield> ',amut
	   open(10,file='amutshield.dat',status='unknown',iostat=ierr)
	   do while (ierr.eq.0)
		read(10,*,iostat=ierr) ig,value
!c We only import a value if the specified group exists and the value is non-zero
            if(ierr.eq.0.and.ig.gt.0.and.ig.le.ngroup.and.value.gt.0.0) amutshield(ig)=value
	   end do
	   close(10)

      end if

	return
	end

	subroutine save_amutshield(amut)
        use run_par
        use groups_routines
    
!c
!c Saves the list of shielding total attenuation coefficients for each group
!c

      real amut
      integer ig,ierr

!c Only save the data if the input attenuation coefficient is > 0

      write(6,*) 'save_amutshield> ',amut
      open(10,file='amutshield.dat',status='unknown',iostat=ierr)
      ig=0
      do while (ig.lt.ngroup.and.ierr.eq.0)
         ig=ig+1
	   write(10,100,iostat=ierr) ig,amutshield(ig)
100   format(1x,i3,1x,f8.3)
	end do
      close(10)


	return
	end 
        
END MODULE amutshield_routines
