!     
! File:   ran1_routines.f90
! Author: aks45
!
! Created on 23 October 2013, 16:18
!
MODULE ran1_routines

    implicit none
    
	INTEGER IA,IM,IQ,IR,NTAB,NDIV
	REAL AM,EPS,RNMX	
	PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	INTEGER mcall,idum,iy,iv(NTAB)

    
CONTAINS

        subroutine initran1(nchan)
	integer i,nchan
!c
!c open file containing random number parameters, or else initialise
!c
	open(nchan,file='ran1.dat',status='old',err=200)
		read(nchan,*,end=201,err=201) mcall,idum,iy,(iv(i),i=1,ntab)
!c
!c if the number of calls to ran1 exceeds 1.0E8, the sequence is reset, as this
!c is approximately 1/20th of the period of the ran1 generator
!c
		if(mcall.gt.100000000) then
			idum=0
			iy=0
			mcall=0
		endif
		close(nchan)
		return
201	close(nchan)
200	continue
!c
!!c if error or end of file occurs, initialise values
!c
		idum=0
		iy=0
		mcall=0
		return
	end subroutine initran1
        
	subroutine saveran1(nchan)
	integer i,nchan
!c
!c open file containing random number parameters, or else initialise
!c
	open(nchan,file='ran1.dat',status='unknown')
		write(nchan,*) mcall,idum,iy,(iv(i),i=1,ntab)
	return
        end subroutine saveran1
        
	real FUNCTION ran1()
            
            integer j,k
!c
!c \Minimal" random number generator of Park and Miller with Bays-Durham shue and
!c added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
!c the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
!c alter idum between successive deviates in a sequence. RNMX should approximate the largest
!c floating value that is less than 1.
!c
!c
!c Initialize if idum.le.0
!c
	if (idum.le.0.or.iy.eq.0) then 
!c
!c Be sure to prevent idum = 0.
!c
		idum=max(-idum,1) 
!c
!c Load the shue table (after 8 warm-ups).
!c
		do j=NTAB+8,1,-1 
			k=idum/IQ
			idum=IA*(idum-k*IQ)-IR*k
			if (idum.lt.0) idum=idum+IM
			if (j.le.NTAB) iv(j)=idum
		enddo
		iy=iv(1)
	endif
!c
!c Start here when not initializing.
!c
	k=idum/IQ 
!c
!c Compute idum=mod(IA*idum,IM) without overflows by Schrage's method.
!c
 	idum=IA*(idum-k*IQ)-IR*k 
	if (idum.lt.0) idum=idum+IM
!c
!c Will be in the range 1:NTAB.
!c
	j=1+iy/NDIV 
	iy=iv(j) 
!c
!c Output previously stored value and rell the shue table.
!c
	iv(j)=idum
!c
!c Because users don't expect endpoint values.
!c
	ran1=min(AM*iy,RNMX) 
	mcall=mcall+1
	return
        END function ran1

end module ran1_routines