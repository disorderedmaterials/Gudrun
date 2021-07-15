!     
! File:   spectrum.f90
! Author: aks45
!
! Created on 01 November 2013, 11:35
!

MODULE spectrum_par
    
        implicit none
        
        integer npar
	integer nrand
	real par(8),parnew(8),parl(8),paru(8)
	real sqrtk,sqrtl
	real ef0,phiepi,et,phimax,w1,w2,alpha,beta
	character(len=256) specname	!name of file containing the spectrum parameters
    
        CONTAINS
        
!c
!c read spectrum parameters
!c
	subroutine get_spectrum_par

            use run_par
            
            character(len=256) text
            integer ivar,i
            open(10,file=specname,status='old',iostat=ivar)
            if(ivar.ne.0) then
                specname=convertfilename(specname)
                open(10,file=specname,status='old',iostat=ivar)
            endif
            if(ivar.ne.0) then
		write(6,200) specname(1:80)
200   format(/'get_spectrum_par> Specified spectrum parameter file: ' &
      ,a,' does not exist')
		stop
            else
!c
!c read parameters
!
		npar=8
		do i=1,npar
			read(10,*) par(i),parl(i),paru(i)
			if(par(i).gt.paru(i)) par(i)=paru(i)
			if(par(i).lt.parl(i)) par(i)=parl(i)
			parnew(i)=par(i)
		end do
		ef0=parnew(1)
		phiepi=parnew(2)
		et=parnew(3)
		phimax=parnew(4)
		w1=parnew(5)
		w2=parnew(6)
		alpha=parnew(7)
		beta=parnew(8)
		sqrtl=sqrt(0.081787)
		sqrtk=sqrt(0.0020717)
            endif
            close(10)
            return
	end

        function spectrum(ntype,val)
            
            real spectrum,val,del,en,ensqrt,expon,rat,sum
            integer ntype,i,ivar


!!c function to calculate incident spectrum as a function of energy (1), 
! wavevector (2) or wavelength (3)
!c  - ntype=1: energy in mev
!c  - ntype=2: wavevector in A**-1
!c  - ntype=3: wavelength in A
!c
!c
!c ntype = 1 means input value is an energy (default value)
!c
            if(val.gt.0.0) then
!c
!c ntype = 2 means input value is a wavevector
!c
		if(ntype.eq.2) then
			ensqrt=sqrtk*val
			en=ensqrt*ensqrt
!c
!c ntype = 3 means input value is a wavelength
!c
		else if(ntype.eq.3) then
			ensqrt=sqrtl/val
			en=ensqrt*ensqrt
		else
			en=0.001*val
			ensqrt=sqrt(en)
		endif
!c
!c epithermal function and window function
!c
		expon=(w1/ensqrt-w2)
		if(expon.lt.20.) then
			DEL=1.+EXP(expon)
			DEL=PHIEPI/(DEL*en**alpha)
		else
			del=0.0
		endif
!c
!c maxwellian term
!c
		rat=phimax*en/(et*et)
		RAT=alog(rat)-en/et
		if(rat.gt.-20.) then
			SUM=EXP(RAT)+DEL
		else
			sum=del
		endif
!c
!c convert to per bin
!c
		spectrum=sum*en/val
            else
!c
!c if no input value is specified, then output a uniform spectrumc
!c
		spectrum=1.0
            endif
            return
	end        

END MODULE
