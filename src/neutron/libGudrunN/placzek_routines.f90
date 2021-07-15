!     
! File:   placzek_routines.f90
! Author: aks45
!
! Created on 15 November 2013, 14:36
!

MODULE pla_corr
    
    implicit none
    
    integer                                 :: nwavpla		!number of wavelengths for Placzek
    integer                                 :: nangpla		!number of angles for Placzek
    real, dimension(:), allocatable         :: wavpla!(mcorrwav)		!wavelengths for Placzek
    real, dimension(:), allocatable         :: angpla!(mgroup)		!angles for Placzek
    real, dimension(:), allocatable         :: lenpla!(mgroup)		!secondary flight paths for Placzek
    real, dimension(:,:), allocatable       :: placzek!(mcorrwav,mgroup)	!Placzek correction

    
END MODULE pla_corr
    
MODULE van_placzek
    
    implicit none
    
    integer                                 :: nvwavpla		!number of wavelengths for Placzek
    real, dimension(:), allocatable         :: vwavpla!(mcorrwav)		!wavelengths for Placzek
    real, dimension(:,:), allocatable       :: vplaczek!(mcorrwav,mgroup)	!single scattering Placzek
    
END MODULE van_placzek

MODULE sam_placzek
    
    implicit none
    
    integer                                 :: nswavpla		!number of wavelengths for Placzek
    real, dimension(:), allocatable         :: swavpla!(mcorrwav)		!wavelengths for Placzek
    real, dimension(:,:), allocatable       :: splaczek!(mcorrwav,mgroup)	!single scattering Placzek
    
END MODULE sam_placzek

MODULE placzek_routines
    
    implicit none
    
    CONTAINS

!***********************************************************************************
!*
!*	get_placzek.FOR
!*
!*	A K Soper, May 2001
!*
!*	gets Placzek correction for any run
!*
!***********************************************************************************
    subroutine get_placzek(ntype,run,nperrq,forcecalc)

        use reallocation_routines
        use inputfilestrings
        use run_par
        use calibration_routines
        use pla_corr
        use van_placzek
        use sam_placzek
!c
!c internal variables
!c
        character*256 fname		!dummy file name
        integer il,ib,ip,ierr			!internal indices
        character*256 run			!filename
        integer nperrq		!period number of data in raw file
	integer ntype			!type of file being read
	integer forcecalc		!= 1 forces calculation of corrs.
	real testup,testdn

        ierr=1
        if(forcecalc.eq.0) then
!c
!c set up filename of file to be read
!c
            fname=run
            call change_ext(fname,'pla')
            write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
!c
!c try to open Placzek correction file. If an error occurs programme
!c redirects to calculate the Placzek correction from input parameters 
!c
            open(10,file=fname,status='old',iostat=ierr)
            if(ierr.eq.0) then
                read(10,*,iostat=ierr) nwavpla,nangpla
                if(ierr.eq.0) then
                    call reallocate1d_r(angpla,nangpla)
                    call reallocate1d_r(wavpla,nwavpla)
                    call reallocate2d_r(placzek,nwavpla,nangpla)
                    ip=0
                    do while (ip.lt.nangpla.and.ierr.eq.0)
                        read(10,*,iostat=ierr) ib,angpla(ip)
                        do il=1,nwavpla
                            read(10,*,iostat=ierr) wavpla(il),placzek(il,ip)
                        end do
                    end do
                end if
            end if
            close(10)
        end if
        if(ierr.ne.0) then
!c
!c Calculate placzek correction correction from input sample parameters
!c
            call calc_placzek(ntype,run,nperrq)
        end if
!c
!c save the results in the appropriate array
!c
	if(ntype.eq.1) then
		nvwavpla=nwavpla
                call reallocate1d_r(vwavpla,nwavpla)
                call reallocate2d_r(vplaczek,nwavpla,nangpla)
		do il=1,nwavpla
			vwavpla(il)=wavpla(il)
		end do
		do ib=1,nangpla
			do il=1,nwavpla
				vplaczek(il,ib)=placzek(il,ib)
			end do
		end do
	else
		nswavpla=nwavpla
                call reallocate1d_r(swavpla,nwavpla)
                call reallocate2d_r(splaczek,nwavpla,nangpla)
		do il=1,nwavpla
			swavpla(il)=wavpla(il)
		end do
		do ib=1,nangpla
			do il=1,nwavpla
				splaczek(il,ib)=placzek(il,ib)
			end do
		end do
	endif
        if(allocated(wavpla)) deallocate(wavpla,placzek)
	return
    
    end subroutine get_placzek
    
!***********************************************************************************
!*
!*	calc_Placzek.FOR
!*
!*	A K Soper, May 2001
!*
!*	calculates the Placzek correction for any run
!*
!***********************************************************************************
    subroutine calc_placzek(ntype,run,nperrq)

        use reallocation_routines
        use inputfilestrings
        use run_par
        use beam_routines
        use calibration_routines
        use groups_routines
        use pla_corr
    	use van_par
    	use sam_par
!c
!c internal variables
!c
	character*256                           :: fname		!file name for m.s.corrections
	integer                                 :: ib,il,ie,is,ig,ic,id    		!internal indices
	character*256                           :: run			!filename
	integer                                 :: nperrq		!period number of data in raw file
	integer                                 :: ntype			!type of file being calculated for
	integer                                 :: nelement		!local number of elements
	integer                                 :: ncmin
	real, dimension(:), allocatable         :: yout!(mcorrwav)		!output placzek correction
	real, dimension(:), allocatable         :: ysum!(mcorrwav)		!temporary arrays
	real, dimension(:), allocatable         :: frac!(melement)	!relative weights for Placzek correction
	real, dimension(:), allocatable         :: mass!(melement)		!masses of elements
	real                                    :: fprat,fracsum			!flight path ratio for each group
	real                                    :: temp			!temperature (K) of sample
	real                                    :: aa,bb,cc			!temporary values
	real                                    :: tthmin,tthmax,tthstep,tthstep2
!c
!c set up the relative weighting arrays for each Placzek correction, based
!c on scattering cross sections and atomic fractions
!c
	fracsum=0
	if(ntype.eq.1) then
            nwavpla=nlambv
            nelement=nvelement
            temp=abs(vtemp)
            call reallocate1d_r(mass,nelement)
            call reallocate1d_r(frac,nelement)
            call reallocate1d_r(wavpla,nwavpla)
            do ie=1,nelement
                mass(ie)=vatwt(ie)
		frac(ie)=vscatcs(ie)*vfrac(ie)
		fracsum=fracsum+frac(ie)
            end do
            do ie=1,nelement
		frac(ie)=frac(ie)/fracsum
            end do
            do il=1,nwavpla
		wavpla(il)=vlamb(il)
            end do
	else
            nwavpla=nlambs
            nelement=nselement(1)
            temp=abs(stemp(1))
            call reallocate1d_r(mass,nelement)
            call reallocate1d_r(frac,nelement)
            call reallocate1d_r(wavpla,nwavpla)
            do ie=1,nelement
		mass(ie)=satwt(ie,1)
		frac(ie)=sscatcs(ie,1)*sfrac(ie,1)
		fracsum=fracsum+frac(ie)
            end do
            do ie=1,nelement
		frac(ie)=frac(ie)/fracsum
            end do
            do il=1,nwavpla
		wavpla(il)=slamb(il)
            end do
	endif
        call reallocate1d_r(yout,nwavpla)
        call reallocate1d_r(ysum,nwavpla)

!c Set up the angles and average secondary flight paths for Placzek correction

        if(inst(1:len_trim(inst)).eq.'D4C') then
!For reactor instrument the wavelength is fixed, so we need to get the full range of required angles
            nangpla=180/ndeg
            call reallocate1d_r(lenpla,nangpla)
            call reallocate1d_r(angpla,nangpla)
            call reallocate2d_r(placzek,nwavpla,nangpla)
            do ic=1,nangpla
                angpla(ic)=(real(ic)-0.5)*real(ndeg)
                lenpla(ic)=0.0
            end do
        else
            nangpla=ngroup
            call reallocate1d_r(lenpla,nangpla)
            call reallocate1d_r(angpla,nangpla)
            call reallocate2d_r(placzek,nwavpla,nangpla)
!c Set the Placzek angles and secondary flight path arrays
        	do ic=1,nangpla
	            lenpla(ic)=lengrp(ic)
	            angpla(ic)=tthgrp(ic)
	        end do
        end if
!c Calculate Placzek correction
	do ib=1,nangpla
!c
!c calculate the appropriate corrections
!c
            do il=1,nwavpla
            	ysum(il)=0.0
            end do
            fprat=lenpla(ib)/lenin
            if(angpla(ib).gt.0.0.and.temp.gt.0.0) then
		do ie=1,nelement
                    call tofideal(3,nwavpla,wavpla,fprat,angpla(ib),temp,mass(ie),yout)
		    do il=1,nwavpla
                        ysum(il)=ysum(il)+frac(ie)*yout(il)
                    end do
                end do
            endif
            do il=1,nwavpla
		placzek(il,ib)=ysum(il)
            end do
	end do
!c
!c set up filename of file to write
!c
        fname=run
        call change_ext(fname,'pla')
        write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
!c
!c Open the Placzek correction file and write it. 
!c
	open(10,file=fname,status='unknown')
	write(10,*) nwavpla,nangpla
	do ib=1,nangpla
            write(10,*) ib,angpla(ib),lenpla(ib)
            do il=1,nwavpla
		aa=wavpla(il)
		bb=placzek(il,ib)
		write(10,102) aa,bb
            end do
	end do
102	format(1x,4(1x,e12.5))
	close(10)
	return

    end subroutine calc_placzek

!c
!c Subroutine TOFIDEAL
!c
!c This program integrates the ideal gas law for a TOF experiment along a path
!c of constant TOF at a series of elastic input parameters energies, wavevectors
!c wavelengths or Q-values
!c
    subroutine tofideal(ntype,lptin,xin,fprat,thetad,temp,amasse,yout)

        use spectrum_par
        
	integer, parameter :: iymax=300000		!maximum number of trials
	integer, parameter :: nst=10
	integer ntype		!Type of x-scale: 1=E, 2=k, 3=lambda, 4=Q
	integer lptin		!no. of wavelengths to calculate corrections
	real, dimension(lptin), intent(in)          :: xin!(mcorrwav)	!wavelengths. It assume these have already been set up
	real fprat		!Flight path ratio
	real thetad		!Detector scattering angle in degrees
	real temp		!Sample temperature
	real amasse		!Effective mass of the sample
	real, dimension(lptin), intent(out)         :: yout!(mcorrwav)	!Integrated dynamic scattering law
        real lambdaconv        !Conversion factor to wavelength from wavevector
        real pi,pi2,piconv,fprat1,xmin,ymin,sigma0,theta,ctheta,theta2,qconv,sumli,sumlf,q
        real ake,ake2,efstar,efe,spece,sigt,yst,sum,y,y2,d,x,aki,speci,x2,omega,qsq,q0sq,dif0,sig,fact,fact1
        real sum1,dif,dif2,vibfc,sumv
        integer ik,iy
        logical test
 	
        pi=4.0*atan(1.0)
	pi2=2.0*pi
	piconv=pi/180
	fprat1=1.0+fprat
	xmin=1.0/fprat1
	ymin=fprat*xmin
	sigma0=2.0*0.041591*temp
	theta=thetad*piconv
	ctheta=cos(theta)
	theta2=0.5*theta
	qconv=2.0*sin(theta2)
!	if(amasse.lt.3.0) then
!            open(10,file='inelasticity.dat',status='unknown')
!            write(10,*)'tofideal beginning ',lptin, thetad
!        endif
        do ik=1,lptin
!            if (amasse.lt.3.0) write(10,*)ik,xin(ik)
            sumli=0.0
            sumlf=0.0
            q=xin(ik)
!c
!c convert input x-scale to wavevector scale
!c
            if(q.gt.0.0) then
!c
!c ntype=1 means incident values are energy
!
		if(ntype.eq.1) then
                    ake2=q/2.717
                    ake=sqrt(ake)
!c
!c ntype=2 means incident values are wavevectors
!
		else if(ntype.eq.2) then
                    ake=q
                    ake2=ake*ake
!c
!c ntype=3 means incident values are wavelengths
!c
                else if(ntype.eq.3) then
                    ake=pi2/q
                    ake2=ake*ake
!c
!c else assume values are Q values
!c
		else
                    ake=q/qconv
                    ake2=ake*ake
		endif
                lambdaconv=pi2/ake
		efstar=ef0/ake
!c detector efficiency at elastic energy
		efe=1.0-exp(-efstar)
!c incident spectrum at elastic energy
		spece=spectrum(2,ake)
		sigt=2.0*sigma0*(1.0-ctheta)/(ake2*amasse)
		sigt=sqrt(sigt)
		yst=sigt/nst
		sum=0.0
		iy=0
                test=.false.
                do while (.not.test.and.iy.lt.iymax)
                    iy=iy+1
                    y=ymin+iy*yst
                    y2=y*y
                    d=1.0/(fprat1*y-fprat)
                    x=y*d
!c incident spectrum at incident energy
                    aki=x*ake
                    speci=spectrum(2,aki)
                    x2=x*x
                    omega=x2-y2
                    qsq=x2+y2-2.0*x*y*ctheta
                    q0sq=qsq/amasse
!c	write(10,999) y,x,omega,qsq,ake2
                    dif0=omega-q0sq
                    sig=sigma0*q0sq/ake2
                    fact=2.0*y2*x*fprat1*yst*(1.0-exp(-efstar/y))
                    fact1=fact*speci/sqrt(pi2*sig)
                    sum1=0.0
                    dif=dif0
                    dif2=dif*dif/(2.0*sig)
999	format(6(1x,e12.5))
                    vibfc=-dif2
                    if(vibfc.lt.20.0) then
			sumv=fact1*exp(vibfc)
                    else
			sumv=0.0
                    endif
!c	write(10,999) float(iv)
                    sum1=sum1+sumv
!c	if (amasse.lt.3.0) then
!c            write(10,999) lambdaconv/x,lambdaconv/y,omega,dif0,omv,sum1
!c      endif
                    sum=sum+sum1
!c Form average incident and final wavelengths.
                    sumli=sumli+sum1*lambdaconv/x
                    sumlf=sumlf+sum1*lambdaconv/y
                    test=sum1.le.0.0001*sum.and.omega.lt.q0sq
                end do
!                write(6,*) 'tofideal> ',iy,iymax,sum,sum1,omega,q0sq 
                sumli=sumli/sum
                sumlf=sumlf/sum
!			if (amasse.lt.3.0) write(10,102) iy
!			if (amasse.lt.3.0) write(10,998) q,sum/(spece*efe),sumli,sumlf
		yout(ik)=sum/(spece*efe)-1.0
            else
		yout(ik)=0.0
            endif
	end do
!	if (amasse.lt.3.0) close(10)
!c	lptout=iy-1
!102	format(1x,'Number of y steps = ',i7)
!998	format(1x,'X value ',e12.5,' Integral = ',e12.5,' Mean incident and final wavelengths:', e12.5,1x,e12.5)
	return

    end subroutine tofideal

END MODULE placzek_routines
