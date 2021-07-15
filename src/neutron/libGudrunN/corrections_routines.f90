!     
! File:   corrections_routines.f90
! Author: aks45
!
! Created on 11 November 2013, 11:14
!

MODULE corrections_routines
    
    !Routines in this collection are mostly standalone - do not make reference to other modules.
    !An exception is the first one!
     
    implicit none
    
CONTAINS

!***********************************************************************************
!*
!*	set_corr_ang.FOR
!*
!*	A K Soper, February 2003
!*
!*	Generates theta and phi values for specified sample geometry. 
!*
!*	For cylinders the z-axis runs along the axis of the cylinder and is assumed
!* 	to be at right angles to the incident beam.
!*	
!*	For flat plates the z-axis is assumed to be perpendicular to the surface of 
!*	slab, with the positive axis on the side of the transmitted beam
!*
!***********************************************************************************
    subroutine set_corr_ang(ngeom,sang,ang,azi,nang,nazi)
        
        use run_par
        use reallocation_routines
        use beam_routines
        use calibration_routines
        use groups_routines
        use bad_detectors
!c
!c internal variables
!c
	integer, intent(in)                                 :: ngeom		!1 = cylinder geometry, 2 = flat plate
	integer, intent(out)                                :: nang,nazi		!no. of theta and phi values for corrections
	integer i,ib,j,index,il,k	!internal indices
	integer is,id		!internal indices
	real, intent(in)                                    :: sang			!Sample rotation angle (flat plate only)
	real, dimension(:), allocatable, intent(out)        :: ang!(mcorrang)	!theta values for corrections
	real, dimension(:), allocatable, intent(out)        :: azi!(mcorrang)	!azimuthal angle values for corrections
	real angmax,angmin	!maximum and minimum theta values
	real azimax,azimin	!maximum and minimum azi values
	real angfind,azifind	!theta and phi values corresponding to specified detector
	real angst,azist	!angular step sizes
	real tthetat,phit

        write(6,*) 'set_corr_ang> ',sang
!c
!c step through good detectors and determine the maximum and minimum values
!c of theta and phi for the respective geometries
!c
	angmax=0.0
	angmin=90.0
	azimax=0.0
	azimin=180.0
	do is=1,nspec
!c
!c check this detector is good
!c
            if(ibad(is).eq.0.and.igrp(is).ge.1.and.igrp(is).le.ngroup) then
!c
!c get detector number for this detector
!c
		id=detno(is)
!c
!c get the corresponding correction angles for this detector
!c
		tthetat=ttheta(id)
		phit=phi(id)
		call get_corr_ang(ngeom,sang,tthetat,phit,angfind,azifind)
		if(angfind.gt.angmax) angmax=angfind
		if(angfind.lt.angmin) angmin=angfind
		if(azifind.gt.azimax) azimax=azifind
		if(azifind.lt.azimin) azimin=azifind
            endif
	end do
!c ensure the case of 2-theta=0 is included in the calculation
        call get_corr_ang(ngeom,sang,0.0,0.0,angfind,azifind)
	if(angfind.gt.angmax) angmax=angfind
	if(angfind.lt.angmin) angmin=angfind
	if(azifind.gt.azimax) azimax=azifind
	if(azifind.lt.azimin) azimin=azifind

	write(6,*) 'set_corr_ang> ',angmax,angmin,azimax,azimin
!c
!c specify the corrections angles in steps of 10 degrees
!c
	nang=(angmax-angmin)/ndeg+1
	nazi=(azimax-azimin)/ndeg+1
	if(nang.gt.1) then
		angst=(angmax-angmin)/(nang-1)
	else
		angst=0.0
	endif
	if(nazi.gt.1) then
		azist=(azimax-azimin)/(nazi-1)
	else
		azist=0.0
	endif
        call reallocate1d_r(ang,nang)
	do i=1,nang
		ang(i)=angmin+(i-1)*angst
	end do
        call reallocate1d_r(azi,nazi)
	do i=1,nazi
		azi(i)=azimin+(i-1)*azist
	end do
	write(6,*) nang,(ang(i),i=1,nang)
	write(6,*) nazi,(azi(i),i=1,nazi)
	return

    end subroutine set_corr_ang
    
!***********************************************************************************
!*
!*	get_corr_ang.FOR
!*
!*	A K Soper, February 2003
!*
!*	Generates theta and phi values for specified sample geometry. 
!*
!*	For cylinders the z-axis runs along the axis of the cylinder and is assumed
!* 	to be at right angles to the incident beam.
!*	
!*	For flat plates the z-axis is assumed to be perpendicular to the surface of 
!*	slab, with the positive axis on the side of the transmitted beam
!*
!***********************************************************************************
    subroutine get_corr_ang(ngeom,sangd,tdetd,pdetd,tfind,pfind)
!c
!c A K Soper, February 2003
!c
!c converts specified detector coordinates (z-axis in downstream beam direction), x-axis
!c in horizontal plane left of z-axis, y-axis vertically upwards) to coordinates appropriate
!c to cylinders (z axis vertically upwards, x axis along downstream beam direction, y-axis
!c in horizontal plane to left of beam direction when facing in same direction as beam) or
!c flat plate (z-axis perpendicular to plane of slab on downstream side, y-axis perpendicular
!c this is vertical plane, x-axis perpendicular to z-axis in horizontal plane.
!c
	integer ngeom		!Sample geometry
	real sangd		!sample angle (deg)
	real sangr		!sample angle (rads)
	real tdetd,pdetd	!Detector coordinates in instrument coordinate axes (deg)
 	real tdetr,pdetr	!Detector coordinates in instrument coordinate axes (rad)
	real tfind,pfind	!Detector coordinates in sample coordinates.
	real piconv		!converts angles to radians
	real costs,xs,ys	!temporary values
        real x,y,z,xprime,yprime,zprime
!c	write(6,*) ngeom,sangd,tdetd,pdetd
	piconv=4.0*atan(1.0)/180.0
!c
!c convert to radians
!c
	sangr=sangd*piconv
	tdetr=tdetd*piconv
	pdetr=pdetd*piconv
!c
!c Determine sample geometry
!c
	if(ngeom.eq.1) then
!c
!c cylindrical geometry - corrections will be symmetric about tfind=90
!c
            costs=abs(sin(tdetr)*sin(pdetr))
            if(abs(costs).gt.1.0) costs=sign(1.0,costs)
            tfind=acos(costs)/piconv
!c
!c only need absolute value of pfind, since corrections will be symmetrical
!c about pfind=0
!c
            xs=sin(tdetr)*cos(pdetr)
            ys=cos(tdetr)
            pfind=abs(atan2(xs,ys)/piconv)	

	else if(ngeom.eq.2) then
!c
!c flat plate geometry, sangr specifies the rotation about the (normally) vertical detector y axis.
!c (If axis of rotation is different from this, then the detectors will have to be rotated
!c until their y-axis coincides with the positive axis of rotation.)
!c
            x=sin(tdetr)*cos(pdetr)
            y=sin(tdetr)*sin(pdetr)
            z=cos(tdetr)
!c Rotated coordinates:
            xprime=x*cos(sangr)-z*sin(sangr)
            yprime=y
            zprime=z*cos(sangr)+x*sin(sangr)        
            if(abs(zprime).gt.1.0) zprime=sign(1.0,zprime)
            tfind=acos(zprime)/piconv
!c 
!c Get the phi value in the rotated coordinate geometry. Remember slab has been rotated about the instrument
!c y-axis (normally the vertical axis).
!c
!c However if the slab rotation is zero, there is no need to calculate all the phi values since their corrections
!c will be the same for all phi. (This is provided sample is a true powder or isotropic liquid - if it has preferred orientation
!c or is a single crystal this will not apply.)
!23/6/2014.No there isn't - this is taken care of in do_abs_corr and do_mul_corr
!For flat plates the corrections are always independent of phi, relative
!to rotated z-axis. The change will occur in the range of theta values that
! are sampled.
!            if(sangr.eq.0.0) then
               pfind=0.0
!            else
!                pfind=atan2(yprime,xprime)/piconv
!            endif

	else
            write(6,*) 'Specified geometry ',ngeom,' is not supported'
            stop
	endif
	return
    end subroutine get_corr_ang

    FUNCTION DIST(R1,RADius,OMEGA)

        real r1,radius,omega,dist
        real b,c,d,t
        
        DIST=0.
        B=R1*SIN(OMEGA)
        IF(ABS(B).GT.RADius) RETURN
        T=R1*COS(OMEGA)
        C=RADius*RADius-B*B
        D=SQRT(C)
        IF(R1.le.RADius) then
            DIST=T+D
            RETURN
        else
            DIST=D*(1+SIGN(1.,T))
        end if
        RETURN
    END function dist
    
!C                                                                       MUL14050
!C       *******************                                             MUL14060
!C       *                 *                                             MUL14070
!C       *  FUNCTION DIST1 *                                             MUL14080
!C       *                 *                                             MUL14090
!C       *******************                                             MUL14100
!C                                                                       MUL14110
    FUNCTION DIST1(R1,R2,RADius,Bb,ASIG)

        real dist1,r1,r2,radius,bb,asig,ra,x
        
        DIST1=0
        IF(ABS(Bb).GT.RADius) RETURN
        RA=R1
        IF(RA.GT.RADius) RA=RADius
        X=RA*RA-Bb*Bb
        if(x.lt.0.) x=0.
        DIST1=DIST1+SQRT(X)
        RA=R2
        IF(RA.GT.RADius) RA=RADius
        X=RA*RA-Bb*Bb
        if(x.lt.0.) x=0.
        DIST1=DIST1+SQRT(X)*ASIG
        DIST1=ABS(DIST1)
        RETURN
    END function dist1    
!C                                                                       MUL05400
!C        ********************                                           MUL05410
!C        *                  *                                           MUL05420
!C        *  FUNCTION DIST2  *                                           MUL05430
!C        *                  *                                           MUL05440
!C        ********************                                           MUL05450
!C                                                                       MUL05460
    FUNCTION DIST2(R1,RADius,OMEGA)

        real r1,radius,omega,dist2
        real b,c,d,t
        
!C                                                                       MUL05480
!C CALCULaTES THE MODULUS OF THE DISTANCE OF THE POINT (R,OMEGA)          MUL05490
!C FROM THE RIGHTMOST SURFACE OF A CIRCLE OF RADIUS RAD ALONG A LINE     MUL05500
!C PARALLEL TO OMEGA=0 AXIS.  IF R.GT. RAD THE VALUE RETURNED IS         MUL05510
!C THE LENGTH OF THE SEGMENT WHICH THE LINE CUTS.  IF                    MUL05520
!C  ABS(R*SIN(OMEGA)) IS .GT. RAD, THERE IS NO PATH THROUGH THE          MUL05530
!C CIRCLE AND A VALUE OF ZERO IS RETURNED                                MUL05540
!C                                                                       MUL05550
        DIST2=0
        D=R1*SIN(OMEGA)
        IF(ABS(D).GT.RADius) RETURN
        T=R1*COS(OMEGA)
        DIST2=RADius*RADius-D*D
        DIST2=SQRT(DIST2)
        IF(R1.le.RADius) then
            DIST2=DIST2-T
            RETURN
        else
            DIST2=DIST2-SIGN(DIST2,T)
            RETURN
        end if
    END function dist2
    
!C                                                                       MUL05690
!C      ********************                                             MUL05700
!C      *                  *                                             MUL05710
!C      *  FUNCTION DIST3  *                                             MUL05720
!C      *                  *                                             MUL05730
!C      ********************                                             MUL05740
!C                                                                       MUL05750
    FUNCTION DIST3(R1,RADius,OMEGA)

        real r1,radius,omega,dist3
        real b,c,d,t
        
!C                                                                       MUL05770
!C SAME COMMENTS AS FOR DIST2 EXCEPT THAT THE DISTANCE CALCULATED        MUL05780
!C IS FROM THE POINT (R,OMEG) TO THE LEFT-MOST SURFACE                   MUL05790
!C                                                                       MUL05800
        DIST3=0
        D=R1*SIN(OMEGA)
        IF(ABS(D).GT.RADius) RETURN
        T=R1*COS(OMEGA)
        DIST3=RADius*RADius-D*D
        DIST3=SQRT(DIST3)
        IF(R1.le.RADius) then
            DIST3=DIST3+T
            RETURN 
        else
            DIST3=DIST3+SIGN(DIST3,T)
            RETURN
        end if
    END function dist3
!C                                                                       MUL14260
!C       *******************                                             MUL14270
!C       *                 *                                             MUL14280
!C       *  FUNCTION AVL2  *                                             MUL14290
!C       *                 *                                             MUL14300
!C       *******************                                             MUL14310
!C                                                                       MUL14320
    FUNCTION AVL2(AA,BB,L0)
        
        REAL avl2,L0,AA,BB
        REAL*8 L1,L2,L2R,L2R2,L2R3,L2R4,L2R5 
        real*8 A,B,AL3,AL2,AR,AR2,AR3,AR4,AR5,BR,BR2,BR3, &
        AL2R,AL2R2,AL2R3,AL2R4,AL2R5,TEMP,CONS,  &
        Z1,Z2,Z3,Z4,A1,A2,BX,THRD,C1,C2,RAT,RAT1,RAT2,RAT3,RAT4, &
        E1,E2,COEFF1,COEFF2,COEFF3,DIFF1,DIFF2,SUM,SUM1,SUM2
        
        SUM=0.
        THRD=1./3.
!C                                                                       MUL14420
!C CHECK THAT A IS .GE. B - IF NOT WE SWITCH VALUES                      MUL14430
!C                                                                       MUL14440
        A=AA
        B=BB
        IF(A.lt.B) then
            A=BB
            B=AA
        end if
!C                                                                       MUL14510
!C SELECT LIMITS FOR INTEGRATION                                         MUL14520
!C                                                                       MUL14530
        L1=L0+B
        L2=L0-B
        AL3=0
        IF(L2.LE.0.) AL3=DABS(L2)
!C                                                                       MUL14580
!C WE REDEFINE THE UPPER LIMIT AS L2 AND THE LOWER LIMIT AS AL2.         MUL14590
!                                                                       MUL14600
        AL2=DABS(L2)
        L2=L1
        IF(L0.gt.0.) then
            AR=A/L0
            BR=B/L0
            IF(BR.le.0.01) then
                AVL2=1
                RETURN
            end if
            L2R=L2/L0
            AL2R=AL2/L0
            AR2=AR*AR
            AR3=AR2*AR 
            AR4=AR3*AR
            AR5=AR4*AR
            BR2=BR*BR
            BR3=BR2*BR
            L2R2=L2R*L2R
            L2R3=L2R2*L2R
            L2R4=L2R3*L2R
            L2R5=L2R4*L2R
            AL2R2=AL2R*AL2R
            AL2R3=AL2R2*AL2R
            AL2R4=AL2R3*AL2R
            AL2R5=AL2R4*AL2R
            TEMP=BR2-1
            CONS=9./(8.*BR3)
            CONS=CONS/AR3 
            Z1=AR2*TEMP
            Z2=THRD*(-TEMP-AR2)
            Z3=0.5*AR4
            A1=0.5*(L2R*Z1+AR2*L2R2+Z2*L2R3-0.5*L2R4+0.2*L2R5-Z3)
            A2=0.5*(AL2R*Z1+AR2*AL2R2+Z2*AL2R3-0.5*AL2R4+0.2*AL2R5-Z3)
            BX=-THRD*(-AR3*TEMP+0.2*AR5)
            Z1=-THRD*(-TEMP+0.2*AR2)*AR
            Z2=0.5*AR
            Z3=-0.2*AR
            Z4=0.5*AR3 
            C1=L2R*Z4+L2R2*Z1+L2R3*Z2+L2R4*Z3
            C2=AL2R*Z4+AL2R2*Z1+AL2R3*Z2+AL2R4*Z3
            RAT1=L2/A
            RAT2=AL2/A
            RAT3=1+RAT1
            RAT4=1+RAT2
            E1=A1+BX
            E2=A2+BX
            COEFF1=E1*DLOG(RAT3)-E2*DLOG(RAT4)
            RAT3=1-RAT1
            RAT4=1-RAT2
            RAT3=DABS(RAT3)
            RAT4=DABS(RAT4)
            DIFF1=BX-A1
            DIFF2=BX-A2
            IF(DABS(DIFF2).LT.1.0E-05) RAT4=1.
            IF(DABS(DIFF1).LT.1.0E-05) RAT3=1.
            COEFF2=DIFF1*DLOG(RAT3)-DIFF2*DLOG(RAT4)
            COEFF3=C1-C2
            SUM1=CONS*(COEFF1+COEFF2+COEFF3)
            SUM=SUM+SUM1
        end if
!C COMPUTE INTEGRAL FOR ENCLOSED REGION - IF APLICABLE                   MUL15210
        IF(AL3.gt.0.) then
            RAT=AL3/A
            IF(RAT.GE.1.) then
                SUM2=2.25
            else
                RAT2=RAT*RAT
                RAT3=RAT2*RAT
                DIFF1=1.-RAT2
                DIFF2=(1+RAT)/(1-RAT)
                A1=-0.5*DIFF1**2*DLOG(DIFF2)
                SUM2=(RAT3+RAT+A1)*9.*L0**2*A/(8*B**3)
            end if
            SUM=SUM+SUM2
        end if
        AVL2=SUM
        RETURN 
    END function avl2

    subroutine smoov_limit(nsmoo,nsmoomax,nchan,detcount,errcount,smocount,ratio,ratiolimit,nchfir,nchlas,alimit,ierr,nchanout)
!c
!c does a square-wave smoothing with width proportional to total number of points in array
!c input in detcount, smoothed result in errcount.

!c If the result is .lt. alimit, then the routine returns the value of ierr as true
!c (i.e. it found an error)
!c
!c 14/12/2010 Introduced a check to avoid over smoothing. The sum of std. dev.s is calculated
!c and the sum of differences (raw-smooth)**2 is expected to be no larger than sum of std. dev.s
!c
        integer nsmoo,nsmoomax,nchan,nchfir,nchlas,nchanout
        real, dimension(:)        :: detcount,errcount,smocount
        real  alimit,ratio,ratiolimit

        integer nfirst,nlast,nup,ndown,nchmin,nchmax,ic,nchansmooth,nsmoov,nsum,nbroad,jref,ic1
        real sum,sumstddev,sumsqrdif,stddevdif
        real smocountmax,smocountlimit,add,dif
        integer ierr
!c
!c determine the first and last non-zero bins of input array
!c Determine the range of non-zero points.
        call trim_zeros(1,nchan,nfirst,nlast,errcount)
!        write(6,*) 'smoov_limit> ',nchmin,nchmax,nfirst,nlast
        !To start with zero all the smoothed channels
        do ic=1,nchan
            if(ic.lt.nfirst.or.ic.gt.nlast) then
                smocount(ic)=0.0
                errcount(ic)=0.0
            else
                smocount(ic)=detcount(ic)
            end if
        end do
!c       write(6,*) nchfir,nchlas,nchmin,nchmax,nfirst,nlast
!c
!c smooth the input array with top hat smoothing

        ierr=0
        if(nfirst.eq.0.or.nlast.eq.0) ierr=1
        ! Obviously if there is no data in the array, then no point in smoothing it
        if(ierr.eq.0) then
            !Also no point in smoothing if the number of smooths requested exceeds the maximum number allowed (typically nchan).
            if(nsmoomax.gt.0) then
                call smootophat(nfirst,nlast,nsmoo,nsmoomax,detcount,smocount,errcount,ratio,ratiolimit,nchanout)
            else
                ratio=1.0
            end if
!c Determine the maximum value of smocount.
            smocountmax=0.0
            do ic=nfirst,nlast
                smocountmax=max(smocount(ic),smocountmax)
            end do
            smocountlimit=alimit*smocountmax
!C Check for ratio .gt. 0 and values above the limit
            if(ratio.le.0.0.and.ratiolimit.gt.0.0) ierr=1
            ic=nfirst-1
            do while(ic.lt.nlast.and.ierr.eq.0)
                ic=ic+1
                if(smocount(ic).lt.smocountlimit) ierr=2
            end do
        end if

        nchfir=nfirst
        nchlas=nlast
            

        return
    
    end subroutine smoov_limit

    subroutine smootophat(nfirst,nlast,nsmoos,nsmoomaxs,detcount,smocount,errcount,ratio,ratiolimit,nchanout)

        integer nfirst,nlast,nsmoo,nsmoos,nsmoomax,nsmoomaxs,nchanout
        real, dimension(:)        :: detcount,errcount,smocount
        real ratio,ratiolimit

        integer ic,ic1,jref,nup,ndown,nbroad,npoints,nsum,nincr,nsmoonew
        integer itimes,ntimes,nsmooold,nsmoonewmin,nsmoonewmax
        real sumsqrdif,sum,add,dif,ratioratio,ratioold,grad

        logical ratiotest

!c Set the initial degree of smoothing

        nincr=2
! Number of points to be smoothed
        npoints=nlast-nfirst+1
        !Set the maximum number of smooths allowed for this range. 
        !The limit of smoothing is determined so that npoints >= 2*nsmoo+1. 
        !Hence spectra with less that 3 channels cannot be smoothed
        !User can also control this limit with nsmoomaxs
        nsmoomax=min(nsmoomaxs,int((npoints-1)/20))
        nsmoo=nsmoos
        if(nsmoo.gt.nsmoomax) nsmoo=nsmoomax
        if(npoints.gt.0) then
!c Initialise the array to be smoothed
            do ic=nfirst,nlast
                smocount(ic)=detcount(ic)
            end do
            if(nsmoo.eq.0) then !Return a ratio of 1.0 if there are real data, but no smoothings are requested or possible
                ratio=1.0
                return
            end if
        else
            ratio=0.0 !Return a ratio of 0.0 if there are no points to smooth.
            return
        end if
!c Start with the minimum of smoothing
        ratio=0.0
        ratiotest=ratio.lt.ratiolimit.and.npoints.gt.0
        itimes=0
        ntimes=100
        do while (ratiotest.and.itimes.lt.ntimes)
!Ensure number of smooothings is compatible with range of values
            if(nsmoo.gt.nsmoomax) nsmoo=nsmoomax
            nup=nsmoo
            ndown=-nsmoo
            nbroad=nup-ndown+1
            sumsqrdif=0.0
            ic=nfirst-1
            do while (ic.lt.nlast)
                ic=ic+1
                sum=0.0
                nsum=0
                do ic1=ndown,nup
                    jref=ic+ic1
!c Try to compensate for end effects in the smoothing
                    if(jref.lt.nfirst) jref=2*nfirst-jref
                    if(jref.gt.nlast) jref=2*nlast-jref
!                    if(jref.lt.nfirst.or.jref.gt.nlast) write(6,*) ic,ic1,nfirst,nlast,nsmoo,jref
                    add=detcount(jref)
                    sum=sum+add
                    nsum=nsum+1
                end do
                smocount(ic)=sum/real(nsum)
                if(errcount(ic).gt.0.0) then
                    dif=detcount(ic)-smocount(ic)
                    sumsqrdif=sumsqrdif+dif*dif/errcount(ic)
                endif
            end do
            itimes=itimes+1
            ratio=sumsqrdif/real(npoints)
!            write(nchanout,*) 'smootophat> ',itimes,ratio,ratiolimit,nfirst,nlast,nsmoo,nsmoomax
!c If the ratio is out of range of the the limit, then alter the number of smoothings accordingly
!c Otherwise if ratio is .le. 0 then abort the smooth.
            if(ratio.gt.0.0) then
                ratiotest=ratio.lt.ratiolimit
                if(ratiotest) then
                    !Exit the loop if we already have used all the available smoothings
                    ratiotest=nsmoo.lt.nsmoomax
                    if(ratiotest) then
                        nsmoonew=nsmoo*1.1
                        if(nsmoonew.eq.nsmoo) then
                            nsmoo=nsmoonew+1
                        else
                            nsmoo=nsmoonew
                        endif
                    end if
                else
!c If over the limit, then reduce the number of smoothings.
!c This will not affect the current result, but will affect
!c subsequent results.
                    if(ratio.gt.ratiolimit) nsmoo=nsmoo*0.5
                endif
            else
                ratiotest=.false.
            endif
        end do
        nsmoos=nsmoo !Save the final value for future smooths.

        return
    end subroutine smootophat
    
    subroutine fit_line(nfirst,nlast,xin,yin,ein,wght,mergepwr,grad,cons)
        
        real, dimension(:)          :: xin,yin,ein,wght
        real                        :: sum,sumx,sumxsq,sumd,sumdx,wt,x,xsq,denom,grad,cons
        integer                     :: nfirst,nlast,i,mergepwr
!c
!c Form the various sums
!c
        sum=0.0
        sumx=0.0
        sumxsq=0.0
        sumd=0.0
        sumdx=0.0
        do i=nfirst,nlast
            if(abs(yin(i)).gt.0.0.and.ein(i).gt.0.0) then
                wt=(wght(i)**mergepwr)/ein(i)
                x=xin(i)
                xsq=x*x
                sum=sum+wt
                sumx=sumx+wt*x
                sumxsq=sumxsq+wt*xsq
                sumd=sumd+wt*yin(i)
                sumdx=sumdx+wt*yin(i)*x
            endif
        end do
!c
!c calculate gradient and constant
!c
        denom=sumx*sumx-sumxsq*sum
        if(denom.ne.0.0) then
            grad=(sumd*sumx-sumdx*sum)/denom
        else
            grad=0.0
        endif
        if(sum.gt.0.0) then
            cons=(sumd-grad*sumx)/sum
        else
            cons=0.0
        endif
!c
!c calculate line to output
!c
        do i=nfirst,nlast
            x=xin(i)
            yin(i)=cons+grad*x
        end do
        return
      
    end subroutine fit_line

    subroutine tophat3d(nsmoo,nchan,detcount,errcount)
!c
!c does a square-wave smoothing with width proportional to total number of points in array
!c input in detcount, smoothed result in errcount
!c
	real, dimension(:)                  :: detcount,errcount
        integer                             :: nsmoo,nchan,nfirst,nlast,ic,n,m1,m,mmn,mpn,mmax,l,lref,mmin
        real                                :: rn,rn2,factv,rm,rm2,rmn1,factt,add,rl,rl2,fact
!c
!c determine the first and last non-zero bins of input array
!c
	call trim_zeros(1,nchan,nfirst,nlast,detcount)
!c
!c smooth the input array with a top hat function in 3D
!c
	n=nsmoo
	rn=n
	rn2=2.0*rn+1
	factv=1.0/(rn2*rn2*rn2)
	do m1=nfirst,nlast
            errcount(m1)=0.0
            m=m1-1
            rm=m
            rm2=2.0*rm   
            rmn1=rm2*rm2-rn2*rn2+1.0
            factt=factv/rm2
            mmn=m-n
            mpn=m+n
!
!c if mmn .lt. 1, then need to include the cases when the whole volume is 
!c included
!c
            if(mmn.lt.1) then
                mmax=iabs(mmn)
!
!c volume of smallest sphere
!c
                l=1
                add=0.0
		if(l.lt.nfirst) then 
                    add=detcount(nfirst)
		else
                    add=detcount(l)
		endif
		errcount(m1)=errcount(m1)+add*factv
	   	if(mmax.gt.0) then
                    do l=1,mmax
			lref=l+1
                        rl=l
                        rl2=12.0*rl*rl+1.0
			fact=2.0*rl2*factv
			if(lref.le.nfirst) then
                            add=detcount(nfirst)
			else if(lref.ge.nlast) then
                            add=detcount(nlast)
			else
                            add=detcount(lref)
			end if
                        errcount(m1)=errcount(m1)+add*fact
                    end do
	   	end if
	   	mmin=mmax+1
	   else
	   	mmin=mmn		
	   endif
	   mmax=mpn
	   do l=mmin,mmax
		lref=l+1
		rl=l
		rl2=rl*rl
		fact=rm2*(12.0*rl2+1.0)-3.0*rl*(4.0*rl2+rmn1)
		fact=fact*factt
		if(lref.lt.nfirst) then
                    add=detcount(nfirst)
		else if(lref.gt.nlast) then
                    add=detcount(nlast)
		else
                    add=detcount(lref)
		end if
		errcount(m1)=errcount(m1)+fact*add
	   end do
	end do
	return
    end subroutine tophat3d

    subroutine tophat3dmod(nchan,delta,rmin,cons,xin,yin,yout,expamp,expdecay,expstretch)
!c
!c does a square-wave smoothing with width proportional to total number of points in array
!c input in detcount, smoothed result in errcount
!c
	real, dimension(:)      :: xin,yin,yout
        real                    :: delta,rmin,cons,expamp,expdecay,expstretch,r,delr,delr3,pofr,expon
        integer                 :: nchan,imin,ic

!c Determine the first channel at or after rmin

	imin=0
	ic=0
	do while (ic.lt.nchan.and.imin.eq.0)
            ic=ic+1
            if(xin(ic).ge.rmin) imin=ic
	end do
	do ic=1,nchan
            r=xin(ic)
            yout(ic)=0.0
            if(r.lt.rmin) then
                yout(ic)=yin(ic)+cons
            endif
            delr=delta*r
            if(delr.gt.0.0) then
                delr3=3.0/(delr*delr*delr)
                pofr=delr3*(sin(delr)-delr*cos(delr))
                if(r.lt.rmin) then
!c		 yout(ic)=yout(ic)-cons*pofr
                else
                    yout(ic)=-yin(ic)*pofr/(1.0-pofr)
                    expon=expdecay*r**expstretch
                    if(expon.gt.0.and.expon.lt.30.0) then
			yout(ic)=yout(ic)+expamp*exp(-expon)
                    endif
                endif
            endif
	end do
	return
        
    end subroutine tophat3dmod
        
    subroutine stogtos(ntype,dens,lptin,xin,yin,lptout,xout,yout)
	
        real, dimension(:)              :: xin,yin,xout,yout
        real                            :: pi,dens,trfac,x1,x2,r1,r2,r,rfac,sum,q1,q2,qdel,qr,sinqr
        integer                         :: ntype,lptin,lptout,i,j
        
        
        pi=4.0*atan(1.0)
       
!c
!c ntype=1 means gtos, otherwise stog
!c
	if(ntype.eq.1) then
            trfac=4.0*pi*dens
	else
            trfac=0.5/pi/pi/dens
	endif
!c
!c note that input array is modified by the program
!c x values are assumed to be bin boundary values throughout
!c
	x1=xin(1)
	do i=1,lptin-1
            x2=xin(i+1)
            yin(i)=0.5*(x1+x2)*yin(i)
            x1=x2
	end do
	yin(lptin)=0.0
	r1=xout(1)
	do i=1,lptout-1
            r2=xout(i+1)
            r=0.5*(r1+r2)
            rfac=trfac/r
            sum=0.
            q1=xin(1)
            do j=1,lptin-1
		q2=xin(j+1)
		qdel=q2-q1
		qr=0.5*r*(q1+q2)
		sinqr=sin(qr)
		sum=sum+yin(j)*sinqr*qdel
                q1=q2
            end do
            yout(i)=rfac*sum
            r1=r2
	end do
	yout(lptout)=0.0
	return
    end subroutine stogtos

    subroutine stogtoslorch(ntype,dens,lptin,xin,yin,lptout,xout,yout,qwindow,broad)

!c Integrate Q sin Qr over range for each bin

	real                            :: broad   !Determines how the broadening changes with r
        real, dimension(:)              :: xin,yin,xout,yout
        real                            :: pi,dens,trfac,rfac,qdel,qr,sinqr
        integer                         :: ntype,lptin,lptout,i,j
        real(kind=8)                    :: q1,q2,q,r1,r2,r,sum
        real(kind=8)                    :: j1q1r,j1q2r
        real                            :: qwindow,bwindow,bwidth

	pi=4.0*atan(1.0)
!c
!c ntype=1 means gtos, otherwise stog
!c
	if(ntype.eq.1) then
            trfac=4.0*pi*dens
	else
            trfac=0.5/pi/pi/dens
	endif
	r1=xout(1)
	do i=1,lptout-1
            r2=xout(i+1)
            r=0.5*(r1+r2)
            rfac=trfac/(r*r*r)
            bwindow = qwindow*(1.0+r**broad)
            sum=0.
            q1=xin(1)
            j1q1r=j1qr(q1,r)
            do j=1,lptin-1
		q2=xin(j+1)
                j1q2r=j1qr(q2,r)
                bwidth=0.5*(q1+q2)*bwindow
!c Form integral of (Q sin Qr)/r from r1 to r2
		sum=sum+yin(j)*(j1q2r-j1q1r)*pofr(bwidth)
                q1=q2
                j1q1r=j1q2r
            end do
            yout(i)=rfac*sum
            r1=r2
 	end do
	yout(lptout)=0.0
	return
        
    end subroutine stogtoslorch

    function j1qr(q,r)

        real(kind=8) j1qr,q,r,qr

        qr=q*r
        j1qr = (sin(qr)-qr*cos(qr))

        return
        
    end function j1qr

    function pofr(delr)
        
        real pofr,delr,delr3

        if (delr.le.0.0) then
            pofr=1.0
        else
            delr3=3.0/(delr*delr*delr)
            pofr=delr3*(sin(delr)-delr*cos(delr))
        endif
        return

    end function pofr

    SUBROUTINE SMOOe(NQMIN,NQMAX,NSMOO,nsmoot,DIFF,smodif,wt,ratio,ratiolimit)

!c NSMOO is the maximum number of smoothings that will be performed - this is to stop
!c infinite loops in the event the specified chisq ratio limit is not reached.

	real DIFF(*),smodif(*),wt(*)
        real chisq,ratiolimit,ratio
        integer nqmin,nqmax,nsmoo,nsmoot,n1,n2,np,ismoo,i,j,nadd
        real sf,sm,sn,dif
!C
!C 3-POINT SMOOTHING keeping integrated area constant
!C
	n1=nqmin
	n2=nqmax
!c
!c put data into the smoothing array
!c
	do i=n1,n2
            smodif(i)=diff(i)
	end do
	np=n2-n1+1
	nsmoot=0
	if(n2.lt.n1.or.nsmoo.eq.0) return
	do ismoo=1,iabs(nsmoo)
            SF=smodif(n1)
            SM=smodif(N1)
            smodif(n2+1)=smodif(n2)
            DO I=N1,N2
                J=I+1
                SN=smodif(J)
                smodif(I)=0.25*(SF+2.*SM+SN)
                SF=SM
                SM=SN
            end do
            chisq=0.0
            nadd=0
            do i=n1,n2
                if(wt(i).gt.0.0) then
                    nadd=nadd+1
                    dif=diff(i)-smodif(i)
                    chisq=chisq+dif*dif/wt(i)
                endif
            end do
!c
!c quit smoothing if chisq ratio is greater than specified limit
!c
            nsmoot=nsmoot+1
            ratio=0
            if(nadd.gt.0) ratio=chisq/real(nadd)
            if(ratio.gt.ratiolimit) return
	end do
	return
    END subroutine smooe

    subroutine polyfit_nofix(ndata,x,ydata,nord,fit,nits,stepmc1,temp,mdata)
        
        use ran1_routines
        
!c Program to fit a polynomial of degree n to a set of data
!c It is assumed the data exist in the range 0 - 1.
        integer ndata,nord,ncoeff,nits,mdata,mord,nrej,noutput
        parameter (mord=50)
        real coeff(0:mord),coeffsave(0:mord)
        real x(mdata),ydata(mdata),fit(mdata),stepmc,stepmc1,temp
        real chisq,chisqsave,chisum,chisumold,frac,frac1,accur,change,diff1,temp1,test,test1,test2
        integer i,it,ichange
        
        stepmc=stepmc1
!c Fractional contribution to chisum
        frac=1.0/3.0
        frac1=1.0-frac
        if(nord.gt.mord) nord=mord
        ncoeff=nord+1
!c Initialise the random number generator
        call initran1(12)
        close(12)
!c Read in existing coefficients if present
!c      open(12,file='polyfitcoeff.txt',status='old',iostat=ierr)
!c      i=0
!c      do while(i.lt.mord.and.ierr.eq.0)
!c         read(12,*,iostat=ierr) i,coeff(i)
!c      end do
!c      close(12)
!c      nord=i
!c Must have at least as many data points as coefficients
        if(ndata.lt.ncoeff) then
            write(6,104)
104         format(/' polyfit_nofix> Insufficient data to fit this order polynomial to')
            ncoeff=ndata
            nord=ncoeff-1
        endif
!c set coefficients to zero
        do i=0,nord
            coeff(i)=0.0
        end do
!c Save the initial coefficients
        do i=0,nord
            coeffsave(i)=coeff(i)
        end do
        chisqsave=0
        chisum=0
        chisumold=100
        accur=0.001
        ichange=1
        nrej=0
        noutput=10000
        it=0
        do while(it.lt.nits.and.(abs(chisum-chisumold)/chisumold).gt.accur)
            it=it+1
!c Calculate the fit to the data and chisq
            chisq=0
            do i=1,ndata
                fit(i)=poly(nord,coeff,x(i))
                diff1=ydata(i)-fit(i)
                chisq=chisq+diff1*diff1
            end do
            if(chisqsave.eq.0.0) chisqsave=chisq
            temp1=temp*chisqsave
            test=(chisq-chisqsave)/temp1
            if(test.le.0) then
                test1=2.0
            else if(test.gt.30.0) then
                test1=0.0
            else
                test1=exp(-test)
            endif
            test2=ran1()
!c Restore the previous coefficients if the change in fit is too large
            if(test1.lt.test2) then
                nrej=nrej+1
                coeff(ichange)=coeffsave(ichange)
            else
!c Otherwise save the new coefficients
                coeffsave(ichange)=coeff(ichange)
                chisqsave=chisq
            endif
!c Choose a coefficient at random and change it. However do not change anything on the last iteration, otherwise it will not
!c be tested
            if(it.lt.nits) then
                test=ran1()*real(nord+1)-0.4999999999
                ichange=nint(test)
                if(ichange.lt.0.or.ichange.gt.nord) then
                    write(6,*) ' polyfit-nofix> Error: ',test,ichange
                    stop
                endif
                change=stepmc*(2.0*ran1()-1.0)
                coeff(ichange)=coeff(ichange)+change
            endif
!c write output every 10000 iterations
            if(noutput*(it/noutput).eq.it.or.it.eq.nits) then
                chisumold=chisum
                chisum=frac1*chisum+frac*chisq
!                write(6,103) it,chisq,chisum,chisumold,real(nrej)/noutput
!103   format(' polyfit> Iteration ',i9,' Chisq values:',3(1x,e12.5),1x,f8.5)
!c aim for 50% rejection rate
                if(nrej.gt.0) then
                    stepmc=stepmc*real(noutput)*0.5/real(nrej)
                else
                    stepmc=stepmc*2.0
                endif
                nrej=0
            endif
        end do
        write(6,*) 'polyfit_nofix> No. of iterations ',it
!c write out coefficients
        open(12,file='polyfitcoeff.text',status='unknown')
        do i=0,nord
            write(12,105) i,coeff(i)
105         format(1x,i5,1x,e12.5)
        end do
        close(12)
!c Get the final fit
        do i=1,ndata
            fit(i)=poly(nord,coeffsave,x(i))
        end do
!c Save the random number generator
        call saveran1(12)
        close(12)
        return
    end subroutine polyfit_nofix

    function poly(nord,coeff,x)
        real coeff(0:*),x,poly,sum
        integer i,nord
        i=nord
        sum=coeff(i)
        do while (i.gt.0)
!c      write(6,*) i,x,sum,coeff(i)
            i=i-1
            sum=coeff(i)+x*sum
        end do
!c      write(6,*) i,x,sum,coeff(i)
        poly=sum
        return
    end function poly
    
    subroutine trim_zeros(nchmin,nchmax,nfirst,nlast,detcount)

        !Determines the first (nfirst) and last (nlast) non-zero values in the array detcount, within the range nchmin,nchmax
        integer, intent(in)                         :: nchmin,nchmax
        integer, intent(out)                        :: nfirst,nlast
        integer                                     :: i,j
        real, dimension(:), intent(in)              :: detcount
        logical                                     :: found,foundfirst
        
        !First check that indices are not out of range
        if((nchmin.lt.lbound(detcount,1)).or.(nchmax.gt.ubound(detcount,1))) then
            write(6,*) 'trim_zeros> ERROR: range indices are out of range. ',nchmin,lbound(detcount),nchmax,ubound(detcount)
            stop
        end if
        nfirst=nchmin-1
        nlast=nfirst-1 !Ensure npoints = nlast-nfirst+1 = 0 at outset
        i=nfirst
        foundfirst=.false.
        do while (i.lt.nchmax)
            i=i+1
            found=detcount(i).ne.0.0
            if(found) then
                if(.not.foundfirst) then
                    nfirst=i
                    foundfirst=.true.
                end if
                nlast=i
            end if
        end do
        return
        
    end subroutine trim_zeros

END MODULE corrections_routines

