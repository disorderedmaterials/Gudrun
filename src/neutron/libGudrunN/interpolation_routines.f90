!     
! File:   interpolation_routines.f90
! Author: aks45
!
! Created on 05 November 2013, 21:41
!

MODULE interpolation_routines
    
    use reallocation_routines
    
    implicit none
    
    CONTAINS
    
    FUNCTION AINTER(ALIN,COR,AL,NPTS)
        real, dimension(:), intent(in) :: ALIN,COR
        real al,ainter
        integer npts
        real x1,x2,y1,y2,grad
        integer i
        logical found

        if(npts.lt.2) then
            ainter=cor(1)
            return
        end if
!C  Interpolate COR at the point X=AL, using linear interpolation.
        x1=alin(1)
        y1=cor(1)
        x2=alin(2)
        y2=cor(2)
        found=al.lt.x2
        i=2
        do while (.not.found.and.i.lt.npts)
            i=i+1
            x1=x2
            y1=y2
            x2=alin(I)
            y2=cor(i)
            found=al.lt.x2
        end do
        if(i.lt.2.or.al.lt.x1) then
            write(6,*) ' ainter> Strange values: ',i,npts,al,x1,x2,y1,y2
!            stop
        end if
        GRAD=(Y2-Y1)/(X2-X1)
        AINTER=Y1+GRAD*(AL-X1)
        RETURN
    END function ainter
    
    subroutine rebinq(nold,xold,yold,nnew,xnew,ynew,iopt,iopt1)
!c
!c Performs a rebin from nold channels with x-values xold, to nnew channels
!c with x-values xnew. Partly overlapping channels are left empty at the 
!c beginning and end of the new array.
!c
!c The supplied x-values are assumed to be BIN BOUNDARIES throughout 
!c Thus there nold values of yold, and nold+1 values of xold. Same with nnew
!c ynew and xnew.
!!c 
!c
!c In this version, where the new bin boundary occurs between two old bin 
!c boundaries, linear interpolation is made between the neighbouring values,
!c rather than simple partitioning between the two.
!c
	integer nold		!number of old values
	integer nnew		!number of new values
	integer iopt		!1 = histograms,2 = rate, 3 = rate std. dev.
	integer iopt1		!0=no printout, 1 = printout
	integer oldref,newref	!intermediate counters
	integer i,inewref,oldmin,oldmax !range of old values to be interpolated over
	real, dimension(:) :: xold!(mold)	!old x- boundary values 
	real, dimension(:) :: yold!(mold)	!old y-values (nold values)
	real, dimension(:) :: xnew!(mnew)	!new x- boundary values 
	real, dimension(:) :: ynew!(mnew)	!new y-values (nnew values)
	real lbold,ubold 	!temporary store for lower and upper boundaries of current bin
	real lold,uold,lnew,unew,llim,ulim	!intermediate limit values
	real binfrac		!accumulates fractional contribution to bin
	real binwidth
	real bincentre
	real oldbinwidth,newbinwidth
	real fracthist,fractrate,fractratestdev,grad,add,cons
        real xoldl,xoldu,xoldbin,xnewl,xnewu,xnewbin,xdiff,xdiffmin
        
        
	if(iopt1.eq.1) then
            open(25,file='junk.dat',status='unknown')
	endif
!c
!c zero all new elements
!c
	do i=1,nnew
            ynew(i)=0.0
	end do

!c If the number of x values is less than 3, then routine will not work. 
!c Therefore simply find the nearest new bin which is closest to this old value

	if(nold.lt.3) then

            xoldl=xold(1)
            xoldu=xold(2)
!c
!c centre of current old bin
!c
            xoldbin=0.5*(xoldu+xoldl)
            xnewl=xnew(1)
            inewref=1
            xdiffmin=0.0
            do i=2,nnew
                xnewu=xnew(i)
!c Centre of new bin
            	xnewbin=0.5*(xnewu+xnewl)
               	xdiff=abs(xnewbin-xoldbin)
		if(i.eq.2.or.xdiff.lt.xdiffmin) then
                    xdiffmin=xdiff
                    inewref=i-1
		endif
		xnewl=xnewu
            end do
!c Only increment the bin if the difference is less than the width of the bin
            if(xdiffmin.le.0.5*(xnew(inewref+1)-xnew(inewref))) then
		ynew(inewref)=yold(1)
            endif
            return

	endif

!c Search the old data for zeros. The first non-zero value determines the
!c first value to be interpolated. The last non-zero value determines the 
!c last value to be interpolated.

	oldmin=1
	oldmax=1
	do i=1,nold
            if(yold(i).ne.0.0.and.oldmin.eq.1) oldmin=i
            if(yold(i).ne.0.0) oldmax=i+1 ! Because these are boundary values we look for the upper boundary for
            !the bin that is non-zero
	end do
!c	write(6,*) oldmin,oldmax
	oldref=oldmin
!c
!c lower and upper boundaries of the current old bin
!c
	lbold=xold(oldmin)
	ubold=xold(oldmin+1)
!c
!c centre of current old bin
!c

	lold=0.5*(lbold+ubold)
!c
!c lower and upper limits on next old bin
!c
	lbold=ubold
	ubold=xold(oldmin+2)
!c
!c centre of next old bin
!c
	uold=0.5*(lbold+ubold)
!c
!c define the boundaries of the first new bin
!c
	newref=1
	lnew=xnew(1)
	unew=xnew(2)
!c
!c Ensure upper limit of current old bin is ABOVE lower limit of current
!c new bin
! AND
!c Ensure lower limit of current old bin is BELOW or EQUAL to lower limit of 
!c current new bin
!c
        do while ((uold.le.lnew.and.oldref.lt.oldmax-2).or.(lold.gt.unew.and.newref.lt.nnew-1))
            if(uold.le.lnew) then
                oldref=oldref+1
                lold=uold
                lbold=ubold
                ubold=xold(oldref+2)
                uold=0.5*(lbold+ubold)
            else if(lold.gt.unew) then
		newref=newref+1
		lnew=unew
		unew=xnew(newref+1)
            endif
        end do
        binfrac=0.0
        do while (oldref.le.oldmax-2)

            grad=(yold(oldref+1)-yold(oldref))/(uold-lold)
            cons=yold(oldref)
!c
!c form the limits to be integrated over
!c
            ulim=min(unew,uold)
            llim=max(lnew,lold)
!c
!c set up bin width and centre
!c
            binwidth=ulim-llim
            bincentre=0.5*(ulim+llim)-lold
            oldbinwidth=uold-lold
            newbinwidth=unew-lnew
!c
!c Form the binning factor for the current old bin and add to new bin
!c	
            fracthist=binwidth/oldbinwidth
            fractrate=binwidth/newbinwidth
            fractratestdev=oldbinwidth/newbinwidth
            binfrac=binfrac+fractrate
!c
!c integrate the straight line between lold and uold
!c between the lower and upper boundaries
!c
            add=(cons+(grad*bincentre))
            if(iopt.eq.1) then
		ynew(newref)=ynew(newref)+add*fracthist
            else if(iopt.eq.2) then
		ynew(newref)=ynew(newref)+add*fractrate
            else
		ynew(newref)=ynew(newref)+add*fractrate*fractratestdev
            endif
            if(iopt1.eq.1) then
                write(25,*) oldref,newref,lold,uold,lnew,unew,llim,ulim
                write(25,*) binwidth,bincentre,cons,grad
                write(25,*) add,fractrate,fractratestdev,binfrac,ynew(newref)
            endif
!c
!c if ulim.eq.unew, time to move on to the next new bin
!c
            if(ulim.ge.unew) then
!c
!c for iopt=2 or 3 divide by bin fraction unless the bin is less than 90% filled
!c
                if(binfrac.gt.0.5) then !Ignore partially filled bins (usually at beginning and end) - set to zero
                    if(iopt.lt.3) then
                        ynew(newref)=ynew(newref)/binfrac
                    else 
                        ynew(newref)=ynew(newref)/(binfrac*binfrac)
                    endif
                else
                    ynew(newref)=0.0
                endif
!c
!c exit when all new bins are exhausted
!c
                if(newref.gt.nnew-1) then
                    if(iopt1.eq.1) close(25)
                    return
                endif
                newref=newref+1
                binfrac=0.0
                lnew=unew
                unew=xnew(newref+1)
            else
!c
!c exit when all old bins are exhausted
!c
                if(oldref.ge.oldmax-2) then
!c
!c zero the last new bin, in case it is only part filled
!c
                    if(binfrac.gt.0.5) then !Ignore partially filled bins (usually at beginning and end) - set to zero
                        if(iopt.lt.3) then
                            ynew(newref)=ynew(newref)/binfrac
                        else 
                            ynew(newref)=ynew(newref)/(binfrac*binfrac)
                        endif
                    else
                        ynew(newref)=0.0
                    endif
                    if(iopt1.eq.1) close(25)
                    return
                endif
!c
!c else move on to the next old bin
!c
                oldref=oldref+1
                lold=uold
                lbold=ubold
                ubold=xold(oldref+2)
                uold=0.5*(lbold+ubold)
            endif
        end do
	
    end subroutine rebinq
    
    function integ(nwav,wavecorr,sammon,mcorrwav)
	!Integrate data in sammon. It is assumed the data are in histogram format
        real integ
	real wavecorr(mcorrwav)
	real sammon(mcorrwav)
        real sum
        integer nwav,mcorrwav,ic
	sum=0.0
	do ic=1,nwav-1
            sum=sum+sammon(ic)*(wavecorr(ic+1)-wavecorr(ic))
	end do
	integ=sum
	return
    end function integ
    
    subroutine inter(nbin,xbin,ybin,nval,xval,yval)
        
!c
!c interpolates/extrapolates set of input values yval on grid xval onto grid
!c xbin, with output in ybin
!c
!c linear interpolation/extrapolation used throughout
!c
        integer nbin,nval      !no. of values in input and output arrays
        real, dimension(:), allocatable, intent(in)                 :: xval,yval,xbin      !input x and y values
        real, dimension(:), intent(out)                             :: ybin      !interpolated results
!        real, dimension(:), allocatable                             :: tbin
        real sumfirst,sumlast,xval1,xval2,xtest,frac
        integer ival,nvalfirst,nvallast,ival1,ival2,ibin

!        write(6,*) 'inter> 1',nbin,size(xbin),size(ybin),nval,size(xval),size(yval)
!        call reallocate1d_r(tbin,nbin)
!c
!c step through output x values. It is assumed both sets of x values are in
!c order of increasing x
!c
! Make sure first bin is not zero
        nvalfirst=0
        nvallast=0
        ival=0
        do while (ival.lt.nval)
            ival=ival+1
            if(yval(ival).ne.0.0.and.nvalfirst.eq.0) nvalfirst=ival
            if(yval(ival).ne.0.0) nvallast=ival
        end do
        if(nvalfirst.eq.0.or.nvallast.eq.0) then
        !If the input data are all zero, simply output zeros
            do ibin=1,nbin
                ybin(ibin)=0.0
            end do
            return
        end if
!        write(6,*) 'inter> ',nvalfirst,nvallast,size(xbin),size(ybin),size(xval),size(yval)
        ival1=nvalfirst  
        ival2=ival1+1
        xval1=xval(ival1)
        xval2=xval(ival2)
        do ibin=1,nbin
            xtest=xbin(ibin)
            ! If the xvalue corresponding to this bin is greater than current xval2, then
            ! step to the next input value if possible
            do while(xtest.gt.xval2.and.ival2.lt.nvallast)
                ival1=ival2
                ival2=ival2+1
                xval1=xval2
                xval2=xval(ival2)
            end do
            frac=(xtest-xval1)/(xval2-xval1)
!            write(6,*) 'inter> 2 ',ibin,nbin,ival1,ival2,nvalfirst,nvallast
!            write(6,*) 'inter> 3 ',xtest,frac,xval(ival1),xval(ival2),yval(ival1),yval(ival2)
!            write(6,*) 'inter> 4 ',ibin,lbound(tbin),ubound(tbin)
!            tbin(ibin)=0.0
            ybin(ibin)=yval(ival1)+frac*(yval(ival2)-yval(ival1))
!            write(6,*) 'inter> 5 ',ibin,nbin,ival1,ival2,xtest,frac,yval(ival1),yval(ival2)
        end do
!        write(6,*) 'inter> 6',nbin,size(xbin),size(ybin),size(tbin)
!        do ibin=1,nbin
!            ybin(ibin)=tbin(ibin)
!        end do
!        write(6,*) 'inter> 7',nbin,size(xbin),size(ybin),nval,size(xval),size(yval)
        return
        
    end subroutine inter
    
    subroutine get_interp_values(ang,nang,azi,nazi,theta,phid,rat,nfirst,nlast,ratp,npfirst,nplast)

!c Gets the theta and phi values at which to interpolate a set of corrections

        real, dimension(:) :: ang,azi
        real theta,phid,angdif,phidif
        real rat,ratp
        integer nang,nazi
        integer ic,nfirst,nlast,npfirst,nplast
!c			
!c determine the nearest pair of correction scattering angles to theta (which
!c are assumed to be in order of increasing angle).
!c
        ic=1
        nfirst=1
        do while (ic.lt.nang.and.theta.gt.ang(ic))
            nfirst=ic
            ic=ic+1
        end do
        nlast=ic
	angdif=ang(nlast)-ang(nfirst)
        if(angdif.gt.0.0) then
            rat=(theta-ang(nfirst))/angdif
        else
            rat=1.0
	endif
!c
!c determine the nearest pair of phi values, which are assumed to be in order of 
!c increasing angle
!c
	ic=1
	npfirst=1
	do while (ic.lt.nazi.and.phid.gt.azi(ic))
            npfirst=ic
            ic=ic+1
	end do
	nplast=ic
	phidif=azi(nplast)-azi(npfirst)
	if(phidif.gt.0.0) then
            ratp=(phid-azi(npfirst))/phidif
	else
            ratp=1.0
	endif
        return
        
    end subroutine get_interp_values
    
    subroutine do_interp_corr(nchan,wavebin,corrbin,nwav,wav,corr,rat,nfirst,nlast,ratp,npfirst,nplast)

!c Interpolates corrections corr onto specified theta and phi

        real, dimension(:,:,:), intent(in)              :: corr!(mcorrwav,mcorrang,mcorrang)
        real, dimension(:), allocatable, intent(in)     :: wav,wavebin
        real, dimension(:), intent(out)    :: corrbin!(mcorrwav)
        real, dimension(:), allocatable                 :: corrtemp
        real cor1,cor2,cordif1,cordif2,rat,ratp
        integer nwav,nchan
        integer ic,nfirst,nlast,npfirst,nplast
!c
!c set up (by interpolation or extrapolation) the correction at this angle
!c
	call reallocate1d_r(corrtemp,nwav)
        do ic=1,nwav
!c
!c interpolate the correction to the current scattering angle
!c
            cordif1=corr(ic,nlast,npfirst)-corr(ic,nfirst,npfirst)
            cordif2=corr(ic,nlast,nplast)-corr(ic,nfirst,nplast)
            cor1=corr(ic,nfirst,npfirst)+rat*cordif1
            cor2=corr(ic,nfirst,nplast)+rat*cordif2
!c
!c interpolate the correction to the current azimuthal angle
!c
            corrtemp(ic)=cor1+ratp*(cor2-cor1)
        end do
!c
!c interpolate correction onto the same wavelength scale as data
!c
	call inter(nchan,wavebin,corrbin,nwav,wav,corrtemp)
        return
        
    end subroutine do_interp_corr

      
END MODULE interpolation_routines
      
      
