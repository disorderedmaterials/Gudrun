!     
! File:   merge_arrays.f90
! Author: aks45
!
! Created on 31 October 2013, 13:48
!

MODULE merge_routines
    
    
    implicit none
    
    integer ngroupt                  !defines no. of groups formed
    integer nq,mq                  !number of Q values
    integer nweighterr            !1 means to use error bars in merge, else 0
    integer nqfirst            !gives the first element of the final merged data that is to be output
    integer outputunitstype     !1 = Q [1/A], 2 = d-space [A], 3 = wavelength [A], 4 = energy [meV], 5 = TOF [musec]
    real, dimension (:,:), allocatable  ::       aggweights      !aggregated weights array
    real, dimension (:,:), allocatable  ::       aggsweights      !aggregated sample weights array
    real, dimension (:,:), allocatable  ::       aggeweights      !aggregated error weights array
    real, dimension (:,:), allocatable  ::       aggdweights      !aggregated deviation weights array
    real, dimension (:), allocatable    ::         qbound            !q boundary values
    real, dimension (:), allocatable    ::         qvalue            !q values
    real grbroad  !Broadening power: 0 = constant broadening, 0.5 = sqrt(r) broadening, 1 = r broadening of g(r)
    real finalgrstep  !Step size for final g(r)
    logical                             :: logbinning,logrbinning,dofr
    real                                :: qval
    
    CONTAINS
    
!***********************************************************************************
!*
!*	init_merge.FOR
!*
!*	A K Soper, May 2001
!*
!*	initiates the merge of data from good detectors
!*
!***********************************************************************************
    subroutine init_merge(qmin,aqstep,qmaxsave,aqpower,iflag)

    use reallocation_routines
    use inputfilestrings
    use run_par
    use calibration_routines
    use groups_routines
    use spec_van
    use spec_sam
!c
!c internal variables
!c
	integer i,ic,ic1,ic2,id,is,iq	!internal indices
	integer ig,j,jf,jl,jref	!internal indices
	integer iflag			!0 = only 1 group, 1 = ngroup of results
	integer nqprime		!used to test the number of q values with logarithmic binning
	real q1,q2,qconv,pi,pi4,tconv	!temporary values
	real qmin,aqstep,qstep,qmax,qmaxsave,qstep2,qpower,aqpower		!qmin, qstep, and qmax for final merge
	real rat,constant,delta,deltaprime	!used to test the number of q values with logarithmic binning
!c
!c half step size
!c
	if(aqstep.lt.0.0) then
            logbinning = .true.
            qstep = abs(aqstep)
            qpower=aqpower
        else
            logbinning = .false.
            qstep = aqstep
            qstep2 = 0.5*qstep
            qpower=aqpower
        endif
      
!c Get Qmax. We store this separately in case it is modified by this routine

	qmax=qmaxsave

!c Values for modified log binning

        delta=qstep
        deltaprime=qpower
!c
!c set up q array
!c
        if(logbinning) then

!c logarithmic binning: we will assume logarithmic binning from qmin to qmax, with the step size increasing
!c according to Q**qpower

            constant=deltaprime/delta
            q1=qmin
            qstep=delta*tanh(constant*q1)
      
	else

!c linear binning

            constant=qstep
            q1=0.0
	
	endif

        nqfirst=1
!c
!c if iflag .eq. 0, then all the data is put in one group (used for SQRAW)
!c if iflag .ne. 0, then data are separated into groups according to 
!c groups_def.dat
!c
	if(iflag.eq.0) then
            ngroupt=1
	else
            ngroupt=ngroup
	endif

!c Set bins and bin boundaries

        mq=100
        call reallocate1d_r(qbound,mq)
        call reallocate1d_r(qvalue,mq)
        qbound(1)=q1
	ic=1
        nqfirst=1
        do while (qbound(ic).lt.qmax)
            i=ic+1
            if(i.gt.mq) then
                mq=mq+100
                call reallocate1d_r(qbound,mq)
                call reallocate1d_r(qvalue,mq)
            end if
            if(logbinning) then
            	q2=q1+qstep
                qstep=delta*tanh(constant*q2)
            else
		q2=(2*ic-1)*qstep2
            endif
            qbound(i)=q2
            qvalue(ic)=0.5*(q1+q2)
            q1=q2
            if(q2.lt.qmin) nqfirst=i+1
!            write(6,*) i,qbound(i),qvalue(ic)
            ic=i
	end do
        nq=ic
        qmax=qbound(nq)
!c
!c set accumulators to zero
!c
	call reallocate2d_r(aggweights,mq,ngroupt)
        call reallocate2d_r(aggsweights,mq,ngroupt)
        call reallocate2d_r(aggeweights,mq,ngroupt)
        call reallocate2d_r(aggdweights,mq,ngroupt)
        do ig=1,ngroupt
            do iq=1,nq
                aggweights(iq,ig)=0.0
                aggsweights(iq,ig)=0.0
                aggeweights(iq,ig)=0.0
		aggdweights(iq,ig)=0.0
            end do
	end do
!c
!c set the checksum accumulators to zero
!c
	call reallocate1d_r(speccheck,nspec)
	call reallocate1d_r(checksum,nspec)
        call reallocate1d_r(errchecksum,nspec)
        do i=1,nspec
            checksum(i)=0.0
            errchecksum(i)=0.0
	end do
	return
    END subroutine init_merge
    
!***********************************************************************************
!*
!*      merge_det.FOR
!*
!*      A K Soper, November 1999
!*
!*      merges the data from good detectors
!*
!***********************************************************************************
    subroutine merge_det(wavemin,wavemax,iflag,nspecwrt)

        use reallocation_routines
        use inputfilestrings
        use run_par
        use calibration_routines
        use groups_routines
        use spec_van
        use spec_sam
        use sam_par
        use local_data
        use corrections_routines
        use interpolation_routines
        use write_routines
!c
!c internal variables
!c
        integer                                 :: i,ic,ic1,ic2,id,is,iq      !internal indices
        integer                                 :: ig,j,jf,jl,jref,isref      !internal indices
        integer                                 :: nfirst,nlast            !indices for range of channels to use
        integer                                 :: iflag                  !flags whether 1 group or several needed
        integer                                 :: nspecwrt,nwrt1,nwrt2  !spectrum to write if needed
        integer                                 :: nsmoo            !number of smoothings for weights array
        integer                                 :: nup,ndown,nbroad,igrpref
        integer                                 :: nterm                  !no. of non-zero terms in checksum
        real, dimension(:), allocatable         :: sambin!(mq)            !binned sample data for a detector
        real, dimension(:), allocatable         :: sambinnosub!(mq)            !binned sample data for a detector with no single atom subtraction
        real, dimension(:), allocatable         :: serrbin!(mq)            !binned sample errors for a detector
        real, dimension(:), allocatable         :: vanbin!(mq)            !binned vanadium data for a detector
        real                                    :: lentot                  !total flight path for detector
        real                                    :: rat,err,add                  !temporary real values
        real                                    :: wave1,wave2            !temporary wavelengths
        real                                    :: wavemin,wavemax            !minimum and maximum wavelengths to use
        real                                    :: q1,q2,qconv,pi,pi4,tconv      !temporary values
        real                                    :: qmin,qstep,qmax            !qmin, qstep, and qmax for final merge
        real                                    :: qup,qlow,qav            !test values to check bins will be full
        real                                    :: sumn,sumx,sumd,sumdx,sumx2 !used to calculate straight line fit
        real                                    :: cons,grad,wtgroup            !parameters of straight line fit
        real                                    :: lentotgrp          !group lambda to Q conversion pars.
        real, dimension(:), allocatable         :: qresf,qresl!(mres)!min and max q values for specified resonances
        real                                    :: tempd,tempe,diffd      !used to check individual spectrum against merge
        real                                    :: sumdata,sumerr !used to calculate average level for each detector
        character(len=256)                      :: fname
        logical                                 :: accept!Used to define whether a Q value is to be used
!c
!c define pi
!c
        pi=4.0*atan(1.0)
        pi4=4.0*pi
        tconv=pi/360.0
        nchan=nchans
        nchanb=nchanbs
        call reallocate1d_r(tcb,nchanb)
        do ic=1,nchanb
              tcb(ic)=tcbs(ic)
        end do
        fname=runs(1,1)
! Initial number of smoothings on standard deviations - this is to remove excessive fluctuations
! in this value
        nsmoo=10
        do is=1,nspecproc
!c
!c Gudrun spectrum number
!c
            isref=specproc(is)
!c
!c detector number of this spectrum
!c
            id=detno(isref)
!c
!c group number of this spectrum
!c
            igrpref=igrp(isref)
            if(iflag.eq.0) then
                  ig=1
            else
                  ig=igrpref
            endif
!c
!c total flight path and Q conversion parameters
!c
            lentot=lenin+lendet(id)
            rat=0.0039554/lentot
!c
!c flight path for group
!c
            lentotgrp=lenin+lengrp(ig)
!
!c Setup the conversion from wavelength to output units
!
            if (outputunitstype.eq.1) then
                qconv=pi4*sin(tconv*ttheta(id))
                qup=qconv/wavemin
                qlow=qconv/wavemax
            else if (outputunitstype.eq.2) then
                qconv=0.5/sin(tconv*ttheta(id))
                qup=qconv*wavemax
                qlow=qconv*wavemin
            else if (outputunitstype.eq.3) then
                qconv=1.0
                qup=qconv*wavemax
                qlow=qconv*wavemin
            else if (outputunitstype.eq.4) then
                qconv=81.787
                qup=qconv/(wavemin*wavemin)
                qlow=qconv/(wavemax*wavemax)
            else if (outputunitstype.eq.5) then
                qconv=lentotgrp*sin(tconv*tthgrp(ig))/(0.0039554*sin(tconv*ttheta(id)))
                qup=qconv*wavemax
                qlow=qconv*wavemin
            endif
!c
!c generate the minimum and maximum q values of the resonances if present
!c
!c                  write(6,*) nreson
            if(nreson.gt.0) then
                call reallocate1d_r(qresf,nreson)
                call reallocate1d_r(qresl,nreson)
                do i=1,nreson
                    if(outputunitstype.eq.1) then
                        qresf(i)=qconv/wavresl(i)
                        qresl(i)=qconv/wavresf(i)
                    else if (outputunitstype.eq.4) then
                        qresf(i)=qconv/(wavresl(i)*wavresl(i))
                        qresl(i)=qconv/(wavresf(i)*wavresf(i))
                    else
                        qresf(i)=qconv*wavresf(i)
                        qresl(i)=qconv*wavresl(i)
                    endif
                end do
            endif

!c set up wavelength boundaries and convert them to the appropriate units

            call reallocate1d_r(wavebound,nchanb)
            if(outputunitstype.eq.1) then    
                wave1=rat*(tcb(nchanb)+deltat(id))      
                wavebound(1)=qconv/wave1
                do ic=2,nchanb
                      ic1=nchanb-ic+1
                      wave2=rat*(tcb(ic1)+deltat(id))
                      wavebound(ic)=qconv/wave2
                      wave1=wave2
                end do
            else if(outputunitstype.eq.4) then
                wave1=rat*(tcb(nchanb)+deltat(id))      
                wavebound(1)=qconv/(wave1*wave1)
                do ic=2,nchanb
                    ic1=nchanb-ic+1
                    wave2=rat*(tcb(ic1)+deltat(id))
                    wavebound(ic)=qconv/(wave2*wave2)
                    wave1=wave2
                end do
            else
                wave1=rat*(tcb(1)+deltat(id))      
                wavebound(1)=qconv*wave1
                do ic=2,nchanb
                    ic1=ic
                    wave2=rat*(tcb(ic1)+deltat(id))
                    wavebound(ic)=qconv*wave2
                    wave1=wave2
                end do
            endif
!c
!c get normalised data for this detector, sorting in order of increasing Q,
!c and saving the range of non-zero values

!c Reverse the order of values if output units are Q or energy

            call reallocate1d_r(detcount,nchanb)
            call reallocate1d_r(dettemp,nchanb)
            call reallocate1d_r(errcount,nchanb)
            if (outputunitstype.eq.1.or.outputunitstype.eq.4) then
                do ic=1,nchan
                    ic1=nchanb-ic
                    detcount(ic)=normsamdet(ic1,is,1)
                    dettemp(ic)=normsamdetnosub(ic1,is)
                    errcount(ic)=errsamdet(ic1,is,1)
                end do
            else
                do ic=1,nchan
                    ic1=ic
                    detcount(ic)=normsamdet(ic1,is,1)
                    dettemp(ic)=normsamdetnosub(ic1,is)
                    errcount(ic)=errsamdet(ic1,is,1)
                end do
            endif
! Remove fluctuations in error array (to avoid excessively weighting some contributions)
            if(nweighterr.eq.2) call simplesmootophat(nchan,nsmoo,errcount)
!c
!c write results if needed
!c
            if(isref.eq.nspecwrt) then
                call change_ext(fname,'premerge')
                call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
                call change_ext(fname,'premergenosub')
                call w_diag_file(fname,nchanb,wavebound,dettemp,errcount)
            endif
!c
!c form checksum for this spectrum
!c
            nterm=0
            sumdata=0.0
            sumerr=0.0
            do ic=1,nchan
                if(errcount(ic).gt.0.0) then
                    nterm=nterm+1
                    sumdata=sumdata+detcount(ic)
                    sumerr=sumerr+errcount(ic)
                endif
            end do
            if(nterm.gt.0) then
                checksum(isref)=sumdata/nterm
                sumerr=sumerr/real(nterm)
                errchecksum(isref)=sumerr
            else
                checksum(isref)=0.0
                sumerr=0.0
                errchecksum(isref)=0.0
            endif

!c
!c rebin data and errors onto a common Q scale, bearing in mind that data
!c are now a ratio rather than histogram
!c
            call reallocate1d_r(sambin,nq)
            call reallocate1d_r(sambinnosub,nq)
            call reallocate1d_r(serrbin,nq)
            call reallocate1d_r(vanbin,nq)
            call rebinq(nchan,wavebound,detcount,nq-1,qbound,sambin,2,0)
            call rebinq(nchan,wavebound,dettemp,nq-1,qbound,sambinnosub,2,0)
            if(isref.eq.nspecwrt) then
                call rebinq(nchan,wavebound,errcount,nq-1,qbound,serrbin,3,1)
            else
                call rebinq(nchan,wavebound,errcount,nq-1,qbound,serrbin,3,0)
            end if
!c
!c determine the first and last non-zero bins
!c
            call trim_zeros(1,nq,nfirst,nlast,sambin)

! Decide how to assign merge weights

            if(nweighterr.eq.2) then
! Assign weights per channel
                do iq=1,nq
                    if(serrbin(iq).gt.0.0) then
                        vanbin(iq)=1.0/serrbin(iq)
                    else
                        vanbin(iq)=0.0
                    endif
                end do
            else if(nweighterr.eq.1) then
! Individual detector weighting
                do iq=1,nq
                    if(serrbin(iq).gt.0.0) then
                        vanbin(iq)=1.0/serrbin(iq)
                    else
                        vanbin(iq)=0.0
                    endif
                end do
            else 
! Uniform weighting
                do iq=1,nq
                    if(serrbin(iq).gt.0.0) then
                        vanbin(iq)=1.0
                    else
                        vanbin(iq)=0.0
                    endif
                end do
            endif
!c
!c write results if needed
!c
            if(isref.eq.nspecwrt) then
                call change_ext(fname,'merge')
                call w_diag_file(fname,nq,qbound,sambin,serrbin)
                call change_ext(fname,'mergenosub')
                call w_diag_file(fname,nq,qbound,sambinnosub,serrbin)
                call change_ext(fname,'mergewts')
                call w_diag_file(fname,nq,qbound,vanbin,serrbin)
            endif
!c
!c Add results to aggregates using defined weighting function (vanbin)
!c
            if(nlast.gt.0) then
                ic=0
                do iq=nfirst,nlast
!c
!c test that this Q value is within the required Q-range, and, if resonances
!c are present, that it is outside the region of a resonance.
!c
                    accept=serrbin(iq).gt.0.0
                    if(accept) accept=qbound(iq+1).gt.qlow.and.qbound(iq).lt.qup
                    if(accept) accept=qbound(iq+1).gt.qmingrp(ig).and.qbound(iq).lt.qmaxgrp(ig)
                    i=0
                    do while(i.lt.nreson.and.accept)
                        i=i+1
!c
!c should reject the bin if the lower limit lies anywhere within the resonance
!c or if the upper limit lies anywhere within the resonance
!c
                        accept=qbound(iq).le.qresf(i).or.qbound(iq).ge.qresl(i)
                        if(accept) accept=qbound(iq+1).le.qresf(i).or.qbound(iq+1).ge.qresl(i)
                    end do
                    if(accept) then
                        aggweights(iq,ig)=aggweights(iq,ig)+vanbin(iq)
                        aggsweights(iq,ig)=aggsweights(iq,ig)+sambinnosub(iq)*vanbin(iq)
                        aggeweights(iq,ig)=aggeweights(iq,ig)+serrbin(iq)*vanbin(iq)*vanbin(iq)
                        aggdweights(iq,ig)=aggdweights(iq,ig)+sambin(iq)*vanbin(iq)
                    endif
                end do
            endif
        end do
        return 
        
    END subroutine merge_det
      
    subroutine simplesmootophat(nchan,nsmoo,detcount)

        use reallocation_routines
        use corrections_routines
        
        integer                             :: nfirst,nlast,nchan,nsmoo
        real, dimension(:)                  :: detcount!(*)
        real, dimension(:), allocatable     :: smocount!(nchan)
        real                                :: sum,add
        integer                             :: ic,ic1,jref,nup,ndown,nbroad,nsum

!c Initialise the array to be smoothed and detect first and last non-zero points

        call reallocate1d_r(smocount,nchan)
        call trim_zeros(1,nchan,nfirst,nlast,detcount)
! nchan must be greater than 2*nsmoo+1, otherwise errors will occur
        do while (nchan.lt.(2*nsmoo+1))
            nsmoo=nsmoo-1
        end do
! Only smooth if requested and there are enough channels
        if(nsmoo.gt.0) then
            nup=nsmoo
            ndown=-nsmoo
            nbroad=nup-ndown+1
            ic=nfirst-1
            do while (ic.lt.nlast)
                ic=ic+1
                sum=0.0
                nsum=0
                do ic1=ndown,nup
                    jref=ic+ic1
! Try to compensate for end effects in the smoothing
                    if(jref.lt.nfirst) jref=2*nfirst-jref-1
                    if(jref.gt.nlast) jref=2*nlast-jref+1
                    add=detcount(jref)
                    sum=sum+add
                    nsum=nsum+1
                end do
                smocount(ic)=sum/real(nsum)
            end do
! Save final result in input array 
            do ic=nfirst,nlast
                detcount(ic)=smocount(ic)
            end do
        end if
        deallocate(smocount)
        return
        
    end subroutine simplesmootophat     

!***********************************************************************************
!*
!*      finish_merge_test.FOR
!*
!*      A K Soper, May 2001
!*
!*      completes the merging of the spectra
!*
!*       Modified 20th December 2001 to calculate and subtract a background line 
!*      (constant and/or log Q) before the final merge of all groups into a 
!*      single spectrum.
!*
!***********************************************************************************
    subroutine finish_merge(acceptance,mergepwr,nbacksub,iflag,qshell,rshell,qwindow,rmax)

        use reallocation_routines
        use inputfilestrings
        use run_par
        use calibration_routines
        use groups_routines
        use spec_van
        use spec_sam
        use van_par
        use sam_par
        use corrections_routines
        use write_routines
!c
!c internal variables
!c
        integer i,ic,ic1,ic2,id,is,iq      !internal indices
        integer ig,j,jf,jl,jref      !internal indices
        integer nperrq            !period number of this run
        integer iflag                  !0 = only 1 group, 1 = ngroup of results
        integer xcode,ycode            !used by Genie load command
        integer nskip                  !how many lines to skip
        integer ngroupw            !no. of groups to write to .dcs file
        integer nfirst,nlast,nmin            !first and last bins for each group
        integer ngroupm            !no. of groups in final merge
        integer mergepwr            !power used to set Q weighting for merge
        integer nbacksub            !=1 if background subtraction is needed
        integer nonzero                  !no. of groups with non zero dcs level
        integer nbroad                  !no. of channels broadening for top hat function
        integer nqfirstsave
        integer nsumhighq           !Used to count number of terms for high Q sum
        integer nqwindow            !Broadening amount for final g(r)
        integer nrfinal             !Counts number of steps in final g(r)
        integer nqlast

        real rat,rats,rate,ratd,err      !temporary real values
        real wave1,wave2            !temporary wavelengths
        real wavemin,wavemax            !minimum and maximum wavelengths to use
        real q1,q2,qconv,pi,pi4,tconv      !temporary values
        real qmin,qstep,qmax            !qmin, qstep, and qmax for final merge
        real qup,qlow                  !test values to check bins will be full
        real lenav,tthav,phiav,wtgroup
        real, dimension(:), allocatable         :: lentemp,tthtemp,phitemp!(mgroup)
        real, dimension(:), allocatable         :: tempdat,backdat,expondat!(mq)     !temporary store of data
        real, dimension(:), allocatable         :: tempsq,tempgr!(mq)      !Stores merged S(Q) and g(r) for final Fourier transform
        real, dimension(:), allocatable         :: qlog!(mq)                  !loq(q) values
        real, dimension(:), allocatable         :: rtemp!(mq)                  !temporary r values for background subtraction
        real, dimension(:), allocatable         :: grad,cons,grplevel!(mgroup) !gradient and constant for each group
        real gradmrg,consmrg          !gradient and constant of merged dcs data.
        real sumcons,sumlogq,add      !relative weighting for cons. and logQ
        real dcslevel,testup,testdn
        real acceptance            !used to select groups to merge
        real bav,bavsq,bsqav,ratiobsqavoverbavsq      !average scattering length and square
        real gmpercc            !density in grams per cm**3
        real avlevel            !average dcs level of all groups
        real sambin,sambinnosub,serrbin,vanbin
        real qshell,rshell,rstep,radius      !minimum radius for background subtraction
        real qwindow,qfac
        real sumhighq, normalisationfactor,csratio,q13,q23,sumq3,sumiq3,term
        real rratio,rfac,rmin,rmax,r1,r2,rlogstep,arg     !used for calculating logarithmic r bins
        real test,threshold                  !used to decide whether a q-bin will be accepted in final merge of groups.
        real mergelevel

        character(len=256)                      :: fname,fnamediag            !name of file to write data to
        character(len=3)                        :: mrg                  !extension on merged file
        character(len=4)                        :: mmrg                  !extension on merged file
        character(len=256)                      :: baddetfname     !bad detector file full filename
        character(len=256)                      :: groupsfilefname !groups file full filename
        character(len=1)                        :: extensionchar     !extra character to extension to signify type of units

        pi=4.0*atan(1.0)

        if (outputunitstype.eq.2) extensionchar = 'd'
        if (outputunitstype.eq.3) extensionchar = 'w'
        if (outputunitstype.eq.4) extensionchar = 'e'
        if (outputunitstype.eq.5) extensionchar = 't'
!c
!c if nbacksub .gt. 0 choose the bin number closest to qshell
!c
        if(qshell.gt.0.0.and.outputunitstype.eq.1.and..not.logbinning) then
            nbroad=0
            i=1
            do while (i.le.nq.and.qbound(i).lt.qshell)
                i=i+1
            end do
            nbroad=i-1
            qshell=qbound(nbroad+2)
            write(6,*) nbroad,qshell
        else
            nbroad=0
        endif
!c
!c run number for sample - used for output
!c
        fname=runs(1,1)
        nperrq=npers(1)
        write(6,*) 'Run filename and period ',fname(1:len_trim(fname)),nperrq
        call reallocate1d_r(grad,ngroup)
        call reallocate1d_r(cons,ngroup)
!c
!c form final averages
!c
!c set up codes for use in GENIE 2 - these signify x-scale is Q, and yscale is
!c a ratio
!c
        xcode=6
        ycode=-1
!c
!c Now read in the full filenames (including directory path for the
!c Bad detector files and groups files
!c
        inquire(file="spec.bad",name=baddetfname)
        inquire(file="groups_def.dat",name=groupsfilefname)
!c
!c Determine the relative weighting of light hydrogen and non-light hydrogen
!c in the sample
!c
        sumcons=0.0
        sumlogq=0.0
!c
!c set the r step value to calculate Fourier Transforms
!c
        rstep=pi/qbound(nq)
!c
!c set up the log(Q) values
!c
        call reallocate1d_r(qlog,mq)
        call reallocate1d_r(rtemp,mq)
        do iq=1,nq
            if(qvalue(iq).gt.0.0) then
                qlog(iq)=alog(qvalue(iq))
            else
                qlog(iq)=0.0
            endif
!
!c Setup a temporary radius scale
!
            rtemp(iq) = rstep*(iq-1)

        end do
!c
!c Expected assymptotic value of DCS
!c
        dcslevel=sscatav(1)/4.0/pi
!c
!c set up limits for the acceptability of each group (currently +-20% of 
!c DCSLEVEL
!c
        testup=dcslevel*1000
        testdn=dcslevel*0.0
!c
!c density of sample in gm/cm**3
!c
        gmpercc=srho(1)*satwtav(1)/0.602217
!c
!c average scattering length of sample and its square
!c
        bav=sscatlenav(1)
        bavsq=bav*bav
        bsqav=sscatlensqav(1)
        ratiobsqavoverbavsq=bsqav/bavsq
!
!c Set up normalisation factor for final merged data
!
        normalisationfactor=1.0
        if(normalisationtype.eq.1.and.bavsq.ne.0.0) normalisationfactor=bavsq
        if(normalisationtype.eq.2.and.bsqav.ne.0.0) normalisationfactor=bsqav
!c Convert to per cm if required
        if(dofr.and.normalisationfactor.eq.1.0) normalisationfactor=normalisationfactor/srho(1)
!c
!c create a file to put the fitted line coefficients in for each group. Also
!c show the expected DCS level for this sample and relative proportion of light
!c hydrogen and other atoms
!c
        call change_ext(fname,'gud')
        if(outputunitstype.eq.1) then
            if(nperrq.gt.1) then
                write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
            end if
        else
            if(nperrq.gt.1) then
                write(fname,'(a,a1,i2.2)') fname(1:len_trim(fname)),extensionchar,nperrq
            else
                write(fname,'(a,a1)') fname(1:len_trim(fname)),extensionchar
            end if
        endif
        open(11,file=fname,status='unknown')
        write(11,98) fname(1:len_trim(fname))
        write(11,98) info(1:len_trim(info))
        write(11,98) uname(1:len_trim(uname))
        write(11,98) sttime(1:len_trim(sttime))
98      format(1x,a/)
        write(11,110) srho(1),gmpercc,bav,bavsq,bsqav,ratiobsqavoverbavsq,dcslevel
110     format(1x,'Number density of this sample (atoms/A**3) = ',e13.6 &
        /1x,'Corresponding density in g/cm**3 = ',f10.5 &
        /1x,'Average scattering length of the sample (10**-12cm) = ',f10.5 &
        /1x,'Average scattering length squared (barns) = ',e13.6 &
        /1x,'Average square of the scattering length (barns) = ',e13.6 &
        /1x,'Ratio of (coherent) single to interference = ',e13.6 &
        //1x,'Expected level of DCS [b/sr/atom] = ',f10.5 &
        //1x,'Group number,  first Q,   last Q,   level [b/sr/atom]',',   gradient in Q (%)'/)
!c
!c average levels and number of non-zero groups
!c
        avlevel=0.0
        nonzero=0.0
        sumhighq=0.0
        nsumhighq=0
        call reallocate1d_r(lentemp,ngroupt)
        call reallocate1d_r(tthtemp,ngroupt)
        call reallocate1d_r(phitemp,ngroupt)
        call reallocate1d_r(grplevel,ngroupt)
        call reallocate1d_r(tempdat,mq)
        call reallocate1d_r(backdat,mq)
        call reallocate1d_r(expondat,mq)
        call reallocate1d_r(tempsq,mq)
        call reallocate1d_r(tempgr,mq)
        do ig=1,ngroupt
            lentemp(ig)=lengrp(ig)
            tthtemp(ig)=tthgrp(ig)
            phitemp(ig)=phigrp(ig)
            do iq=1,nq
                rat=aggweights(iq,ig)
                if(rat.gt.0.0) then
                    rats=aggsweights(iq,ig)/rat
                    rate=aggeweights(iq,ig)/rat/rat
                    ratd=aggdweights(iq,ig)/rat
                else
                    rat=0.0
                    rats=0.0
                    rate=0.0
                    ratd=0.0
                endif
                err=sqrt(rate)
                aggweights(iq,ig)=rat
                aggsweights(iq,ig)=rats
                aggeweights(iq,ig)=err
                aggdweights(iq,ig)=ratd
            end do
!
!c Fit a line to these data. This line is weighted by Q**mergepwr to emphasize the large Q region
!c
!c find the first and last non-zero bin numbers for this group
!c
            do iq=1,nq
                tempdat(iq)=aggsweights(iq,ig)
            end do
            call trim_zeros(1,nq,nfirst,nlast,tempdat)
            if(nfirst.gt.0) then
                do iq=nfirst,nlast
                    tempdat(iq)=aggsweights(iq,ig)
                    backdat(iq)=tempdat(iq)
                    expondat(iq)=aggeweights(iq,ig)
                end do
                call fit_line(nfirst,nlast,qbound,backdat,expondat,qbound,mergepwr,grad(ig),cons(ig))
!If gradient is > 0, report the high Q limit
                if(grad(ig).gt.0.0) then
                    grplevel(ig)=cons(ig)+qbound(nlast)*grad(ig)
                else
                    grplevel(ig)=cons(ig)
                end if
                if(abs(grplevel(ig)).gt.0.0) then
                    grad(ig)=100.0*grad(ig)/abs(grplevel(ig))
                else
                    grad(ig)=0.0
                endif
                if(cons(ig).gt.0) then
                    nsumhighq=nsumhighq+1
                    sumhighq=sumhighq+cons(ig)
                end if
                write(11,599) ig,qbound(nfirst),qbound(nlast),grplevel(ig),grad(ig)
                write(6,599) ig,qbound(nfirst),qbound(nlast),grplevel(ig),grad(ig)
            else
                do iq=1,nq
                    aggweights(iq,ig)=0.0
                    aggsweights(iq,ig)=0.0
                    aggeweights(iq,ig)=0.0
                    aggdweights(iq,ig)=0.0
                enddo
                cons(ig)=0.0
                grad(ig)=0.0
                write(11,599) ig,0.0,0.0,cons(ig),grad(ig)
                write(6,599) ig,0.0,0.0,cons(ig),grad(ig)
            endif
599         format(1x,i4,9x,f9.4,1x,f9.4,5x,f10.5,12x,f8.4)
        end do
!c
!c write the DCS for individual groups
!c
        write(6,*) 'Formed mean DCS for groups'

!c Temporarily set nqfirst to 1 so that all the q values are output for the individual groups

        nqfirstsave=nqfirst
        nqfirst=1

!c Write out the unsubtracted groups - w_int

        mrg='dcs'
        call change_ext(fname,mrg)
        if(outputunitstype.eq.1) then
            if(nperrq.gt.1) then
                write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
            end if
        else
            if(nperrq.gt.1) then
                write(fname,'(a,a1,i2.2)') fname,extensionchar,nperrq
            else
                write(fname,'(a,a1)') fname(1:len_trim(fname)),extensionchar
            end if
        endif
!        write(6,*) fname(1:len_trim(fname))
        call w_int(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)

!c Calculate exponential background to subtract if required

        do iq=1,nq
            expondat(iq)=0.0
        end do
        if(nexpon.gt.0.and.outputunitstype.eq.1) then

!c The problem is that we want to subtract the exponential BEFORE performing the top hat smoothing
!c to avoid the unphysical behaviour at low Q being incorporated into the smoothed version of the data.
!c However we cannot do this to individual groups since they all have different Q dependencies. So instead
!c we generate the exponential function and put IT through the smoothing process. The difference between it
!c and the original exponential function will then be subtracted from the difference stored in aggsweights.

            do iq=1,nq
                do i=1,nexpon
                    if(expondecay(i).gt.0.0) then
                        qval=qvalue(iq)/expondecay(i)
                        if(qval.ge.0.0.and.qval.lt.40.0) then
                            expondat(iq)=expondat(iq)+exponamp(i)*exp(-qval)
                        endif
                    endif
                end do
            end do
        endif
        !Mean high Q limit
        if(nsumhighq.gt.0) then
            sumhighq=sumhighq/nsumhighq
        else
            sumhighq=0
        end if
!c
!c now merge the groups (except for TOF output units where this merge is meaningless)
!c
        do ig=1,ngroupt
            do iq=1,nq-1
                tempdat(iq)=aggsweights(iq,ig)
            end do
            call trim_zeros(1,nq-1,nfirst,nlast,tempdat)
!Only process between the first and last non-zero points
            if(nfirst.gt.0) then
                do iq=1,nq
                    if(iq.ge.nfirst.and.iq.le.nlast) then
                        tempdat(iq)=aggdweights(iq,ig)-expondat(iq)
                    else
                        tempdat(iq)=0.0
                    endif
                    backdat(iq)=tempdat(iq)
                end do
!c
!c Calculate the top hat background function and subtract it if needed.  
!c Then reincorporate the result into aggsweights
!c
                if(qshell.ne.0.0) then
                    if(nbroad.gt.0) then
                        call change_ext(fname,'sub')
                        call w_diag_file(fname,nq,qbound,tempdat,expondat)
                        call tophat3d(nbroad,nq,tempdat,backdat)
                        call change_ext(fname,'smo')
                        call w_diag_file(fname,nq,qbound,backdat,expondat)
                        do iq=nfirst,nlast
                            tempdat(iq)=tempdat(iq)-backdat(iq)
                        end do
                    else
!Simply use the last few non-zero points of this group to generate a constant background
                        sumhighq=0
                        nsumhighq=0
                        iq=(nfirst+nlast)/2-1
                        do while (iq.lt.nlast)
                            iq=iq+1
                            nsumhighq=nsumhighq+1
                            sumhighq=sumhighq+tempdat(iq)
                        end do
                        if(nsumhighq.gt.0) sumhighq=sumhighq/nsumhighq
                        do iq=nfirst,nlast
                            tempdat(iq)=tempdat(iq)-sumhighq
                        end do
                    end if                
                end if
                do iq=1,nq
                    rat=aggweights(iq,ig)
                    if(rat.gt.0.0) then
                        aggdweights(iq,ig)=tempdat(iq)
                    else
                        aggdweights(iq,ig)=0.0
                    endif
                enddo
            endif
        enddo

!c Write out the self scattering subtracted groups - w_intd

        mrg='sub'
        call change_ext(fname,mrg)
        if(outputunitstype.eq.1) then
            if(nperrq.gt.1) then
                write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
            end if
        else
            if(nperrq.gt.1) then
                write(fname,'(a,a1,i2.2)') fname(1:len_trim(fname)),extensionchar,nperrq
            else
                write(fname,'(a,a1)') fname(1:len_trim(fname)),extensionchar
            end if
        endif
!        write(6,*) fname(1:len_trim(fname))
        call w_intd(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
!c
!c form the merge of non-subtracted and background subtracted groups
!c
        lenav=0.0
        tthav=0.0
        phiav=0.0
        ngroupm=0

!c Set a threshold at which to include a particular group and Q value. This is determined
!c by the ratio of error bar to (non-background subtracted) value. Any non-background subtracted
!c values <= 0.0 will be ignored.
        threshold=acceptance
        do ig=1,ngroupt
!c
!c first check that the DCS level for this group is acceptable (currently set
!c to +- acceptance*DCSLEVEL)
!c
            avlevel=avlevel+grplevel(ig)
            ngroupm=ngroupm+1
            lenav=lenav+lengrp(ig)
            tthav=tthav+tthgrp(ig)
            phiav=phiav+phigrp(ig)
!
!c Determine the range of Q values for this group. This determined to be from the 
!c largest value below the minimum for which the threshold is not satisfied to the smallest value 
!c above the minimum for which the threshols is also not satisfied.
!c Only do this if the threshold is less than 1.0
!
            test=threshold
            nmin=0
!c First find the value which has the smallest error bar
            do iq=1,nq
                sambinnosub=aggsweights(iq,ig)
                serrbin=aggeweights(iq,ig)
                if(sambinnosub.gt.0.0) then
                    rat=serrbin/sambinnosub
                    if(iq.eq.1) test=rat
                    if(rat.lt.threshold.or.threshold.ge.1.0) then

                        if(rat.lt.test) then
                            test=rat
                            nmin=iq
                        endif
                    endif
                endif
            end do
!
!c Now find the nearest values either side of this minimum below the threshold
!
            nfirst=0
            nlast=0
            if(nmin.gt.0.and.threshold.lt.1.0) then
                nfirst=1
                nlast=nq
                do iq=1,nq
                    sambinnosub=aggsweights(iq,ig)
                    serrbin=aggeweights(iq,ig)
                    if(sambinnosub.gt.0.0) then

                        rat=serrbin/sambinnosub
                        if(rat.gt.threshold) then

                            if(iq.lt.nmin.and.iq.ge.nfirst) nfirst=iq+1
                            if(iq.gt.nmin.and.iq.le.nlast) nlast=iq-1

                        endif
                    else
                        if(iq.lt.nmin.and.iq.ge.nfirst) nfirst=iq+1
                        if(iq.gt.nmin.and.iq.le.nlast) nlast=iq-1
                    endif
                end do
            else
                nfirst=1
                nlast=nq
            endif
!            write(6,*) ig,nfirst,nlast,nmin,test
            if(nmin.gt.0.and.nlast.ge.nfirst) then
                wtgroup=1.0
                if(nweighterr.eq.1) then
                    if(lengrp(ig).gt.0.0) then
                        wtgroup=(lengrp(ig)**2)/ndetgr(ig) !Attempt to weight group according to solid angle and number of detectors in group
                    else
                        if(ndetgr(ig).gt.0) wtgroup=1.0/ndetgr(ig)
                    end if
                end if
                do iq=nfirst,nlast
                    rat=aggweights(iq,ig)*wtgroup
                    rats=aggsweights(iq,ig)
                    rate=aggeweights(iq,ig)
                    rate=rate*rate
                    ratd=aggdweights(iq,ig)
                    if(ngroupm.eq.1) then
                        aggweights(iq,1)=rat
                        aggsweights(iq,1)=rats*rat
                        aggeweights(iq,1)=rate*rat*rat
                        aggdweights(iq,1)=ratd*rat
                    else 
                        aggweights(iq,1)=aggweights(iq,1)+rat
                        aggsweights(iq,1)=aggsweights(iq,1)+rats*rat 
                        aggeweights(iq,1)=aggeweights(iq,1)+rate*rat*rat
                        aggdweights(iq,1)=aggdweights(iq,1)+ratd*rat
                    endif
                end do
            endif
        end do
!c
!c normalise sums to number of groups in merge
!c
        if(ngroupm.gt.0) then
            avlevel=avlevel/ngroupm
            lentemp(1)=lenav/ngroupm
            tthtemp(1)=tthav/ngroupm
            phitemp(1)=phiav/ngroupm
            do iq=1,nq
                rat=aggweights(iq,1)
                if(rat.gt.0.0) then
                    rats=aggsweights(iq,1)/rat
                    rate=aggeweights(iq,1)/rat/rat
                    ratd=aggdweights(iq,1)/rat
                else
                    rats=0.0
                    rate=0.0
                    ratd=0.0
                endif
                err=sqrt(rate)
                tempdat(iq)=rat
                backdat(iq)=rats
                expondat(iq)=err
                aggsweights(iq,1)=rats/normalisationfactor
                aggeweights(iq,1)=err/normalisationfactor
                aggdweights(iq,1)=ratd/normalisationfactor
                tempsq(iq)=aggdweights(iq,1)
            end do

!c Fit a line to the merged dcs data

            call trim_zeros(1,nq,nfirst,nlast,backdat)
            call fit_line(nfirst,nlast,qbound,backdat,expondat,qbound,mergepwr,gradmrg,consmrg)
!If gradient is > 0, report the high Q limit
            if(gradmrg.gt.0.0) then
                mergelevel=consmrg+qbound(nlast)*gradmrg
            else
                mergelevel=consmrg
            end if

!c Save the merged unsubtracted dcs

            ngroupt=1
!c
!c write out the merged unsubtracted data, over the requested q-range
!c
            if(outputunitstype.lt.5) then
                call change_ext(fname,'mdcs')
                if(outputunitstype.eq.1) then
                    if(nperrq.gt.1) then
                        write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
                    end if
                else
                    if(nperrq.gt.1) then
                        write(fname,'(a,a1,i2.2)') fname(1:len_trim(fname)),extensionchar,nperrq
                    else
                        write(fname,'(a,a1)') fname(1:len_trim(fname)),extensionchar
                    end if
                end if
!                write(6,*) fname(1:len_trim(fname))
                call w_int(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
                nqfirst=nqfirstsave
!c 
!c write out the merged subtracted data over the requested Q-range
!c
                call change_ext(fname,'msub')
                if(outputunitstype.eq.1) then
                    if(nperrq.gt.1) then
                        write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
                    end if
                else
                    if(nperrq.gt.1) then
                        write(fname,'(a,a1,i2.2)') fname(1:len_trim(fname)),extensionchar,nperrq
                    else
                        write(fname,'(a,a1)') fname(1:len_trim(fname)),extensionchar
                    end if
                end if
!                write(6,*) fname(1:len_trim(fname))
                call w_intd(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
            endif
!
!c if nbroad.gt.0 then correct for effect of smoothing function used to remove background
!
            if(outputunitstype.eq.1) then

            !Choose a constant for the background so that g(0) = -csratio
                csratio=bavsq/normalisationfactor
                sumq3=0.0
                sumiq3=0.0
                nqfirst=1
                nqlast=nq-1
                q1=qbound(1)
                q13=q1*q1*q1
                do iq=nqfirst,nqlast
                    q2=qbound(iq+1)
                    q23=q2*q2*q2
                    term=q23-q13
                    sumq3=sumq3+term
                    sumiq3=sumiq3+term*aggdweights(iq,1)
                    q1=q2
                    q13=q23
                end do
                sumhighq=(sumiq3+6.0*pi*pi*srho(1)*csratio)/sumq3
!
!c Form the average high Q limit and subtract it.
!
                if(nsumhighq.gt.0) then
                    sumhighq=sumhighq/nsumhighq
                    do iq=1,nq-1
                        if(aggdweights(iq,1).ne.0.0) then
                            aggdweights(iq,1)=aggdweights(iq,1)-sumhighq
                        endif
                    end do
                endif
            endif

!c Now finalise the interference function if requested

            if(outputunitstype.eq.1.and.(nbacksub.eq.1.or.qshell.ne.0.0)) then

!c First copy the subtracted data into the correct array

                do iq=1,nq
                    aggsweights(iq,1)=aggdweights(iq,1)
                end do
                do iq=1,nq
                    tempdat(iq)=aggsweights(iq,1)
                    backdat(iq)=0.0
                end do

!c
!c Fourier transform to r-space
!c
                call stogtos(2,srho(1),nq,qbound,tempdat,nq,rtemp,backdat)
!c
!c write diagnostic file
!c
                call change_ext(fname,'gr1')
                call w_diag_file(fname,nq,rtemp,backdat,tempdat)
!c
!c correct F.T. for effect of smoothing function, and apply a minimum radius
!c (recommended but only if specified)
!c
                call tophat3dmod(nq,qshell,rshell,csratio,rtemp,backdat,tempdat,0.0,1.0,1.0)
                call change_ext(fname,'gr2')
                call w_diag_file(fname,nq,rtemp,tempdat,backdat)
!c
!c Finally back transform to Q space and subtract the background function, applying normalisation factor
!c
                call stogtos(1,srho(1),nq,rtemp,tempdat,nq,qbound,backdat)
                call change_ext(fname,'bak')
                call w_diag_file(fname,nq,qbound,backdat,tempdat)
                do iq=1,nq
                    if(aggsweights(iq,1).ne.0.0) then
                        aggsweights(iq,1)=(aggsweights(iq,1)-backdat(iq))
                        tempsq(iq)=aggsweights(iq,1)
                    endif
                end do
!c
!c now Fourier transform difference to r space, and output results
!c
!c           First setup the new radius scale
                if(logrbinning) then
                    rmin = finalgrstep/log(2.0)
                    rlogstep=0.01
                    rratio = rlogstep*log(2.0)/finalgrstep
                    rtemp(1)=0.0
                    rtemp(2)=finalgrstep
                    i=2
                    do while (i.lt.mq.and.rtemp(i).lt.rmax)
                        arg=rratio*rtemp(i)
                        if(arg.lt.10.0) then
                            rstep=finalgrstep+rmin*log(cosh(arg))
                        else
                            rstep=rtemp(i)*rlogstep
                        endif
                        i=i+1
                        rtemp(i)=rtemp(i-1)+rstep
                    end do
                    nrfinal=i
                 else
                    rstep = finalgrstep
                    rtemp(i)=0.0
                    i=1
                    do while (i.lt.mq.and.rtemp(i).lt.rmax)
                        i=i+1
                        rtemp(i) = rstep*(i-1)
                    end do
                    nrfinal=i
                endif
                call stogtoslorch(2,srho(1),nq,qbound,tempsq,nrfinal,rtemp,tempgr,qwindow,grbroad)
!
!c Write g(r)
!
                call change_ext(fname,'mgor01')
                if(nperrq.gt.1) then
                    write(fname,'(a,i2.2)') fname,nperrq
                end if
                call w_diag_file(fname,nrfinal,rtemp,tempgr,tempsq)
!
!c Write out d(r)
!
                r1 = rtemp(1)
                rfac=4*pi*srho(1)
                do i = 2,nrfinal
                    r2 = rtemp(i)
                    tempgr(i-1)=0.5*(r1+r2)*rfac*tempgr(i-1)
                    r1=r2
                end do
                call change_ext(fname,'mdor01')
                if(nperrq.gt.1) then
                    write(fname,'(a,i2.2)') fname,nperrq
                end if

!c write this out as a diagnostic file

                call w_diag_file(fname,nrfinal,rtemp,tempgr,tempsq)


!c Convert the units of the output dcs file to cm**-1 if required. dofr signifies whether to do this transformation

                if(dofr) then
                    do i = 1,nq-1
                        aggsweights(iq,1)=aggsweights(iq,1)*srho(1)
                        aggeweights(iq,1)=aggeweights(iq,1)*srho(1)*srho(1)
                    end do
                endif
!c
!c create new file name for the merged interference data
!c
                call change_ext(fname,'mint')
                if(nperrq.gt.1) then
                    write(fname,'(a,i2.2)') fname,nperrq
                end if
!                write(6,*) fname(1:len_trim(fname))
!c
!c write out the interference data, over the requested q-range
!c
                nqfirst=nqfirstsave
                call w_int(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
            end if
            !If the gradient is negative report the level at the low Q end,
            !Otherwise report the level at the high Q end
            if(abs(mergelevel).gt.0.0) then
                gradmrg=100.0*gradmrg/abs(mergelevel)
            else
                gradmrg=0.0
            endif
! Calculate a suggested new tweak factor based on consmrg and expected level
            if(consmrg.gt.0.0) then
                tweak(1)=tweak(1)*dcslevel/mergelevel
            end if
            write(11,"(/1x,'No. of groups accepted for merge = ',i4)") ngroupm
            write(11,297) mergelevel,gradmrg
            write(6,'(/" No. of groups accepted for merge = ",i4)') ngroupm
            write(6,297) mergelevel,gradmrg
297         format(/1x,'Average level of merged dcs is ',f10.5,' b/sr/atom;' &
            ,//1x,'Gradient of merged dcs: ',f7.4,'% of average level.')
            avlevel=100.0*(mergelevel/dcslevel-1.0)
            if(avlevel.gt.10.0) then
                write(11,298) avlevel
                write(11,300)
                write(6,298) avlevel
                write(6,300)
            else if(avlevel.lt.-10.0) then
                write(11,299) (-avlevel)
                write(11,300)
                write(6,299) (-avlevel)
                write(6,300)
            else
                write(11,301) (100.0+avlevel)
                write(6,301) (100.0+avlevel)
            endif
298         format(/1x,'WARNING! This DCS level is ',f6.1,'% ABOVE expected level.')
299         format(/1x,'WARNING! This DCS level is ',f6.1,'% BELOW expected level.')
300         format(/1x,'Please check sample density, size or thickness, and composition.' &
            /1x,'If all is in order, then refer to your local contact for further advice')
301         format(/1x,'This DCS level is ',f6.1,'% of expected level')
            if(tweak(1).lt.100.0.and.tweak(1).ge.0.01) then
                write(11,302) tweak(1)
                write(6,302) tweak(1)
            else
                write(11,303) tweak(1)
                write(6,303) tweak(1)
            endif
302         format(/1x,'Suggested tweak factor: ',f10.5)
303         format(/1x,'Suggested tweak factor: ',e12.5)
        else
            write(11,"(/1x,'No. of groups accepted for merge = ',i4)") ngroupm
        end if
        close (11)
!c
!c write out the detector checksums
!c
        call change_ext(fname,'chksum')
        if(nperrq.gt.1) then
            write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
        endif
        do i=1,nspec
            speccheck(i)=i
        end do
        call w_diag_file(fname,nspec,speccheck,checksum,errchecksum)
        return      
    
    END subroutine finish_merge

    subroutine w_intd(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
        
        use run_par
        use calibration_routines
        
        integer, parameter :: nwrmax=100
        integer xcode,ycode            !used by Genie load command
        integer nskip      !how many lines to skip
        integer lrec,lrecf,ilen,iwr,iwr1,iwr2,ir,iq,it,nt,nblock     !record lengths
        real lentemp(*),tthtemp(*),phitemp(*),xo
        character*256 fname      !name of file to write data to
        character*256 baddetfname     !bad detector file full filename
        character*256 groupsfilefname !groups file full filename
        character*32 formatout
        character*2 blk
!c
! After the initial 3 lines nskip is the number of lines of info that need to be skipped
! to get to the actual information
!c
        nskip=10
!c
! set up number of output files to write to
!c
      nt=(ngroupt+(nwrmax-1))/nwrmax
      ilen=index(fname,' ')-1
      iwr1=1
      do it=1,nt
!c
! number of elements to write for each x value
!c
            iwr=ngroupt-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
            iwr2=iwr1+iwr-1
!c
! setup output format specification
!c
            nblock=2*iwr+1
            write(formatout,332) '(a1,',nblock,'(1x,e14.7))'
332      format(a4,i3.3,a11)
            write(6,*) formatout
            lrec=nblock*15+1
            if(lrec.lt.82) lrec=82
            write(blk,331) it
331      format(i2.2)
            fname=fname(1:ilen)//blk(1:2)
            write(6,101) nq-nqfirst+1,ngroupt,fname(1:len_trim(fname))
101            format(1x,'Writing ',i6,' points in ',i3,' groups to ',a)
            open(10,file=fname,status='unknown',form='formatted',recl=lrec)
                  lrecf=lrec-3
                  write(10,99) fname(1:index(fname,'   '))
                  write(10,99) info(1:len_trim(info))
                  write(10,100) ngroupt,nq
                  write(10,100) nskip
                  write(10,99) uname(1:len_trim(uname))
                  write(10,99) sttime(1:len_trim(sttime))
                  write(10,99) baddetfname(1:len_trim(baddetfname))
                  write(10,99) groupsfilefname(1:len_trim(groupsfilefname))
                  write(10,100) xcode
                  write(10,100) ycode
                  write(10,102) lenin
99     format('# ',a)
100    format('#',2(1x,i5))
102    format('#',6(1x,e14.7))
                  xo=0.0
!c
! write out group secondary flight paths
!c
                  write(10,formatout) '#',xo,(lentemp(ir),xo,ir=iwr1,iwr2)
!c103      format(a1,037(1x,e14.7))
!c
! write out group scattering angles
!c
                  write(10,formatout) '#',xo,(tthtemp(ir),xo,ir=iwr1,iwr2)
!c
! write out group azimuthal angles
!c
                  write(10,formatout) '#',xo,(phitemp(ir),xo,ir=iwr1,iwr2)
!c
! write the data for this range of groups
!c
                  do iq=nqfirst,nq
                        write(10,formatout) ' ',qbound(iq),(aggdweights(iq,ir),aggeweights(iq,ir),ir=iwr1,iwr2)
                  end do
                  iwr1=iwr2+1
            close(10)
      end do
      return
      end

      subroutine w_int(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
        
        use run_par
        use calibration_routines
        
      integer, parameter :: nwrmax=100
      integer xcode,ycode            !used by Genie load command
      integer nskip      !how many lines to skip
      integer lrec,lrecf,ilen,iwr,iwr1,iwr2,ir,iq,it,nt,nblock     !record lengths
      real lentemp(*),tthtemp(*),phitemp(*),xo
      character*256 fname      !name of file to write data to
      character*256 baddetfname     !bad detector file full filename
      character*256 groupsfilefname !groups file full filename
      character*32 formatout
      character*2 blk
!c
! After the initial 3 lines nskip is the number of lines of info that need to be skipped
! to get to the actual information
!c
        nskip=10
!c
! set up number of output files to write to
!c
      nt=(ngroupt+nwrmax-1)/nwrmax
      ilen=index(fname,' ')-1
      iwr1=1
      do it=1,nt
!c
! number of elements to write for each x value
!c
            iwr=ngroupt-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
            iwr2=iwr1+iwr-1
!c
! setup output format specification
!c
            nblock=2*iwr+1
            write(formatout,332) '(a1,',nblock,'(1x,e14.7))'
332      format(a4,i3.3,a11)
            write(6,*) formatout
            lrec=nblock*15+1
            if(lrec.lt.82) lrec=82
            write(blk,331) it
331      format(i2.2)
            fname=fname(1:ilen)//blk(1:2)
            write(6,101) nq-nqfirst+1,ngroupt,fname(1:len_trim(fname))
101            format(1x,'Writing ',i6,' points in ',i3,' groups to ',a)
            open(10,file=fname,status='unknown',form='formatted',recl=lrec)
                  lrecf=lrec-3
                  write(10,99) fname(1:index(fname,'   '))
                  write(10,99) info(1:len(info))
                  write(10,100) ngroupt,nq
                  write(10,100) nskip
                  write(10,99) uname(1:len_trim(uname))
                  write(10,99) sttime(1:len_trim(sttime))
                  write(10,99) baddetfname(1:len_trim(baddetfname))
                  write(10,99) groupsfilefname(1:len_trim(groupsfilefname))
                  write(10,100) xcode
                  write(10,100) ycode
                  write(10,102) lenin
99     format('# ',a)
100    format('#',2(1x,i5))
102    format('#',6(1x,e14.7))
                  xo=0.0
!c
! write out group secondary flight paths
!c
                  write(10,formatout) '#',xo,(lentemp(ir),xo,ir=iwr1,iwr2)
!c103      format(a1,037(1x,e14.7))
!c
! write out group scattering angles
!c
                  write(10,formatout) '#',xo,(tthtemp(ir),xo,ir=iwr1,iwr2)
!c
! write out group azimuthal angles
!c
                  write(10,formatout) '#',xo,(phitemp(ir),xo,ir=iwr1,iwr2)
!c
! write the data for this range of groups
!c
                  do iq=nqfirst,nq
                        write(10,formatout) ' ',qbound(iq),(aggsweights(iq,ir),aggeweights(iq,ir),ir=iwr1,iwr2)
                  end do
                  iwr1=iwr2+1
            close(10)
      end do
      return
      end
      
      subroutine w_intw(fname,baddetfname,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
        
        use run_par
        use calibration_routines
        
      integer, parameter :: nwrmax=100
      integer xcode,ycode            !used by Genie load command
      integer nskip      !how many lines to skip
      integer lrec,lrecf,ilen,iwr,iwr1,iwr2,ir,iq,it,nt,nblock     !record lengths
      real lentemp(*),tthtemp(*),phitemp(*),xo
      character*256 fname      !name of file to write data to
      character*256 baddetfname     !bad detector file full filename
      character*256 groupsfilefname !groups file full filename
      character*32 formatout
      character*2 blk
!c
! After the initial 3 lines nskip is the number of lines of info that need to be skipped
! to get to the actual information
!c
        nskip=10
!c
! set up number of output files to write to
!c
      nt=(ngroupt+nwrmax-1)/nwrmax
      ilen=index(fname,' ')-1
      iwr1=1
      do it=1,nt
!c
! number of elements to write for each x value
!c
            iwr=ngroupt-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
            iwr2=iwr1+iwr-1
!c
! setup output format specification
!c
            nblock=2*iwr+1
            write(formatout,332) '(a1,',nblock,'(1x,e14.7))'
332      format(a4,i3.3,a11)
            write(6,*) formatout
            lrec=nblock*15+1
            if(lrec.lt.82) lrec=82
            write(blk,331) it
331      format(i2.2)
            fname=fname(1:ilen)//blk(1:2)
            write(6,101) nq-nqfirst+1,ngroupt,fname(1:len_trim(fname))
101            format(1x,'Writing ',i6,' points in ',i3,' groups to ',a)
            open(10,file=fname,status='unknown',form='formatted',recl=lrec)
                  lrecf=lrec-3
                  write(10,99) fname(1:index(fname,'   '))
                  write(10,99) info(1:len(info))
                  write(10,100) ngroupt,nq
                  write(10,100) nskip
                  write(10,99) uname(1:len_trim(uname))
                  write(10,99) sttime(1:len_trim(sttime))
                  write(10,99) baddetfname(1:len_trim(baddetfname))
                  write(10,99) groupsfilefname(1:len_trim(groupsfilefname))
                  write(10,100) xcode
                  write(10,100) ycode
                  write(10,102) lenin
99     format('# ',a)
100    format('#',2(1x,i5))
102    format('#',6(1x,e14.7))
                  xo=0.0
!c
! write out group secondary flight paths
!c
                  write(10,formatout) '#',xo,(lentemp(ir),xo,ir=iwr1,iwr2)
!c103      format(a1,037(1x,e14.7))
!c
! write out group scattering angles
!c
                  write(10,formatout) '#',xo,(tthtemp(ir),xo,ir=iwr1,iwr2)
!c
! write out group azimuthal angles
!c
                  write(10,formatout) '#',xo,(phitemp(ir),xo,ir=iwr1,iwr2)
!c
! write the data for this range of groups
!c
                  do iq=nqfirst,nq
                        write(10,formatout) ' ',qbound(iq),(aggweights(iq,ir),aggeweights(iq,ir),ir=iwr1,iwr2)
                  end do
                  iwr1=iwr2+1
            close(10)
      end do
      return
      end

END MODULE merge_routines
    
    
