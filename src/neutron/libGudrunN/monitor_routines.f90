!     
! File:   monitor.f90
! Author: aks45
!
! Created on 27 April 2012, 13:12
!

MODULE monitor_routines
    
    implicit none

    integer                             :: ninmon        !no. of spectra for incident monitor
    integer                             :: ntrmon        !no. of spectra for transmission monitor
    integer, dimension(:), allocatable  :: incid_mon    !spectrum no.s of incident monitor
    integer, dimension(:), allocatable  :: trans_mon    !spectrum no.s of transmission monitor
    integer                             :: ncoeff        !no. of polynomial coeffs. for smovan
    integer, dimension(2)               :: oldmondims,mondims !Arrays to store the dimensions of the monitor arrays
    integer                             :: amontsumsize   !Size of the array amontsum
    real, dimension(:), allocatable     :: waveimon        !incident monitor wavelengths
    real, dimension(:), allocatable     :: wavetmon        !transmission monitor wavelengths
    real                                :: fraci,fract        !fractions of monitor for background
    real                                :: stepr            !step size for Monte Carlo fitting
    real                                :: acceptr            !acceptance factor for MC fitting
    real                                :: moncons,monthick        !monitor constant and foil thickness
    real                                :: monmus,monmuc        !monitor foil mus and muc
    real                                :: detcons            !expected detector constant
    real, dimension(:), allocatable     :: edmr        !fitted detector constants
    real, dimension(:), allocatable     :: dfactor           !detector normalisation factors
    real, dimension(:,:), allocatable   :: aexp    !polynomial coefficients for smovan
    real                                :: monwave1,monwave2    !wavelength range for monitor integration
    real, dimension(:), allocatable     :: amontsum    !sum of monitors for all periods
    
    CONTAINS
            
!***********************************************************************************
!*
!*      get_mon.FOR
!*
!*      A K Soper, May 2001
!*
!*      Gets the monitor data channel by channel for each spectrum, adding the
!*      data for a series of runs
!*
!***********************************************************************************
    subroutine get_mon(ntype,monref)

        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines
        use local_data
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use get_data_routines
        use spike_routines
        use write_routines
        use bad_detectors
        
        IMPLICIT NONE
!c
!c internal variables
!c
        character(len=256) fname            !general purpose filename
        character(len=20) handle            !handle to refer to each raw file
        integer ntype                  !type of data. 1=vanadium
        integer monref            !1 = incident mon., 2 = transmission mon. 
        integer i,ic,is,isref,iff,il,id,ip      !internal indices
        integer j,jf,jl,jref,ind,ic1      !internal indices
        integer nfile                  !no. of files to sum
        character(len=256), dimension(:), allocatable ::      runno!(mruns)            !filenames of files to sum
        real, dimension(:), allocatable                 :: amonsum !Temporary sum of monitor counts 
        integer nperrq            !period number to access
        integer ndum1,ndum2,nuse      !dummy variables
        integer errcode            !errcode on exit from reading data
        real lentot                  !total flight path
        real wave1,wave2            !wavelength boundaries
        real rat,fir,sec            !temporary real value
        real tsum                  !temporary ratio sum
        real err,ratiosum                  !temporary error
        real muamphrs
        integer goodframes,npersum
        real sum                    !temporary sum
        
        !Reallocate the amontsum array if necessary
        if(.not.allocated(amontsum)) then
            call reallocate1d_r(amontsum,ntype)
            amontsumsize=size(amontsum)
        else
            amontsumsize=size(amontsum)
            if(ntype.gt.amontsumsize) call reallocate1d_r(amontsum,ntype)
        end if
!        write(6,*) 'get_mon> Size of amontsum: ',size(amontsum)
!c
!c set up the time channel boundaries
!c
        if(ntype.lt.3) then
            nchan=nchanv
            nchanb=nchanbv
            call check_and_allocate_local_data(nchanb)
            call reallocate1d_r(tcb,nchanb)
            call reallocate1d_r(tcw,nchan)
            do ic=1,nchanb
                  tcb(ic)=tcbv(ic)
            end do
            do ic=1,nchan
                  tcw(ic)=tcwv(ic)
            end do
        else
            nchan=nchans
            nchanb=nchanbs
            call check_and_allocate_local_data(nchanb)
            call reallocate1d_r(tcb,nchanb)
            call reallocate1d_r(tcw,nchan)
            do ic=1,nchanb
                  tcb(ic)=tcbs(ic)
            end do
            do ic=1,nchan
                  tcw(ic)=tcws(ic)
            end do
        endif
        call reallocate1d_r(waveimon,nchanb)
        call reallocate1d_r(wavetmon,nchanb)


!c For D4C there are several monitors. Since these do not need to be smoothed or stored
!c we can skip most of this routine altogether.

        if(index(inst,'D4C').gt.0) then

!c calculate the wavelength scales of incident monitor

!c incident monitor: total flight path

            id=detno(incid_mon(1))
            lentot=lenin+lendet(id)
            rat=0.0039554/lentot
            wave1=rat*(tcb(1)+deltat(id))      
            waveimon(1)=wave1

!c set up average wavelength and width of each channel

            do ic=2,nchanb
                  ic1=ic-1
                  wave2=rat*(tcb(ic)+deltat(id))
                  waveimon(ic)=wave2
                  wave1=wave2
            end do
            amontsum(ntype)=1.0
            return

        endif

!c
!c first assign the run numbers
!c
        if(ntype.eq.1) then
            nfile=nfilev
            call reallocate1d_c(runno,len(runno),nfile)
            do i=1,nfile
                  runno(i)=runv(i)
            end do
            nperrq=nperv
        else if(ntype.eq.2) then
            nfile=nfilebv
            call reallocate1d_c(runno,len(runno),nfile)
            do i=1,nfile
                  runno(i)=runbv(i)
            end do
            nperrq=nperbv
        else if(ntype.eq.3) then
            nfile=nfilebs
            call reallocate1d_c(runno,len(runno),nfile)
            do i=1,nfile
                  runno(i)=runbs(i)
            end do
            nperrq=nperbs
        else 
            ind=ntype-3
            nfile=nfiles(ind)
            call reallocate1d_c(runno,len(runno),nfile)
            do i=1,nfile
                  runno(i)=runs(i,ind)
            end do
            nperrq=npers(ind)
        end if
!Get the number of periods for the first run. It is assumed the remaining runs have the same number of periods
        call get_periods(runno(1),npersum)
        !Set up and zero the monitor sums.
        call reallocate1d_r(amonsum,npersum)
        do ip=1,npersum
            amonsum(ip)=0.0
        end do
!        write(6,*) 'get_mon> Size of amonsum: ',npersum,size(amonsum)
!c
!c define total no. of channels being read. For monitors we need to retrieve the data one spectrum at at time
!c
        nspecrq=1
        nchanrq=nspecrq*nchan
        oldcountsdims=0
        oldtcountsdims=0
        countsdims(1)=nchanb
        countsdims(2)=nspecrq
        countsdims(3)=1
        tcountsdims(1)=countsdims(1)
        tcountsdims(2)=countsdims(2)
        call reallocate2d_r(tcounts,nchanb,nspecrq)

!Step through all the periods!c
!c step through requested files, adding to accumulator errcount if there
!c is more than one spectrum and more than one file
!c
        do iff=1,nfile
            do is=1,nspecproc
                if(specproc(is).ge.0) then
                    if(monref.eq.1) then
                        do ip=1,npersum
                            isref=nspecb*(ip-1)+specproc(is)
                            nspecf=isref
                            nspecl=isref
                            call get_data(runno(iff),nperrq)
!c
!c check for spikes
!c
                            if(spikedev.gt.0.0) call find_spike()

                            !Save the required period and sum the counts
                            if(ip.eq.nperrq) then
                                do ic=1,nchan
                                    amonsum(ip)=amonsum(ip)+tcounts(ic,1)
                                    detcount(ic)=tcounts(ic,1)/tcw(ic)
                                end do
                            else
                                !Otherwise, just sum the counts
                                do ic=1,nchan
                                    amonsum(ip)=amonsum(ip)+tcounts(ic,1)
                                end do
                            end if
                        end do
                    else
                        isref=nspecb*(nperrq-1)+specproc(is)
                        nspecf=isref
                        nspecl=isref
                        call get_data(runno(iff),nperrq)
!c
!c check for spikes
!c
                        call find_spike()

                        do ic=1,nchan
                            detcount(ic)=tcounts(ic,1)/tcw(ic)
                        end do
                    end if
                else
                    if(specproc(is).eq.-1) then
!c Get the number of microamp hours for this run
                        call get_muamphrs(runno(iff),muamphrs)
                        do ic=1,nchan
                            detcount(ic)=muamphrs
                        end do
                    else
                        call get_goodframes(runno(iff),goodframes)
                        do ic=1,nchan
                            detcount(ic)=real(goodframes)
                        end do
                    endif
                    !Summing the monitor counts is pointless in this case
                end if
!c put data into store
                if(is.eq.1.and.iff.eq.1) then
                    do ic=1,nchan
                        errcount(ic)=detcount(ic)
                    end do
                else
                    do ic=1,nchan
                        errcount(ic)=errcount(ic)+detcount(ic)
                    end do
                endif
            end do
        end do
        errcount(nchanb)=0.0
!c
!c save data to appropriate array
!c
        if(monref.eq.2) then
!c
!c calculate the wavelength scales of transmission monitors
!c
!c
!c transmission monitor: total flight path
!c
            id=detno(abs(trans_mon(1)))  !Make sure the selected spectrum exists, even if monitor is muamphrs or goodframes
            lentot=lenin+lendet(id)
            write(6,*) 'Trans> ',trans_mon(1),id,specno(id),lentot
            rat=0.0039554/lentot
            wave1=rat*(tcb(1)+deltat(id))      
            wavetmon(1)=wave1
!c
!c set up average wavelength and width of each channel
!c
            do ic=2,nchanb
                ic1=ic-1
                wave2=rat*(tcb(ic)+deltat(id))
                wavetmon(ic)=wave2
                wave1=wave2
            end do
            if(ntype.eq.1) then
                call reallocate1d_r(smovantrans,nchanb)
                do ic=1,nchanb
                    smovantrans(ic)=errcount(ic)
                end do      
            else if(ntype.eq.2) then
                call reallocate1d_r(smobvantrans,nchanb)
                do ic=1,nchanb
                    smobvantrans(ic)=errcount(ic)
                end do      
            else if(ntype.eq.3) then
                call reallocate1d_r(smobsamtrans,nchanb)
                do ic=1,nchanb
                    smobsamtrans(ic)=errcount(ic)
                end do      
            else 
                ind=ntype-3
                call reallocate2d_r(smosamtrans,nchanb,ind)
                do ic=1,nchanb
                    smosamtrans(ic,ind)=errcount(ic)
                end do
            endif
            fname=runno(1)
            call change_ext(fname,'rawtrans')
            do ic=1,nchanb
                detcount(ic)=errcount(ic)
                errcount(ic)=0.0
            end do
            call w_diag_file(fname,nchanb,wavetmon,detcount,errcount)
        else
!c
!c calculate the wavelength scales of incident monitors
!c
!c
!c incident monitor: total flight path
!c
            id=detno(abs(incid_mon(1))) !Make sure the selected spectrum exists, even if monitor is muamphrs or goodframes
            lentot=lenin+lendet(id)
            write(6,*) 'Incident> ',incid_mon(1),id,specno(id),lentot
            rat=0.0039554/lentot
            wave1=rat*(tcb(1)+deltat(id))      
            waveimon(1)=wave1
!c
!c set up average wavelength and width of each channel
!c
            do ic=2,nchanb
                  ic1=ic-1
                  wave2=rat*(tcb(ic)+deltat(id))
                  waveimon(ic)=wave2
                  wave1=wave2
            end do
            if(ntype.eq.1) then
                call reallocate1d_r(smovanmon,nchanb)
                do ic=1,nchanb
                    smovanmon(ic)=errcount(ic)
                end do      
            else if(ntype.eq.2) then
                call reallocate1d_r(smobvanmon,nchanb)
                do ic=1,nchanb
                    smobvanmon(ic)=errcount(ic)
                end do      
            else if(ntype.eq.3) then
                call reallocate1d_r(smobsammon,nchanb)
                do ic=1,nchanb
                    smobsammon(ic)=errcount(ic)
                end do      
            else 
                ind=ntype-3
                call reallocate2d_r(smosammon,nchanb,ind)
                do ic=1,nchanb
                    smosammon(ic,ind)=errcount(ic)
                end do
            endif
            fname=runno(1)
            call change_ext(fname,'rawmon')
            do ic=1,nchanb
                detcount(ic)=errcount(ic)
                errcount(ic)=0.0
            end do
            call w_diag_file(fname,nchanb,waveimon,detcount,errcount)
            !Finally output the relative monitor sum for the requested period
            sum=0.0
            do ip=1,npersum
                sum=sum+amonsum(ip)
            end do
!            write(6,*) 'get_mon> sum: ',sum,npersum,ntype,nperrq,size(amontsum),size(amonsum),amonsum(nperrq)/sum
            if(sum.gt.0.0) then
                ratiosum=amonsum(nperrq)/sum
            else
                ratiosum=1.0
            end if
            amontsum(ntype)=ratiosum
            write(6,*) 'get_mon> Relative monitor sums ',nperrq,amontsum(ntype)
        endif
        return
    end subroutine get_mon

!***********************************************************************************
!*
!*	form_smoo_mon.FOR
!*
!*	A K Soper, May 2001
!*
!*	Forms a smoothed incident monitor spectrum.
!*	Current version simply smooths the incident monitor using a 3-point
!*	 routine and uses this smoothed version.
!*
!***********************************************************************************
    subroutine form_smoo_mon(ntype,nsmoo,wavemin,wavemax)
	
        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines
        use local_data
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use write_routines
        use corrections_routines
        
        IMPLICIT NONE

!c
!c internal variables
!c
	character*256 fname		!name of file to write data to
	integer ntype			!type of data. 1=vanadium
	integer nsmoo			!maximum number of 3-point smoothings
	integer nsmoo1		!number of 3-point smoothings used
	integer i,ic,ic1,is,if,il	!internal indices
	integer j,jf,jl,jref,ind	!internal indices
	integer id			!stores the detector number of monitor
	character*256 run			!filename for the file to be smoothed
	integer naccept		!1 if spectrum fit is successful, else 0
	real, dimension(:), allocatable :: moncount		!monitor counts in this array
	real lentot			!total flight path for monitor
	real rat,err,ratio			!temporary real values
	real wave1,wave2		!temporary wavelengths
	real wavemin,wavemax		!minimum and maximum wavelengths to use

!c For D4C there are several monitors. Since these do not need to be smoothed or stored
!c we can skip this routine altogether.

	if(index(inst,'D4C').gt.0) then
		return
	endif
! Allocate some space
        call reallocate1d_r(moncount,nchanb)
        call check_and_allocate_local_data(nchanb)
!c
!c calculate counts per bin width and get the monitor wavelength boundaries
!c
	if(ntype.eq.1) then
		run=runv(1)
		do ic=1,nchan
			moncount(ic)=smovanmon(ic)
		end do
	else if(ntype.eq.2) then
		run=runbv(1)
		do ic=1,nchan
			moncount(ic)=smobvanmon(ic)
		end do
	else if(ntype.eq.3) then
		run=runbs(1)
		do ic=1,nchan
			moncount(ic)=smobsammon(ic)
		end do
	else 
		ind=ntype-3
		run=runs(1,ind)
		do ic=1,nchan
			moncount(ic)=smosammon(ic,ind)
		end do
	endif
!c
!c set up the time channel boundaries
!c
	if(ntype.lt.3) then
		nchan=nchanv
		nchanb=nchanbv
		do ic=1,nchanb
                    tcb(ic)=tcbv(ic)
		end do
                do ic=1,nchan
                    tcw(ic)=tcwv(ic)
                end do
	else
		nchan=nchans
		nchanb=nchanbs
		do ic=1,nchanb
                    tcb(ic)=tcbs(ic)
		end do
                do ic=1,nchan
                    tcw(ic)=tcws(ic)
                end do
	endif
!c
!c set up the weights array
!c
	do ic=1,nchan
		errcount(ic)=moncount(ic)/tcw(ic)
	end do
!c
!c smooth these data
!c
!c	call smoov(nsmoo,nchan,moncount,detcount,1,nchan)
	call smooe(1,nchan,nchan,nsmoo1,moncount,detcount,errcount,ratio,1.0)
	write(6,*) 'form_smoo_mon> ',nsmoo1,' smoothings used for monitor in run ',run(1:len_trim(run)),' with chisq ratio ',ratio
	fname=run
        call change_ext(fname,'smomon')
        call w_diag_file(fname,nchanb,waveimon,detcount,errcount)
!c
!c save the smoothed monitor counts
!c
	if(ntype.eq.1) then
		do ic=1,nchan
			smovanmon(ic)=detcount(ic)
		end do
	else if(ntype.eq.2) then
		do ic=1,nchan
			smobvanmon(ic)=detcount(ic)
		end do
	else if(ntype.eq.3) then
		do ic=1,nchan
			smobsammon(ic)=detcount(ic)
		end do
	else 
		ind=ntype-3
		do ic=1,nchan
			smosammon(ic,ind)=detcount(ic)
		end do
	endif
101	format(3(1x,e13.6))
	return	
    END subroutine form_smoo_mon
    
END MODULE monitor_routines
