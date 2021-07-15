!     
! File:   module_sums.f90
! Author: aks45
!
! Created on 15 November 2013, 20:50
!

MODULE module_sums
    
    implicit none
    
    real, dimension(:,:,:), allocatable         :: modulecnts     !sum of counts in each channel and module, for each sample type.
    real, dimension(:,:), allocatable           :: sol!sums of counts per pulse in whole detector array and in modules
    integer, dimension(:), allocatable          :: nspecpermodule
    integer                                     :: nmod1                  !reference for sum of counts in all modules
    
    CONTAINS
    
!***********************************************************************************
!*
!*      get_module_sum.FOR
!*
!*      Copyright A K Soper, November 2003
!*       derived from the ISIS GET routines by Freddie Akeroyd
!*
!*      Gets the data channel by channel for each spectrum, adding the
!*      data for a series of runs, and combining the spectra into the DAE input 
!*     'modules'.
!*
!*      It is assumed that the monitor sums have already been read in with get_mon 
!*      so that the total monitor counts in each period are available
!*
!***********************************************************************************
    subroutine get_module_sum(ntype)

    use reallocation_routines
    use run_par
    use calibration_routines
    use local_data
    use monitor_routines
    use groups_routines
    use summed_data_routines
    use deadtime_routines
    use bad_detectors
    
!c
!c internal variables
!c
        character(len=256) fname
        character*32 formatout
        character*1 number
        integer ntype                  !type of data. 1=vanadium
        
        integer i,ic,is,iff,il,icref      !internal indices
        integer id,im,iflag,it,niter
        integer j,jf,jl,jref,ind,lrec      !index range limits
        integer nblock,nfirst,nlast,nmodtime
        real rat,rat1,rat2,ratsave,rat1save            !temporary real value
        real tsum                  !temporary ratio sum
        real err                  !temporary error 
        real wave,wfac
        real lentot
        real runfac            !run factor for a given run number
        real get_runfactor      !function which returns the factor for a specified run number
        real sumdae,summodule
        real sumchan,summod,sumdiff,sumdiffold,facdae,solnew,diff
        real prod,proddae,prodmax,tcbmax             !temporary product
        real, dimension(:), allocatable     :: runtime !Use for D4C to store the total runtime for a set of modules
        integer immax

        nmod1=nmod+1
        if(nmod.le.0) then
            write(6,*) 'get_module_sum> No modules defined: cannot proceed!'
            stop
        endif
        if(nchan.le.0) then
            write(6,*) 'get_module_sum> No time channels defined: cannot proceed!'
            stop
        end if
        if(nspec.le.0) then
            write(6,*) 'get_module_sum> No spectra defined: cannot proceed!'
            stop
        end if
        !Get the summed data, number of good frames, and time channel boundaries
        call get_summed_data(ntype)
        write(6,*) 'get_module_sum> Total of ',ngoodframes,' goodframes read from ',runno(1)(1:len_trim(runno(1)))
        call reallocate3d_r(modulecnts,nchanb,nmod1,ntype)
        call reallocate2d_r(sol,nchanb,nmod1)
        call reallocate1d_i(nspecpermodule,nmod)
        if(index(inst,'D4C').gt.0) then
            call reallocate1d_r(runtime,nchanb)
        end if
!c
!c set the module sums for this ntype to 1 or zero, depending on whether deadtimes have been defined
!c
        if(adeadtime.eq.0.0.and.mdeadtime.eq.0.0.and.ddeadtime.eq.0.0) then
            do im=1,nmod1
                do ic=1,nchan
                    modulecnts(ic,im,ntype)=1.0
                    sol(ic,im)=0.0
                end do
            end do
            return
        else
            do im=1,nmod1
                do ic=1,nchan
                    modulecnts(ic,im,ntype)=0.0
                    sol(ic,im)=0.0
                end do
            end do
        end if            
!c
!c Zero the number of spectra per module
!c
        do im=1,nmod
            nspecpermodule(im)=0
        end do
!c
!c Count the number of spectra which will contribute to each module count
!c
        do is=1,nspec
!c
!c set up the module number for this spectrum
!c
            im=imodule(detno(is))
            nspecpermodule(im)=nspecpermodule(im)+1
        end do
!c 
!c store the data in the appropriate module counts array
!c
        do is=1,nspec

!c Do not include this detector if it is bad

            if(ibad(is).eq.0) then
!c
!c set up the module number for this spectrum
!c
                im=imodule(detno(is))
                do ic=1,nchan
!c
!c save the counts in the appropriate array
!c
                    modulecnts(ic,im,ntype)=modulecnts(ic,im,ntype)+summed_data(ic,is)
                end do
            end if

        end do
!c
!c calculate the number of good frames for this period
!c
        ngoodframes=ngoodframes*amontsum(ntype)
        write(6,*) 'get_module_sum> ' &
        ,ngoodframes,' good frames found for run ' &
        ,runno(1)(1:index(runno(1),'  ')-1),' period ',nperrq
!c
!c Now estimate the module deadtime correction
!c
        sumdae=0.0
        do im=1,nmod
            iflag=0
            if(index(inst,'D4C').gt.0) then !For D4C need to normalise the module counts to the total run time
                nmodtime=crateno(im)+1 !crateno stores the module that stores the monitor, one before the module storing the run time
                call get_summed_data_spec(ntype,nmodtime,runtime)
                rat1=1.0/(runtime(1)*1000.0) !Convert time in milliseconds to microseconds
            else if(ngoodframes.gt.0) then
                rat1=1.0/real(ngoodframes)
            else
                rat1=0.0
            end if
            if(nspecpermodule(im).gt.0) then
                rat2=1.0/real(nspecpermodule(im))
            else
                rat2=0.0
            endif
            summodule=0.0
            do ic=1,nchan

!c Normalise counts to number of good frames

                rat=modulecnts(ic,im,ntype)*rat1

!c Accumulate DAE and module sums

                sumdae=sumdae+rat
                summodule=summodule+rat
!c
!c Count rate in channel
!c
                rat=rat/tcw(ic)
!c
!c sum the counts in all modules
!c
                modulecnts(ic,nmod1,ntype)=modulecnts(ic,nmod1,ntype)+rat
!c Average countrate per detector in each module as a function of time of flight
                modulecnts(ic,im,ntype)=rat*rat2
            end do
!c Pulse averaged count rate per detector in this module
            modulecnts(nchanb,im,ntype)=summodule*rat2/(tcb(nchanb)-tcb(1))
        end do
!c Pulse averaged count rate in whole detector array
        modulecnts(nchanb,nmod1,ntype)=sumdae/(tcb(nchanb)-tcb(1))
!c Find a solution for the dead time correction
        nfirst=1
        nlast=nmod
        niter=1
!c Set up the initial solution
        do ic=1,nchan
            do im=1,nmod
                sol(ic,im)=modulecnts(ic,im,ntype)
            end do
        end do
!c Sum the deadtime corrected counts in all the DAE
        sumchan=0.0
        do ic=1,nchan
            summod=0.0
            do im=nfirst,nlast
                summod=summod+sol(ic,im)*nspecpermodule(im)/(1.0+sol(ic,im)*ddeadtimem(im))
            end do
            sumchan=sumchan+summod*tcw(ic)
        end do
!c Overall data rate entering the DAE
        write(6,103) sumchan
103     format(/' get_module_sum> Overall DAE count rate (cts/musec): ',e13.6)
        immax=0
        tcbmax=0
        prodmax=0
        do ic=1,nchan
!c Ignore channels where the deadtime correction doesn't work
            proddae=modulecnts(ic,nmod1,ntype)*adeadtime
            if(proddae.ge.0.99) then
                tcbmax=tcb(ic)
                proddae=0.0
                write(6,156) adeadtime,modulecnts(ic,nmod1,ntype),ic
156   format(/' get_module_sum> DAE deadtime: ',f10.5,' is too large for observed count rate',e13.6,' at time channel ',i4)
                stop
            endif
!c DAE factor for this time channel
            facdae=1.0/(1.0-proddae)
            if(tcb(ic).gt.tcbmax) then
                summod = 0.0
                do im=nfirst,nlast
                    prod = modulecnts(ic,im,ntype)*ddeadtimem(im)
                    if(prod.ge.0.99) then
                        immax=im
                        prodmax=prod
                        tcbmax=tcb(ic)
                        prod=0.0
                        write(6,157) ddeadtimem(im),modulecnts(ic,im,ntype),ic,im
157   format(/' get_module_sum> Detector deadtime: ',f10.5 &
        ,' is too large for observed count rate',e13.6,' at time channel ',i4,' in module ',i4)
                        stop
                    endif
                    sol(ic,im)=modulecnts(ic,im,ntype)*facdae/(1.0-prod)
                    summod=summod+sol(ic,im)
                end do
!c Calculate the overall countrate in the DAE for this time channel
                sol(ic,nmod1)=summod
            endif
        end do
        write(6,*) immax,ntype,tcbmax,prodmax
        if(it.eq.1) sumdiffold=sumdiff
!c Finally convert the result to a ratio of solution divided by input data
        do ic=1,nchan
            do im=1,nmod
                if(modulecnts(ic,im,ntype).gt.0.0) then
                    modulecnts(ic,im,ntype)=sol(ic,im)/modulecnts(ic,im,ntype)
                else
                    modulecnts(ic,im,ntype)=1.0
                endif
            end do
        end do
!c
!c write diagnostic file
!c
        if(ntype.eq.1) fname='vanadium.module'
        if(ntype.eq.2) fname='vanadium_back.module'
        if(ntype.eq.3) fname='sample_back.module'
        if(ntype.eq.4) fname='sample.module'
        if(ntype.gt.4) then
            write(number,200) ntype-4
200         format(i1)
            fname='container_'//number(1:1)//'.module'
        endif
!c
!c setup output format specification
!c
        nblock=nmod+2
        write(formatout,332) '(a1,',nblock,'(1x,e14.7))'
332     format(a4,i3.3,a11)
        lrec=nblock*15+1
        if(lrec.lt.82) lrec=82
        open(10,file=fname,status='unknown',form='formatted',recl=lrec)
        write(10,99) fname
        write(10,99)
        write(10,99)
        write(10,99)
99      format('# ',a)
        do ic=1,nchanb
            write(10,formatout) '  ',tcb(ic),(modulecnts(ic,im,ntype),im=1,nmod1)
        end do      
        close(10)
!c      do im=1,nmod
!c         write(6,*) im,nspecpermodule(im)
!c      end do
        if(allocated(modulecnts)) write(6,*) 'get_module_sum> modulecnts(1,nmod,1:ntype)',(modulecnts(1,nmod,i),i=1,ntype)
        return
        
    END subroutine get_module_sum
    
!***********************************************************************************
!*
!*	get_sum.FOR
!*
!*	Copyright A K Soper, January 2003
!* 	derived from the ISIS GET routines by Freddie Akeroyd
!*
!*	Gets the data channel by channel for each spectrum, adding the
!*	data for a series of runs
!*
!***********************************************************************************
    subroutine get_sum(ntype,nspecwrt)

        use reallocation_routines
        use run_par
        use local_data
        use summed_data_routines
        use inputfilestrings
        use calibration_routines
        use groups_routines
        use monitor_routines
        use write_routines
        use spec_van
        use spec_bakv
        use spec_baks
        use spec_sam
        use deadtime_routines
!c
!c internal variables
!c
	character(len=256) fname		!general purpose filename
	integer ntype			!type of data. 1=vanadium
	integer i,ic,is,iff,il,icref	!internal indices
	integer id
	integer j,jf,jl,jref,ind	!index range limits
	integer nfile			!no. of files to sum
	integer nspecwrt		!spectrum no. for diagnostic file
	integer isref			!Gudrun spectrum number
	real rat,fir,sec			!temporary real value
	real tsum			!temporary ratio sum
	real err			!temporary error 
	real wave,wfac
	real lentot
	real runfac		!run factor for a given run number
!c
!c Set up the time channel boundaries
!c
	if(ntype.lt.3) then
            nchan=nchanv
            nchanb=nchanbv
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
            call reallocate1d_r(tcb,nchanb)
            call reallocate1d_r(tcw,nchan)
            do ic=1,nchanb
		tcb(ic)=tcbs(ic)
            end do
            do ic=1,nchan
                tcw(ic)=tcws(ic)
            end do
	endif
!
! Get the relevant data and allocate arrays
!
        call reallocate2d_r(tcounts,nchanb,nspec)
        if(ntype.eq.1) then
            fname=runv(1)
            ngoodframes=van_goodframes
            do is=1,nspec
                do ic=1,nchan
                    tcounts(ic,is)=van_summed_data(ic,is)
                end do
            end do
            call reallocate2d_r(smovandet,nchanb,nspecproc)
            call reallocate2d_r(errvandet,nchanb,nspecproc)
        else if(ntype.eq.2) then
            fname=runbv(1)
            ngoodframes=vanback_goodframes
            do is=1,nspec
                do ic=1,nchan
                    tcounts(ic,is)=vanback_summed_data(ic,is)
                end do
            end do
            call reallocate2d_r(normbakdet,nchanb,nspecproc)
            call reallocate2d_r(errbakdet,nchanb,nspecproc)
        else if(ntype.eq.3) then
            fname=runbs(1)
            ngoodframes=samback_goodframes
            do is=1,nspec
                do ic=1,nchan
                    tcounts(ic,is)=samback_summed_data(ic,is)
                end do
            end do
            call reallocate2d_r(normbakdet,nchanb,nspecproc)
            call reallocate2d_r(errbakdet,nchanb,nspecproc)
	else 
            ind=ntype-3
            fname=runs(1,ind)
            ngoodframes=sam_goodframes(ind)
            do is=1,nspec
                do ic=1,nchan
                    tcounts(ic,is)=sam_summed_data(ic,is,ind)
                end do
            end do
            call reallocate3d_r(normsamdet,nchanb,nspecproc,ind)
            call reallocate2d_r(normsamdetnosub,nchanb,nspecproc)
            call reallocate3d_r(errsamdet,nchanb,nspecproc,ind)
	end if
!
!Do deadtime correction
!
        if(ngoodframes.gt.0) call do_ddeadtime_corr(nspecwrt,ntype)
!c
!c save data to appropriate array
!c
	do is=1,nspecproc
            isref=specproc(is)
            if(ntype.eq.1) then
		do ic=1,nchan
                    smovandet(ic,is)=tcounts(ic,isref)/tcw(ic)
                    errvandet(ic,is)=smovandet(ic,is)/tcw(ic)
		end do	
            else if(ntype.eq.2) then
		do ic=1,nchan
                    normbakdet(ic,is)=tcounts(ic,isref)/tcw(ic)
                    errbakdet(ic,is)=normbakdet(ic,is)/tcw(ic)
		end do	
            else if(ntype.eq.3) then
		do ic=1,nchan
                    normbakdet(ic,is)=tcounts(ic,isref)/tcw(ic)
                    errbakdet(ic,is)=normbakdet(ic,is)/tcw(ic)
 		end do	
            else 
		ind=ntype-3 
		do ic=1,nchan
                    normsamdet(ic,is,ind)=tcounts(ic,isref)/tcw(ic)
                    errsamdet(ic,is,ind)=normsamdet(ic,is,ind)/tcw(ic)
                end do	
            endif
!c
!c write data if diagnostic file is required for this spectrum
!c
            if(isref.eq.nspecwrt) then
                call change_ext(fname,'cnt')
		id=detno(isref)
		lentot=lenin+lendet(id)
		wfac=0.0039554/lentot
                call reallocate1d_r(wavebound,nchanb)
                call reallocate1d_r(detcount,nchanb)
                call reallocate1d_r(errcount,nchanb)
!                write(6,*) 'get_sum> ngoodframes: ',ngoodframes,tcw(nchan)
		do ic=1,nchanb
                    wavebound(ic)=wfac*(tcb(ic)+deltat(id))
                    if(ic.lt.nchanb) then
!                        rat=real(tcounts(ic,isref))
                        rat=real(tcounts(ic,isref))/(ngoodframes*tcw(ic))
                        detcount(ic)=rat
			errcount(ic)=rat/(ngoodframes*tcw(ic))
                    else
                        detcount(ic)=0.0
			errcount(ic)=0.0
                    endif
		end do
		call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
            endif
	end do
!
!c For D4C the monitor spectrum is the first spectrum in the series, provided
!c the number of spectra processed at one time, nproc = 9 * 64 +1. This value is
!c set in get_run_par_D4C

	if(index(inst,'D4C').gt.0) then
            icref=1
            if(ntype.eq.1) then
                 call reallocate1d_r(smovanmon,nchanb)
                 do ic=1,nchan
                     smovanmon(ic)=tcounts(icref,1)/tcw(ic)
                     icref=icref+1
                 end do
            else if(ntype.eq.2) then
                call reallocate1d_r(smobvanmon,nchanb)
                do ic=1,nchan
                    smobvanmon(ic)=tcounts(icref,1)/tcw(ic)
                    icref=icref+1
                end do
            else if(ntype.eq.3) then
                call reallocate1d_r(smobsammon,nchanb)
                do ic=1,nchan
                    smobsammon(ic)=tcounts(icref,1)/tcw(ic)
                    icref=icref+1
                end do
            else 
                ind=ntype-3
                call reallocate2d_r(smosammon,nchanb,ind)
                do ic=1,nchan
                    smosammon(ic,ind)=tcounts(icref,1)/tcw(ic)
                    icref=icref+1
                end do
            endif
!c Write diagnostic file
            call change_ext(fname,'rawmon')
            open(10,file=fname,status='unknown')
            icref=1
            call reallocate1d_r(wavebound,nchanb)
            call reallocate1d_r(detcount,nchanb)
            call reallocate1d_r(errcount,nchanb)
            do ic=1,nchanb
            	wavebound(ic)=waveimon(ic)
		if(ic.lt.nchanb) then
                    rat=real(tcounts(icref,1))/tcw(ic)
                    detcount(ic)=rat
                    errcount(ic)=rat/tcw(ic)
                    icref=icref+1
		else
                    detcount(ic)=0.0
                    errcount(ic)=0.0
		endif
            end do
            call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
	endif
	return	
    END subroutine get_sum

!***********************************************************************************
!*
!*	do_ddeadtime_corr.FOR
!*
!*	Copyright A K Soper, November 2003
!* 	derived from the ISIS GET routines by Freddie Akeroyd
!*
!*	Performs a simple scattering detector deadtime correction on the data in tcounts
!*	It is assumed the time channel widths, TCW, have already been defined.
!*
!***********************************************************************************
    subroutine do_ddeadtime_corr(nspecwrt,ntype)
            
        use run_par
        use local_data
        use reallocation_routines
        use calibration_routines
        use monitor_routines
        use groups_routines
        use summed_data_routines
        use deadtime_routines
!c
!c internal variables
!c
	integer ic,is,icref,ierr	!internal indices
	integer im,ntype,iflag,nspecwrt,isref	!module number and type of data file (vanadium, background, etc.)
	real rat,ratm,factor			!temporary real value
	real deadtot,deadcorr,deadcorrmax,get_datarate,deadcorrfac
!c
!c Step through the data, performing the appropriate correction
!c
	deadcorrmax=0.0
	do is=1,nspecproc
            isref=specproc(is)
            im=imodule(detno(isref))
!            write(6,*) 'do_ddeadtime_corr> ',is,isref,im,ntype,modulecnts(1,im,ntype)
            do ic=1,nchan
!c Get deadtime correction
                deadcorr=modulecnts(ic,im,ntype)
                if(deadcorr.gt.deadcorrmax) deadcorrmax=deadcorr
                tcounts(ic,isref)=tcounts(ic,isref)*deadcorr
            end do	
	end do
	write(6,*) 'do_ddeadtime_corr> Maximum deadtime correction: ',deadcorrmax,ntype
	return	

    END subroutine do_ddeadtime_corr

END MODULE module_sums

