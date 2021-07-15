!     
! File:   summed_data.f90
! Author: aks45
!
! Created on 15 November 2013, 19:53
!

MODULE summed_data_routines
    
    implicit none
    
    real, dimension(:,:), allocatable               :: van_summed_data
    real, dimension(:,:), allocatable               :: vanback_summed_data
    real, dimension(:,:), allocatable               :: samback_summed_data
    real, dimension(:,:,:), allocatable             :: sam_summed_data
    integer                                         :: van_goodframes
    integer                                         :: vanback_goodframes
    integer                                         :: samback_goodframes
    integer, dimension(:), allocatable              :: sam_goodframes
    real, dimension(:,:), allocatable               :: summed_data
    integer                                         :: ngoodframes,nrunno,nperrq
    character(len=256), dimension(:), allocatable   :: runno            !filename to obtain calibration
    
    CONTAINS
    
    subroutine get_summed_data(ntype)
        
        use reallocation_routines
        use run_par
        use local_data
        use get_data_routines
        use spec_van
        use spec_bakv
        use spec_baks
        use spec_sam
        use bad_detectors
        use spike_routines
        use runfactor_list
        
        integer, intent(in)                                 :: ntype
        integer                                             :: i,is,ic,iff,ind,nfile
        character(len=256) :: fname        !name of file to read
        integer                                             :: goodframes
        real                                                :: runfac
        
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
!                  write(6,*) 'get_summed_data> Sample background file: ',runno(i)(1:len_trim(runno(i)))
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
!c
!c now set up the time channel boundaries
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
            nspch=nspchv
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
            nspch=nspchs
        endif
        !Allocate the required array
        if(ntype.eq.1) then
            call reallocate2d_r(van_summed_data,nchanb,nspec)
        else if(ntype.eq.2) then
            call reallocate2d_r(vanback_summed_data,nchanb,nspec)
        else if(ntype.eq.3) then
            call reallocate2d_r(samback_summed_data,nchanb,nspec)
        else 
            ind=ntype-3
            call reallocate3d_r(sam_summed_data,nchanb,nspec,ind)
            call reallocate1d_i(sam_goodframes,ind)
        end if
        !Get the temporary store
        call reallocate2d_r(summed_data,nchanb,nspec)
!c
!c ensure period requested is sensible - if not use period 1
!c
        if(nperrq.le.0.or.nperrq.gt.nper) nperrq=1
!c
!c set up the spectrum numbers to read.
!c
        nspecf=1
        nspecl=nspec
!c
!c step through the files summing the data into a single array
!c
        ngoodframes=0
        do iff=1,nfile
!c
!c get the factor for this runno.
!c
            runfac=get_runfactor(runno(iff))
!c
!c get the data
!c
            write(6,*) 'get_summed_data> Reading data file: ',runno(iff)(1:len_trim(runno(iff))),' period ',nperrq,' factor ',runfac
            call get_data(runno(iff),nperrq)
!            write(6,*) 'get_summed_data> Read data'
            if (spikedev.gt.0.0) call find_spike()
!           write(6,*) 'get_summed_data> Checked for spikes'
            call get_goodframes(runno(iff),goodframes)
!            write(6,*) 'get_summed_data> Goodframes: ',goodframes
            if(runfac.ne.1.0) then
                if (iff.eq.1) then
                    !Initialise the sums
                    ngoodframes=goodframes
                    do is=1,nspec
                        do ic=1,nchan
                            summed_data(ic,is)=runfac*tcounts(ic,is)
                        end do
                        summed_data(nchanb,is)=0 !Last value set to zero
                    end do
                else
                    !Add to existing sum
                    if(index(inst,'D4C').eq.0) ngoodframes=ngoodframes+goodframes
                    do is=1,nspec
                        do ic=1,nchan
                            summed_data(ic,is)=summed_data(ic,is)+runfac*tcounts(ic,is)
                        end do
                        summed_data(nchanb,is)=0 !Last value set to zero
                    end do
                end if
            else
                if (iff.eq.1) then
                    !Initialise the sums
                    ngoodframes=goodframes
                    do is=1,nspec
                        do ic=1,nchan
                            summed_data(ic,is)=tcounts(ic,is)
                        end do
                        summed_data(nchanb,is)=0 !Last value set to zero
                    end do
                else
                    !Add to existing sum
                    if(index(inst,'D4C').eq.0) ngoodframes=ngoodframes+goodframes
                    do is=1,nspec
                        do ic=1,nchan
                            summed_data(ic,is)=summed_data(ic,is)+tcounts(ic,is)
                        end do
                        summed_data(nchanb,is)=0 !Last value set to zero
                    end do
                end if
            end if
        end do
        !Save to the appropriate array
        if(ntype.eq.1) then
            do is=1,nspec
                do ic=1,nchanb
                    van_summed_data(ic,is)=summed_data(ic,is)
                end do
                if(is.eq.2886.or.is.eq.2890) write(6,*) 'get_summed_data> ',is,summed_data(nchan,is)
            end do
            van_goodframes=ngoodframes
        else if(ntype.eq.2) then
            do is=1,nspec
                do ic=1,nchanb
                    vanback_summed_data(ic,is)=summed_data(ic,is)
                end do
                if(is.eq.2886.or.is.eq.2890) write(6,*) 'get_summed_data> ',is,summed_data(nchan,is)
            end do
            vanback_goodframes=ngoodframes
        else if(ntype.eq.3) then
            do is=1,nspec
                do ic=1,nchanb
                    samback_summed_data(ic,is)=summed_data(ic,is)
                end do
            end do
            samback_goodframes=ngoodframes
        else 
            ind=ntype-3
            do is=1,nspec
                do ic=1,nchanb
                    sam_summed_data(ic,is,ind)=summed_data(ic,is)
                end do
            end do
            sam_goodframes(ind)=ngoodframes
        end if
!        if(allocated(tcounts)) deallocate(tcounts)
        return
        
    end subroutine get_summed_data
    
!***********************************************************************************
!*
!*	get_summed_data_spec.f90
!*
!*	A K Soper, February 2014
!*
!*	puts a temporary array back into the spectrum array
!*
!* This is needed by D4C to get the correct monitor counts, since these do not appear in the usual arrays, normvandet, etc.
!*
!***********************************************************************************
    subroutine get_summed_data_spec(ntype,is,tarray)

        use run_par
!c
!c internal variables
!c
	integer ic,is,index		!internal indices
	integer ntype			!vanadium, sample, can, furnace, etc.
	real, dimension(:)              :: tarray
!c	write(6,*) nfilebv,runbv(1)
!c
!c get normalised data for this spectrum
!c
	if(ntype.eq.1) then
            do ic=1,nchanv
		tarray(ic)=van_summed_data(ic,is)
            end do
	else if(ntype.eq.2) then
            do ic=1,nchanv
		tarray(ic)=vanback_summed_data(ic,is)
            end do
	else if(ntype.eq.3) then
            do ic=1,nchans
		tarray(ic)=samback_summed_data(ic,is)
            end do
        else if(ntype.gt.3) then
            index=ntype-3
            do ic=1,nchans
            	tarray(ic)=sam_summed_data(ic,is,index)
            end do
	end if
	return	
    END subroutine get_summed_data_spec
    
END MODULE summed_data_routines
