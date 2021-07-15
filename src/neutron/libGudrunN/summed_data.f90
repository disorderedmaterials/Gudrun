!     
! File:   summed_data.f90
! Author: aks45
!
! Created on 15 November 2013, 19:53
!

MODULE summed_data
    
    use reallocation_routines
    use run_par
    use local_data
    use get_data_routines
    use spec_van
    use spec_bakv
    use spec_baks
    use spec_sam
    
    implicit none
    
    real, dimension(:,:), allocatable               :: van_summed_data
    real, dimension(:,:), allocatable               :: vanback_summed_data
    real, dimension(:,:), allocatable               :: samback_summed_data
    real, dimension(:,:,:), allocatable             :: sam_summed_data
    integer                                         :: van_goodframes
    integer                                         :: vanback_goodframes
    integer                                         :: samback_goodframes
    integer                                         :: sam_goodframes
    
    CONTAINS
    
    subroutine get_summed_data(ntype)
        
        
        integer, intent(in)                                 :: ntype
        character(len=256), dimension(:), allocatable       :: runno            !filename to obtain calibration
        integer                                             :: nrunno,nperrq,i,is,ic,iff,ind,nfile,nchar
        integer                                             :: goodframes,ngoodframes
        real, dimension(:,:), allocatable                   :: summed_data
        character(len=256) :: fname        !name of file to read
        
!c
!c first assign the run numbers
!c
        nchar=256
        if(ntype.eq.1) then
            nfile=nfilev
            call reallocate1d_c(runno,nchar,nfile)
            do i=1,nfile
                  runno(i)=runv(i)
            end do
            nperrq=nperv
        else if(ntype.eq.2) then
            nfile=nfilebv
            call reallocate1d_c(runno,nchar,nfile)
            do i=1,nfile
                  runno(i)=runbv(i)
            end do
            nperrq=nperbv
        else if(ntype.eq.3) then
            nfile=nfilebs
            call reallocate1d_c(runno,nchar,nfile)
            do i=1,nfile
                  runno(i)=runbs(i)
            end do
            nperrq=nperbs
        else 
            ind=ntype-3
            nfile=nfiles(ind)
            call reallocate1d_c(runno,nchar,nfile)
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
            do ic=1,nchanb
                  tcb(ic)=tcbv(ic)
                  tcw(ic)=tcwv(ic)
            end do
            nspch=nspchv
        else
            nchan=nchans
            nchanb=nchanbs
            do ic=1,nchanb
                  tcb(ic)=tcbs(ic)
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
        end if
        !Get the temporary store
        call reallocate2d_r(tcounts,nchanb,nspec)
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
!c total number of spectra being read
!c
        nspecrq=nspec
!c
!c total number of channels requested
!c
        nchanrq=nspecrq*nchanb
!c
!c step through the files summing the data into a single array
!c
        ngoodframes=0
        do iff=1,nfile
!c
!c get the data
!c
            write(6,*) 'get_summed_data> Reading data file: ',runno(iff)(1:len_trim(runno(iff))),' period ',nperrq
            call get_data(runno(iff),nperrq)
            call get_goodframes(runno(iff),goodframes)
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
                ngoodframes=ngoodframes+goodframes
                do is=1,nspec
                    do ic=1,nchan
                        summed_data(ic,is)=summed_data(ic,is)+tcounts(ic,is)
                    end do
                    summed_data(nchanb,is)=0 !Last value set to zero
                end do
            end if
        end do
        !Save to the appropriate array
        if(ntype.eq.1) then
            do is=1,nspec
                do ic=1,nchanb
                    van_summed_data(ic,is)=summed_data(ic,is)
                end do
            end do
        else if(ntype.eq.2) then
            do is=1,nspec
                do ic=1,nchanb
                    vanback_summed_data(ic,is)=summed_data(ic,is)
                end do
            end do
        else if(ntype.eq.3) then
            do is=1,nspec
                do ic=1,nchanb
                    samback_summed_data(ic,is)=summed_data(ic,is)
                end do
            end do
        else 
            ind=ntype-3
            do is=1,nspec
                do ic=1,nchanb
                    sam_summed_data(ic,is,ind)=summed_data(ic,is)
                end do
            end do
        end if
        if(allocated(tcounts)) deallocate(tcounts,summed_data)
        return
        
    end subroutine get_summed_data
    
END MODULE summed_data
