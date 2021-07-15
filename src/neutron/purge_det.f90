!     
! File:   purge_det_1.f90
! Author: aks45
!
! Created on 27 April 2012, 13:23
!
!***********************************************************************************
!*
!*    Purge_det_NEW.FOR
!*
!*    COPYRIGHT A K Soper, JANUARY 2003
!*
!*    This is a stand-alone version of PURGED, reading data directly from 
!*  ISIS RAW files
!*
!*     Further modified February 2005 for input from GudrunGUI_1
!*
!********************************************************************************
MODULE spectrum_ratio

    implicit none
    
    real, dimension(:,:), allocatable  :: specrat               !integrated spectrum ratios
    real, dimension(:,:), allocatable  :: errrat                !rms deviation on ratio

    CONTAINS

    subroutine get_ratio(runs,nfiles,nperrq)
        
        use reallocation_routines
        use inputfilestrings
        use run_par
        use local_data
        use get_data_routines
        use spike_routines
        use monitor_routines
        
!***********************************************************************************
!*
!*        get_ratio.FOR
!*
!*        Copyright A K Soper, January 2003 
!*
!*        Gets the integrated ratio, spectrum by spectrum for two runs
!*
!***********************************************************************************
        character(len=256),dimension(:),intent(in)  :: runs     !Filenames
        integer, intent(in)                         :: nfiles   !No. of files to calculate ratios for
        integer, dimension(:), intent(inout)        :: nperrq   !period numbers to access
        !c internal variables
        character(len=256)                          :: fname                !name of file to write ratios to
        integer                                     :: nblock
        integer                                     :: i,is,j,isref,itref,ierr        !internal indices
        integer                                     :: iff,jf,jl        !index range limits
        real, dimension(:,:,:), allocatable         :: scounts        !stores counts to be ratioed
        real                                        :: rat,ratsum,eratsum,tsum,err  !temporary values
        integer                                     :: nc1,nc2,nspecwrt,ilen
        integer                                     :: ms1,ms2,mssum  !Used to store time information
        integer, dimension(2)                       :: oldsize,newsize

        call system_clock(ms1)
        !Allocate the spectrum and error ratios
        oldsize(1)=0
        oldsize(2)=0
        newsize(1)=nspec
        newsize(2)=nfiles
        call reallocate2d_r(specrat,nspec,nfiles)
        call reallocate2d_r(errrat,nspec,nfiles)
        !c decide how many spectra to read at mone time
        nspecrq=nspec
!c
!c decide how many blocks are needed
!c
        nblock=nspec/nspecrq
        if(nblock*nspecb.lt.nspec) nblock=nblock+1
        write(6,*) ' Total no. of spectra and time channel boundaries = ',nspec,nchanb
        write(6,*) ' No. of blocks of spectra = ',nblock
        write(6,*) ' No. of spectra per block = ',nspecrq
        write(6,*)
!c
!c set first spectrum number to retrieve
!c
        nspecf=1
!c 
!c step through blocks
!c
        if(allocated(scounts)) deallocate(scounts)
        if(allocated(tcounts)) deallocate(tcounts)
        mssum=0
        do i=1,nblock
!        i=1
!c
!c open a dump file to check data
!c
            open(10,file='get_ratio_new.dat',status='unknown')
!c
!c set last spectrum number to retrieve
!c
            nspecl=nspecf+nspecb-1
            if(nspecl.gt.nspec) nspecl=nspec
            write(6,*)
            write(6,*) 'Spectra being read: ',nspecf,nspecl
            write(6,*)
! Counts arrays need to be allocated BEFORE calling get_data
            call reallocate3d_r(scounts,nchanb,nspec,2)
!c
!c get the data for the first file
!c
            call get_data(runs(1),nperrq(1),mssum)
!c
!c check for spikes
!c
            call find_spike()
!c
!c save these counts
!c
            do is=1,nspecrq
                do j=1,nchan
                    scounts(j,is,2)=tcounts(j,is)
                end do
!                if(is.eq.2050) then
!                    open(20,file='testcounts.txt',status='unknown',iostat=ierr)
!                    do j=1,nchan
!                        write(20,*) j,tcb(j),tcounts(j,is)
!                    end do
!                    close(20)
!                end if
            end do
!c
!c now repeat for each sample file
!c
            do iff=2,nfiles
!c
!c get data for this file
!c
                call get_data(runs(iff),nperrq(iff),mssum)
!c
!c check for spikes
!c
                call find_spike()
!c
!c save these counts
!c
                do is=1,nspecrq
                    do j=1,nchan
                        scounts(j,is,1)=tcounts(j,is)
                    end do
                end do
!c
!c define starting reference integer for counts array
!c - ignore channel 1 which corresponds to channel 0 on DAE
!c
                do is=nspecf,nspecl
!c
!c define initial and final reference integers for counts array
!c
                    isref=is-nspecf+1
                    jf=2
                    jl=jf+nchan-1
!c
!c integrate the respective channels
!c
                    ratsum=0.0
                    eratsum=0.0
                    tsum=0.0
                    do j=jf,jl
                        itref=j-jf+1
                        if(scounts(j,isref,2).gt.0.0) then
                            tsum=tsum+tcw(itref)
                            rat=scounts(j,isref,1)/scounts(j,isref,2)
                            ratsum=ratsum+rat*tcw(itref)
                            eratsum=eratsum+(1.0+rat)*tcw(itref)
                        endif
                    end do
!                    write(6,*) is,tsum
                    if(tsum.gt.0.0) then
!c
!c average ratio and std. dev.
!c
                        specrat(is,iff)=ratsum/tsum
                        errrat(is,iff)=eratsum/tsum/nchan
                    else
                        specrat(is,iff)=0.0
                        errrat(is,iff)=0.0
                    endif
!c
!c update starting reference integer for counts array
!c
                    jf=jl+1
                end do
            end do
!c
!c write out 1 spectrum of data
!c
            nspecwrt=max(int(nspec/2),1)
            nc1=1
            nc2=nchanb
            do j=nc1,nc2
                write(10,202) j,scounts(j,nspecwrt,1),scounts(j,nspecwrt,2)
202             format(1x,i8,2(1x,f12.3))
            end do
            close(10)
            nspecf=nspecl+1
            deallocate(scounts)
            deallocate(tcounts)
        end do
!c
!c normalise the ratios to that of the incident beam monitor and write to .rat 
!c ASCII files
!c
        do iff=2,nfiles
            ilen=len_trim(runs(iff))
            fname=runs(iff)(1:ilen)
            call change_ext(fname,'rat')
            open(10,file=fname,status='unknown')
            write(10,'(a,a)') '#',fname
            write(10,*) '#'
            write(10,*) '#'
            write(10,*) '#'
            rat=specrat(max(1,incid_mon(1)),iff)
            do is=1,nspec
                if(rat.gt.0.0) then
                    specrat(is,iff)=specrat(is,iff)/rat
                    errrat(is,iff)=errrat(is,iff)/rat/rat
                    if(errrat(is,iff).gt.0.0) then
                        err=sqrt(errrat(is,iff))
                    else
                        err=0.0
                    endif                        
                endif
                write(10,101) is,specrat(is,iff),err
101             format(1x,i5,2(1x,e13.6))
            end do
        end do
        call system_clock(ms2)
        write(6,'(a,1x,i6,a)') 'get_ratio> Time spent retrieving data = ',mssum,' ms'
        write(6,'(a,1x,i6,a)') 'get_ratio> Total time taken = ',(ms2-ms1),' ms'
        return        
    end subroutine get_ratio
    
    subroutine check_ratio(runs,irun)

        use run_par
        use bad_detectors
        use groups_routines
        !c
!c CHECK_RATIO.FOR
!c
!!c program to histogram and correlate spectrum ratio data, and use this 
!c information to generate bad detector file, and generate appropriate groups
!c file
!c
!c A.K.Soper, October 1999. 
!c
!c This is a revised version of the routine PASSRATE which was a Genie transform
!c command programme
!c
    character(len=256), intent(in)          :: runs             !run number being compared
    integer, intent(in)                     :: irun            !reference integer for run number
    integer, dimension(:), allocatable      :: nsumgrp
    real, dimension(:), allocatable         :: data        !internal array to access ratios
    real, dimension(:), allocatable         :: err        !internal array to access errors
    real                                    :: yhiin,yloin        !anticipated maximum and minimum ratios
    real, dimension(:), allocatable         :: grpsum,grpsum2,stdup,stddown,stdevo
    real, dimension(:), allocatable         :: errmin,errmax,yhi,ylo
    real                                    :: ymax,emax,ylimit,yfac
    real                                    :: rms,sum1,sum2,sum3,ratlim,diff
    integer                                 :: i,lpt,ngood,ipass,ig,iflag
    
    ymax=0.0
    emax=0.0
    ylimit=1000.0
    if(allocated(data)) deallocate(data)
    if(allocated(err)) deallocate(err)
    allocate(data(nspec),err(nspec))
!c
! determine the maximum ratio and error values
!c
    do lpt=1,nspec
!c
!c ignore very large values
!c
        data(lpt)=specrat(lpt,irun)
        if(errrat(lpt,irun).gt.0.0) then
            err(lpt)=sqrt(errrat(lpt,irun))
        else
            err(lpt)=0.0
        endif
        if(data(lpt).lt.ylimit) then
            if(data(lpt).ge.ymax) ymax=data(lpt)
        endif
        if(err(lpt).gt.emax) emax=err(lpt)
    end do
    ngood=nspec
    write(6,75) ngood
75    format(1x,i4,' good spectra at start')
!c
!c set the initial (broad)limits
    if(allocated(yhi)) deallocate(yhi)
    if(allocated(ylo)) deallocate(ylo)
    if(allocated(stdevo)) deallocate(stdevo)
    if(allocated(stdup)) deallocate(stdup)
    if(allocated(stddown)) deallocate(stddown)
    allocate(yhi(ngroup),ylo(ngroup),stdevo(ngroup),stdup(ngroup),stddown(ngroup))
    if(allocated(grpsum)) deallocate(grpsum)
    if(allocated(grpsum2)) deallocate(grpsum2)
    if(allocated(errmin)) deallocate(errmin)
    if(allocated(errmax)) deallocate(errmax)
    if(allocated(nsumgrp)) deallocate(nsumgrp)
    allocate(grpsum(ngroup),grpsum2(ngroup),errmin(ngroup),errmax(ngroup),nsumgrp(ngroup))
!c
!c we don't expect sample to scatter more than 100 times the calibration or less 
!c than 0.00001 times the calibration
!c
    yloin=0.00001
    yhiin=100
    do ig=1,ngroup
        yhi(ig)=yhiin
        ylo(ig)=yloin
        stdevo(ig)=0.0
        stdup(ig)=0.0
        stddown(ig)=0.0
    end do
    ipass=1
500    continue
!c
!c calculate mean ratio and r.m.s. deviation
!c
    ngood=0
    sum1=0.0
    sum2=0.0
    sum3=0.0
    do ig=1,ngroup
        grpsum(ig)=0.0
        grpsum2(ig)=0.0
        errmin(ig)=0.0
        errmax(ig)=0.0
        nsumgrp(ig)=0
    end do
    write(6,*) data(13),err(13)
    do i=1,nspec
!c
!c only update ibad if it is currently good (=0)
!c
        ig=igrp(i)
        if(ibad(i).eq.0.and.ig.ge.1.and.ig.le.ngroup) then
            if(data(i).lt.ylo(ig)) then
                ibad(i)=-1
            else if(data(i).gt.yhi(ig)) then
                ibad(i)=+1
            else if(ipass.gt.1.and.err(i).gt.stdup(ig)) then
                ibad(i)=-2
            else if(ipass.gt.1.and.err(i).lt.stddown(ig)) then
                ibad(i)=2
            else
!c
!c find maximum and minimum values of non-zero errorbars for all detectors
!c
                if(err(i).gt.errmax(ig)) errmax(ig)=err(i)
                if(err(i).gt.0.0.and.errmin(ig).eq.0.0) errmin(ig)=err(i)
                if(err(i).lt.errmin(ig)) errmin(ig)=err(i)
!c
!c sum ratios
!c
                grpsum(ig)=grpsum(ig)+data(i)
                grpsum2(ig)=grpsum2(ig)+data(i)*data(i)
                nsumgrp(ig)=nsumgrp(ig)+1
                ngood=ngood+1
            endif
        endif
    end do
!c
!c flag to decide whether to make another pass through detectors
!c
    iflag=0
!c
!c if there are any good detectors(!) see if the accepted ones have the
!c expected statistical fluctuations
!c
    do ig=1,ngroup
        if(nsumgrp(ig).gt.0) then
!c
!c mean ratio
!c
            sum1=grpsum(ig)/nsumgrp(ig)
            sum2=grpsum2(ig)/nsumgrp(ig)
            diff=sum2-sum1*sum1
            if(diff.gt.0.0) then
                rms=sqrt(diff)
            else
                rms=0.0
            endif
            grpsum(ig)=sum1
            grpsum2(ig)=rms
!c
!c set the ratio of upper and lower limits on rms deviations to that specified
!c in input file (stdfac)
!c
            if(errmin(ig).gt.0.0) then
                ratlim=errmax(ig)/errmin(ig)
            else
                ratlim=0.0
            endif
!c
!c if this ratio is equal to or smaller than that requested, use these new limits
!c 
            if(ratlim.le.stdfac) then
                stdup(ig)=errmax(ig)
                stddown(ig)=errmin(ig)
            else
!c
!c otherwise ignore bins with progressively higher numbers of entries
!c
                stdup(ig)=errmax(ig)*0.8
                stddown(ig)=errmin(ig)*1.02
            endif
!c 
!c check ratio range and improve if possible
!c
!c if the average value changes by more than 1% attempt a new range of 
!c values based on current average fluctuation about the mean.
!c
!c ignore cases when the scattering is very weak, e.g. sum1 .lt. 20%. In this
!c case the routine simply eliminates any noisy or non-counting detectors
!c
            if(abs(sum1-stdevo(ig)).gt.0.01*stdevo(ig) &
            .and.ipass.lt.10.and.sum1.ge.0.2) then
                iflag=1
!c
!c calculate new limits +- 6 rms deviation
!c
                yfac=1.0+rmsfac*rms/sum1
                ylo(ig)=sum1/yfac    
                yhi(ig)=sum1*yfac
                if(yhi(ig).gt.yhiin) yhi(ig)=yhiin    
                if(ylo(ig).lt.yloin) ylo(ig)=yloin
                stdevo(ig)=sum1
            endif
        endif
    end do
    write(6,200) ipass,ngood,(ig,nsumgrp(ig) &
    ,grpsum(ig),grpsum2(ig),stddown(ig),stdup(ig),ig=1,ngroup)
200    format(1x,'Pass ',i2,'; ',i5,' good spectra. ' &
     ,/(1x,'Grp: ',i3,1x,i4,' Dets, mean ratio: ',e12.5 &
     ,', and mean dev.: ',e12.5 &
     ,/,' Next RMS test values: ',e12.5,1x,e12.5))
!c
!c go for another pass if iflag has been set
!c
    if(iflag.ne.0.and.ipass.lt.10) then
        ipass=ipass+1
        go to 500
    end if
!c
!c write out the new groups, bad spectrum and spike files
!c
    call write_groups(runs)
    return
        
    end subroutine check_ratio

END MODULE spectrum_ratio

PROGRAM purge_det
         
    use inputfilestrings
    use reallocation_routines
    use run_par
    use beam_routines
    use calibration_routines
    use groups_routines
    use monitor_routines
    use get_data_routines
    use spectrum_ratio
    use bad_detectors

    implicit none
!c
!c internal variables
!c
    integer :: i,is,j,k,l,m,n,ierr    !internal indices
    integer :: if,il,jf,jl,kf,kl    !index range limits
    integer :: nfiles,mfiles            !number of runs to compare
    character(len=256), dimension(:), allocatable :: runs        !sample files
    integer, dimension(:), allocatable :: nperrq        !period number to access
    integer :: nspike        !no. of spectra read from spike.dat
    integer :: debug
    integer :: inform,outputunitstype,goodframes,inxs
    integer :: ignore,ms1,ms2        !=1 to ignore existing spec.bad and spike.dat, else 0
    character(len=256) :: fname_nxs        !name of nexus definitions file and detector calibration file
    real :: wavemin,wavemax,muamphrs
    logical :: found

    call system_clock(ms1)
    debug = 1    ! Set to 1 for debugging
    inform = 1    ! Set to 1 for informational messages

!c
!c open input file
!c
    open(15,file='purge_det.dat',status='old',iostat=ierr)

!c Ignore empty lines or lines that begin with spaces. The smallest allowable
!c character is a blank - ASCII decimal 32. Anything less than or equal to this should be ignored.

    nchar=0
    do while(nchar.eq.0.and.ierr.eq.0)
        read(15,'(a)',iostat=ierr) line
        nchar=len_trim(line)
    end do
    if(ierr.ne.0) then
        write(6,'(a,1x,i4)') 'purge_det> ierr = ',ierr
        stop
    end if

!c The first valid line read must contain the strings which signify the spacing between
!c values (spc2), the spacing between the last value and the following comment (spc5), and the 
!c pathseparator character. These values are delineated by apostrophes as the number of 
!c spaces is important for these values. Hence a special reading method is required.
!c
!c Find the first and second occurrence of a '
!c
    index1=index(line(1:nchar),'''')
    index2=index1+index(line(index1+1:nchar),'''')
    spc2=line(index1+1:index2-1)
    lenspc2=1
    index1=index2+index(line(index2+1:nchar),'''')
    index2=index1+index(line(index1+1:nchar),'''')
    spc5=line(index1+1:index2-1)
    lenspc5=index2-index1-1
    index1=index2+index(line(index2+1:nchar),'''')
    index2=index1+index(line(index1+1:nchar),'''')
    pathseparator=line(index1+1:index2-1)
    lenpathseparator=index2-index1-1
    write(6,*) spc2(1:lenspc2),';',spc5(1:lenspc5),';' &
     ,pathseparator(1:lenpathseparator)

!c In the Fortran program it is assumed there may be only a minimum of 1 space between
!c values in a list.

!    spc2=' '
!    lenspc2=1

    call getaline(15)
    inst=line(ncf(1):ncl(1))
    write(6,*) inst(1:index(inst,' '))
    call getaline(15)
    directinp=line(ncf(1):ncl(nwords))
    lendirectinp=ncl(nwords)-ncf(1)+1
    write(6,999) directinp(1:lendirectinp)
999     format(1x,a)
    call getaline(15)
    direct=line(ncf(1):ncl(nwords))
    lendirect=ncl(nwords)-ncf(1)+1
    write(6,999) direct(1:lendirect)
!c
!c detector calibration file - use ' ' if using calibration contained in RAW file
!c
    call getaline(15)
    fnamed=line(ncf(1):ncl(nwords))
    write(6,*) fnamed(1:len_trim(fnamed))

!c
!c Read name of detector groups file. Use ' ' to use default (groups_def.dat)
!c
    call getaline(15)
    fname_grps=line(ncf(1):ncl(nwords))
    ilenfname_grps=ncl(nwords)-ncf(1)+1
    write(6,*) fname_grps(1:ilenfname_grps)
!c
!c read spectrum number of incident beam monitor
!c
    call getaline(15)
    
    allocate(incid_mon(nwords),trans_mon(nwords))
    do i=1,nwords
        read(line(ncf(i):ncl(i)),*) incid_mon(i)
    end do
    write(6,*) incid_mon(1)
    trans_mon(1)=2
!c
!c read first and last channel numbers to check for spikes (0,0 to use all)
!c
    call getaline(15)
    nchfir=0
    nchlas=0
!    read(line(ncf(1):ncl(2)),*) nchfir,nchlas
    write(6,*) nchfir,nchlas
!c
!c read flag to decided whether to check for spikes in the data and allowed number
!c of std's to determine whether spike has occurred.
!c
    call getaline(15)
    read(line(ncf(1):ncl(1)),*) spikedev
    write(6,*) spikedev
    if(spikedev.gt.0.0) then
        flagsp=1
    else
        flagsp=0
    endif
    spikedev=spikedev*spikedev
!Get the name of the NeXus information file if is present
    call getaline(15)
    inxs=index(line(ncf(1):ncl(nwords)),'nexus_txt')
    if(inxs.gt.0) then
        fname_nxs=line(ncf(1):ncl(nwords))
        call getaline(15)
    end if
!c
!c read number of rms deviations to allow in CHECK_RATIO, and maximum allowed
!c variation of standard deviation of ratio ABOVE the average value
!c
    read(line(ncf(1):ncl(2)),*) rmsfac,stdfac
    write(6,*) rmsfac,stdfac

!c Read ignore: this determines whether spec.bad and spike.dat will be ignored or not. 
!c 1 means spec.bad and spike.dat will be ignored. Anything else and it will be read.

    call getaline(15)
    read(line(ncf(1):ncl(1)),*) ignore
    write(6,*) ignore

    nfiles=0
    mfiles=0
    call getaline(15)
    do while(nwords.gt.0)
        if(nwords.gt.1) then
            !First check that this filename has not already been included
            !It is assumed the filename occupies words 1 - nwords-1 with the period number in the last word
            i=0
            found=.false.
            do while(i.lt.nfiles.and..not.found)
                i=i+1
                found=line(ncf(1):ncl(nwords-1)).eq.runs(i)(1:len_trim(runs(i)))
            end do
            if(.not.found) then
                nfiles=nfiles+1
                do while(nfiles.gt.mfiles)
                    !allocate more space
                    mfiles=mfiles+10
                    call reallocate1d_c(runs,len(runs),mfiles)
                    call reallocate1d_i(nperrq,mfiles)
                end do
                runs(nfiles)=line(ncf(1):ncl(nwords-1))
                read(line(ncf(nwords):ncl(nwords)),*) nperrq(nfiles)
!        write(6,'(1x,i4,1x,a,1x,i4)') nfiles,runs(nfiles)(1:len_trim(runs(nfiles))),nperrq(nfiles)
            end if
        endif
        call getaline(15)
!        write(6,*) nwords,nfiles,mfiles
    end do
    if(nfiles.eq.0) then
        write(6,'(a)') 'purget_det> No data files found!'
        stop
    end if
    do i=1,nfiles
        write(6,*) runs(i)(1:len_trim(runs(i)))
    end do
    
!c For nexus files we need to read in the list of paths needed to find particular data.
!c The filename will be called 'instrument'_nexus.txt

    ext=get_ext(runs(1))
    write(6,'(1x,a)') ext(1:len_trim(ext))
    if(ext.eq.'nxs'.or.ext.eq.'NXS') then
        if(inxs.gt.0) then
            call get_pathnames_nxs(fname_nxs)
        else
            write(6,*) 'purge_det> Valid NeXus pathnames definition file not specified'
            stop
        end if
    endif
!c
!c initialise run parameters
!c
    call get_run_par(runs,nfiles)
    write(6,*) 'purge_det> Done get_run_par'
!    call get_title(runs(1))
!    call get_goodframes(runs(1),goodframes)
!    call get_muamphrs(runs(1),muamphrs)
    
!c
!c Get the time channel boundaries for this run - will assume they are the same for all runs
!c
    call get_tcbs(runs(1))

! Initialise the time channel boundaries and widths

    nchanb=nchanbs
    nchan=nchans
    if(allocated(tcb)) deallocate(tcb)
    if(allocated(tcw)) deallocate(tcw)
    allocate(tcb(nchanb),tcw(nchan))
    do i=1,nchanb
        tcb(i)=tcbs(i)
    end do
    do i=1,nchan
        tcw(i)=tcws(i)
    end do

!c initialise calibration parameters
!c
!c    WRITE(6,*) 
!!c    WRITE(6,*) DIRECT
!c    WRITE(6,*)
!c preset the incident flight path to zero so that the value will be read from
!c header
!Notional incident flight to allow get_calibration to proceed - not used in purge_det
    incidentfp=100.0
    call get_calibration(runs,nfiles)
!c    WRITE(6,*) 
!!c    WRITE(6,*) DIRECT
!c    WRITE(6,*)
!c
!c read existing groups, bad spectra and spike files if they exist
!c Provide some dummy wavelengths - these will not be used
    wavemin=0.05
    wavemax=5.0
    outputunitstype=1
    call get_groups(ignore,wavemin,wavemax,outputunitstype)
    write(6,*) 'Done groups'
!c
!c form integrated ratio for each spectrum
!c
    call get_ratio(runs,nfiles,nperrq)
    write(6,*) 'Done get ratio'
!c
!c form histogram of ratios and error bars and use this information to
!c determine good and bad detectors. 
!c Write groups file accordingly.
!c
    do i=2,nfiles
        call check_ratio(runs(i),i)
    end do
    call system_clock(ms2)
    write(6,'(a,1x,i10,a)') 'purge_det> Total run time = ',(ms2-ms1), ' ms'
    stop
    
    END program purge_det
