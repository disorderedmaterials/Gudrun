!     
! File:   get_ratio.f90
! Author: aks45
!
! Created on 30 April 2012, 13:39
!

MODULE get_data_routines

    implicit none
      
CONTAINS

    subroutine get_pathnames_nxs(fname)
        
        USE inputfilestrings
        USE run_par
        USE local_data
    
        
        character(len=256), intent(inout) :: fname
        integer i,ierr
        
        open(16,file=fname,status='old',iostat=ierr)
        if(ierr.ne.0) then
            fname=convertfilename(fname)
            open(16,file=fname,status='old',iostat=ierr)
        endif
        if(ierr.ne.0) then 
            write(6,'(a,a,a)') 'File ',fname(1:len_trim(fname)),' does not exist!'
            stop
        endif
        !Initialise the path names
        titlepath=' '
        usernamepath=' '
        starttimepath=' '
        tcbpath=' '
        periodspresent=.true. !Default for ISIS data. Others may be not, so make sure specify 'no' or anything other than 'yes'!
        goodframespath=' '
        muamphrspath=' '
!c Read the list of paths
        nspectrumpaths=0
        mspectrumpaths=0
        periodspresent=.false.
        call getaline(16)
        do while(nwords.gt.0)
            if(line(ncf(1):ncl(1)).eq.'title') titlepath=line(ncf(3):ncl(3))
            if(line(ncf(1):ncl(1)).eq.'username') usernamepath=line(ncf(3):ncl(3))
            if(line(ncf(1):ncl(1)).eq.'starttime') starttimepath=line(ncf(3):ncl(3))
            if(line(ncf(1):ncl(1)).eq.'timechannels') tcbpath=line(ncf(3):ncl(3))
            if(line(ncf(1):ncl(1)).eq.'periodspresent') periodspresent=line(ncf(3):ncl(3)).eq.'yes'
            if(line(ncf(1):ncl(1)).eq.'goodframes') goodframespath=line(ncf(3):ncl(3))
            if(line(ncf(1):ncl(1)).eq.'muamphrs') muamphrspath=line(ncf(3):ncl(3))
            if(line(ncf(1):ncl(1)).eq.'spectrum') then
                nspectrumpaths=nspectrumpaths+1
                do while (nspectrumpaths.gt.mspectrumpaths)
                    mspectrumpaths=mspectrumpaths+10
                    call reallocate1d_c(spectrumcountspath,len(spectrumcountspath),mspectrumpaths)
                    call reallocate1d_c(spectrummodulevalues,len(spectrummodulevalues),mspectrumpaths)
                end do
                spectrumcountspath(nspectrumpaths)=line(ncf(3):ncl(3))
                if(nwords.gt.3) then
                    spectrummodulevalues(nspectrumpaths)=line(ncf(4):ncl(nwords))
                else
!c If no module values are read in, then it is assumed all the detectors in this data group belong to
!c the same module.
                    spectrummodulevalues(nspectrumpaths)=' '
                endif
            endif
            call getaline(16)
        end do
        close (16)
        write(6,'(1x,a,a)') 'title : ',titlepath(1:50)
        write(6,'(1x,a,a)') 'timechannels : ',tcbpath(1:50)
        write(6,'(1x,a,l)') 'periodspresent : ',periodspresent
        i=0
        do while(i.lt.nspectrumpaths)
            i=i+1
            write(6,'(a,a)') 'spectrum : ',spectrumcountspath(i)(1:50)
        end do
        
    end subroutine get_pathnames_nxs
    

    subroutine get_data(runno,nperrq,mssum)
        
        use reallocation_routines
        use inputfilestrings
        USE run_par
        USE local_data
    
!***********************************************************************************
!*
!*        GET_DATA.FOR
!*
!*        Copyright A K Soper, January 2003
!*         derived from the ISIS GET routines by Freddie Akeroyd
!*
!*        This program reads some data from a standard ISIS raw data file 
!*        according to predefined criteria.
!*
!***********************************************************************************
        character(len=256), intent(in)          :: runno                        !filename of dataset
        integer, intent(inout)                  :: nperrq                !period number to extract
        integer, intent(inout), optional        :: mssum
!c
!c internal variables
!c
        character(len=256)                      :: fname                !name of file to read
        integer                                 :: i,j,iref,jf,jl,ndum1,ndum2,nuse        !dummy variables
        integer, dimension(:), allocatable      :: buffer
        integer                                 :: errcode                !errcode on exit from reading data
        integer                                 :: ms1,ms2  !Records the number of clock pulses
!c
!c define total no. of spectra being read
!c
        nspecrq=nspecl-nspecf+1
!c
!c define total no. of channels being read
!c
        nchanrq=nspecrq*nchanb
        countsdims(1)=nchanb
        countsdims(2)=nspecrq
        countsdims(3)=1
!        write(6,*) 'get_data> ',nspecf,nspecl,nspecrq,nchanrq,size(tcounts,1),size(tcounts,2)
! Counts arrays need to be allocated BEFORE calling get_data
!        if(allocated(tcounts)) deallocate(tcounts)
        call reallocate2d_r(tcounts,nchanb,nspecrq)
!        write(6,*) 'get_data> ',nspecf,nspecl,nspecrq,nchanrq,countsdims(1),countsdims(2),countsdims(3)
        call reallocate3d_i(counts,nchanb,nspecrq,1)
!        write(6,*) 'get_data> ',nspecf,nspecl,nspecrq,nchanrq,countsdims(1),countsdims(2),countsdims(3)
        call system_clock(ms1)
!c
!c ensure period requested is sensible - if not use period 1
!c
        if(nperrq.le.0.or.nperrq.gt.nper) nperrq=1
!c
!c determine number of characters in filename
!c
        nchar=len_trim(runno)
!c
!c generate filename to read
!c
        fname=direct(1:lendirect)//runno(1:nchar)
        write(6,*) 'get_data> Filename: ',fname(1:len_trim(fname))

        if(inst(1:len_trim(inst)).eq.'D4C') then

            call get_data_D4C(fname)
            do i=1,nspecrq
                do j=1,nchan
                    tcounts(j,i)=real(counts(j,i,1))
                end do
                tcounts(nchanb,i)=0
                if(i.eq.2886.or.i.eq.2890) write(6,*) 'get_data> ',i,tcounts(nchan,i)
            end do

        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

            call get_data_nxs(fname,nperrq)
            do i=1,nspecrq
                do j=1,nchan
                    tcounts(j,i)=real(counts(j,i,1))
                end do
                tcounts(nchanb,i)=0
            end do
        else

!c
!c open the raw file
!c
            call open_data_file(fname,ndum1,ndum2,nuse,errcode)
            if(errcode.ne.0) then
                write(6,'(a,a)') 'Error code',errcode,' when opening file',fname
                write(6,'(4(1x,i4))') nchan,ndum1,ndum2,errcode
            endif
            mcount=countsdims(1)*countsdims(2)*countsdims(3)
 !                if(allocated(buffer)) deallocate(buffer)
!                allocate(buffer(mcount))
!c
!c define first spectrum number to access
!c
            nf=(nperrq-1)*nspecb+nspecf
            call getdat(fname,nf,nspecrq,counts,mcount,errcode)
            if(errcode.ne.0) then
                write(6,201) errcode
201             format(1x,'Error code ',i3,' in get_data')
                stop
            endif
            do i=1,nspecrq
                do j=1,nchan
!Don't forget to ignore the infamous channel 0!
                    tcounts(j,i)=real(counts(j+1,i,1))
                end do
                tcounts(nchanb,i)=0
            end do

        endif
        call system_clock(ms2)
        if(present(mssum)) mssum=ms2-ms1
        if(allocated(counts)) deallocate(counts)
        return

    end subroutine get_data
    
    subroutine get_goodframes(runno,goodframes)
!***********************************************************************************
!*
!*        GET_DATA.FOR
!*
!*        Copyright A K Soper, January 2003
!*         derived from the ISIS GET routines by Freddie Akeroyd
!*
!*        This program reads some data from a standard ISIS raw data file 
!*        according to predefined criteria.
!*
!***********************************************************************************
        USE run_par
        USE local_data
        USE datasetlist
        character(len=256), intent(in)          :: runno                        !filename of dataset
        integer, intent(out)                    :: goodframes                 !no. of good frames for this run
!c
!c internal variables
!c
        character(len=256)                      :: fname                !name of file to read
        integer                                 :: i,j,iref,jf,jl,ndum1,ndum2,nuse        !dummy variables
        integer, dimension(:), allocatable      :: buffer
        integer                                 :: errcode                !errcode on exit from reading data
        INTEGER                                     :: idataset0d
        INTEGER, DIMENSION(:), ALLOCATABLE, TARGET     :: idataset1d
        INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET   :: idataset2d
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: idataset3d
        integer                                         :: msdiff,mssum,ms1,ms2
!c
!c determine number of characters in filename
!c
        nchar=len_trim(runno)
!c
!c generate filename to read
!c
        fname=direct(1:lendirect)//runno(1:nchar)

        if(inst(1:len_trim(inst)).eq.'D4C') then

            goodframes=1

        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

            !retrieve the good frame count
            call retrieve_data_i(fname,goodframespath &
            ,idataset0d,idataset1d,idataset2d,idataset3d)
            goodframes=idataset0d
            
        else

!c
!c open the raw file
!c
            call open_data_file(fname,ndum1,ndum2,nuse,errcode)
            if(errcode.ne.0) then
                write(6,'(a,a)') 'Error code',errcode,' when opening file',fname
                write(6,'(4(1x,i4))') nchan,ndum1,ndum2,errcode
            endif
            allocate(buffer(32))
            call getpari(fname,'RPB',buffer,32,nuse,errcode)
            goodframes=buffer(10)
            call close_data_file()
            deallocate(buffer)
        
        endif
        write(6,'(a,1x,i10,a,1x,a)') ' get_goodframes> Retrieved ' &
        ,goodframes,'  goodframes from:-',fname(1:len_trim(fname))
        return
    end subroutine get_goodframes
    
    subroutine get_muamphrs(runno,muamphrs)
!***********************************************************************************
!*
!*        GET_DATA.FOR
!*
!*        Copyright A K Soper, January 2003
!*         derived from the ISIS GET routines by Freddie Akeroyd
!*
!*        This program reads some data from a standard ISIS raw data file 
!*        according to predefined criteria.
!*
!***********************************************************************************
        USE run_par
        USE local_data
        USE datasetlist
        character(len=256), intent(in)                  :: runno                        !filename of dataset
        real, intent(out)                               :: muamphrs                 !no. of good frames for this run
!c
!c internal variables
!c
        character(len=256)                              :: fname                !name of file to read
        integer                                         :: i,j,iref,jf,jl,ndum1,ndum2,nuse        !dummy variables
        real, dimension(:), allocatable                 :: buffer
        integer                                         :: errcode                !errcode on exit from reading data
        real(kind=4), DIMENSION(:), ALLOCATABLE, TARGET         :: rdataset1d
        real(kind=8), DIMENSION(:), ALLOCATABLE, TARGET         :: rrdataset1d
!c
!c determine number of characters in filename
!c
        nchar=len_trim(runno)
!c
!c generate filename to read
!c
        fname=direct(1:lendirect)//runno(1:nchar)

        if(inst(1:len_trim(inst)).eq.'D4C') then

            muamphrs=1

        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

            !retrieve the good frame count
            call retrieve_data_r(fname,muamphrspath,rdataset1d,rrdataset1d)
            if(inst(1:len_trim(inst)).eq.'NOVA') then
                    muamphrs=real(rrdataset1d(1))
            else
                    muamphrs=rdataset1d(1)
            end if

            
        else

!c
!c open the raw file
!c
            call open_data_file(fname,ndum1,ndum2,nuse,errcode)
            if(errcode.ne.0) then
                write(6,'(a,a)') 'Error code',errcode,' when opening file',fname
                write(6,'(4(1x,i4))') nchan,ndum1,ndum2,errcode
            endif
            allocate(buffer(32))
            call getparr(fname,'RRPB',buffer,32,nuse,errcode)
            muamphrs=buffer(8)
            call close_data_file()
            deallocate(buffer)
        
        endif
        write(6,'(a,1x,e12.5,a,a)') ' get_muamphrs> Retrieved ' &
        ,muamphrs,'  micro-amp-hrs from:-',fname(1:len_trim(fname))
        return
    end subroutine get_muamphrs
    
    subroutine get_title(runno)
!***********************************************************************************
!*
!*        GET_DATA.FOR
!*
!*        Copyright A K Soper, January 2003
!*         derived from the ISIS GET routines by Freddie Akeroyd
!*
!*        This program reads some data from a standard ISIS raw data file 
!*        according to predefined criteria.
!*
!***********************************************************************************
        USE run_par
        USE local_data
        USE datasetlist
        character(len=256), intent(in)          :: runno                        !filename of dataset
!c
!c internal variables
!c
        character(len=256)                      :: fname,sometext                !name of file to read
        integer                                 :: i,j,iref,jf,jl,ndum1,ndum2,nuse,ifirst,ilast       !dummy variables
        integer                                 :: errcode                !errcode on exit from reading data
        integer :: foundoldformat,foundnewformat,istart
        integer :: iyear,imonth,iday,ihour,imin,isec
        INTEGER, TARGET                                :: idataset0d
        INTEGER, DIMENSION(:), ALLOCATABLE, TARGET     :: idataset1d
        INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET   :: idataset2d
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: idataset3d
        real(kind=4), DIMENSION(:), ALLOCATABLE, TARGET     :: rdataset1d
        real(kind=8), DIMENSION(:), ALLOCATABLE, TARGET     :: rrdataset1d
!c
!c determine number of characters in filename
!c
        nchar=len_trim(runno)
!c
!c generate filename to read
!c
        fname=direct(1:lendirect)//runno(1:nchar)

        if(inst(1:len_trim(inst)).eq.'D4C') then

            open(10,file=fname,status='old',iostat=errcode)
            if(errcode.ne.0) write(6,*) 'Error code',errcode,' when opening file',fname(1:lendirect+nchar)

            foundoldformat=0
            foundnewformat=0
            do while((foundoldformat.eq.0.and.foundnewformat.eq.0).and.errcode.eq.0)
                read(10,'(a)',iostat=errcode) sometext
                foundoldformat=index(sometext,'SAMPLE             :')
                foundnewformat=index(sometext,'Subtitle           :')
            end do

!c Search for the ':' - title is the text that follows this.

            istart=index(sometext,':')+1
            info=sometext(istart:istart+79)
            write(6,*) 'get_run_title> ',info(1:80)
! Rewind to get the user name
            rewind(10)
            foundoldformat=0
            foundnewformat=0
            do while((foundoldformat.eq.0.and.foundnewformat.eq.0).and.errcode.eq.0)
                read(10,'(a)',iostat=errcode) sometext
                foundoldformat=index(sometext,'Experimentalist    :')
                foundnewformat=index(sometext,'User Name          :')
            end do
            istart=index(sometext,':')+1
            uname=sometext(istart:istart+19)
! Rewind to get the start time
            rewind(10)

!c Search for the word 'Start time'

            sometext=' '
            do while(index(sometext,'Start time of numor:').lt.1.and.errcode.eq.0)
                read(10,'(a)',iostat=errcode) sometext
            end do
            istart=index(sometext,':')+1
            sttime=sometext(istart:istart+19)
!            write(6,*) 'get_run_title> ',info(1:len_trim(info)

        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

            !retrieve the run title
            if(len_trim(titlepath).gt.0) then
                info=retrieve_data_c(fname,titlepath)
            else
                info='No title provided'
            end if
            if(len_trim(usernamepath).gt.0) then
                uname=retrieve_data_c(fname,usernamepath)
            else
                uname='No username provided'
            end if
            if(len_trim(starttimepath).gt.0) then
                i=retrieve_dims(fname,starttimepath)
                if(class.eq.H5T_STRING_F) then
                    sttime=retrieve_data_c(fname,starttimepath)
                else if ((class.eq.H5T_FLOAT_F.or.class.eq.H5T_INTEGER_F) &
                    .and.ndims.eq.1.and.dims(1).gt.5) then
                    if(class.eq.H5T_FLOAT_F) then
                        call retrieve_data_r(fname,starttimepath,rdataset1d,rrdataset1d)
                        if(precision.eq.32) then
                            iyear=int(rdataset1d(1))
                            imonth=int(rdataset1d(2))
                            iday=int(rdataset1d(3))
                            ihour=int(rdataset1d(4))
                            imin=int(rdataset1d(5))
                            isec=int(rdataset1d(6))
                        else if(precision.eq.64.and.dims(1).gt.5) then
                            iyear=int(rrdataset1d(1))
                            imonth=int(rrdataset1d(2))
                            iday=int(rrdataset1d(3))
                            ihour=int(rrdataset1d(4))
                            imin=int(rrdataset1d(5))
                            isec=int(rrdataset1d(6))
                        end if
                    else if(class.eq.H5T_INTEGER_F) then
                        call retrieve_data_I(fname,starttimepath, &
                        idataset0d,idataset1d,idataset2d,idataset3d)
                        iyear=int(idataset1d(1))
                        imonth=int(idataset1d(2))
                        iday=int(idataset1d(3))
                        ihour=int(idataset1d(4))
                        imin=int(idataset1d(5))
                        isec=int(idataset1d(6))
                    end if
                else
                    iyear=0
                    imonth=0
                    iday=0
                    ihour=0
                    imin=0
                    isec=0
                end if
                write(sttime,'(i4.4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2)') &
                iyear,imonth,iday,ihour,imin,isec
            else
                sttime='No start time provided'
            end if
        else

!c
!c open the raw file
!c
            call open_data_file(fname,ndum1,ndum2,nuse,errcode)
            if(errcode.ne.0) then
                write(6,'(a,a)') 'Error code',errcode,' when opening file',fname
                write(6,'(4(1x,i4))') nchan,ndum1,ndum2,errcode
            endif
  	    call getparc(fname,'TITL',info,1,nuse,errcode)
	    call getparc(fname,'HDR',headinfo,1,nuse,errcode)
            call close_data_file()
!c
!c Now extract the username and start time/date from the
!c headinfo variable
!c
            uname=headinfo(9:28)
            sttime=headinfo(53:72)
        
        endif
        write(6,'(a,1x,a,1x,a)') ' get_title> Title retrieved =',info(1:len_trim(info)),'from:-',fname(1:len_trim(fname))
        return
    end subroutine get_title
    
    subroutine get_periods(runno,nper_return)
        
        !Determines the number of periods in a dataset
        
        USE run_par
        USE local_data
        USE datasetlist

        character(len=256), intent(in)          :: runno                        !filename of dataset
        integer                                 :: nper_return
        !c
!c internal variables
!c
        character(len=256)                      :: fname,sometext                !name of file to read
        integer                                 :: i,j,iref,jf,jl,ndum1,ndum2,nuse,ifirst,ilast,icountpaths       !dummy variables
        integer                                 :: errcode                !errcode on exit from reading data
        integer, dimension(1)                   :: idumaks         !dummy variables
        logical                                 :: found
!c
!c generate filename to read
!c
        fname=direct(1:lendirect)//runno

        if(inst(1:len_trim(inst)).eq.'D4C') then
            
            nper_return=1
            
        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then
            icountpaths=0
            nper_return=0
            do while (icountpaths.lt.nspectrumpaths)
                icountpaths=icountpaths+1
                found=.false.
                !Find this spectrum path in the dataset list
                found=retrieve_dims(fname,spectrumcountspath(icountpaths)).gt.0
                if(found) then
                    if(periodspresent.and.ndims.gt.2) then

!c The number of spectra will the product of dimensions for ranks 2 or more. If only rank 1 is
!c present then this path contributes just one spectrum. If periods are present, the largest rank
!c is the period number and so should not be used to assess the number of spectra.

                        nper_return=max(nper_return,dims(ndims))
                    else
                        nper_return=max(nper_return,1)
                    endif
                else
                    write(6,'(1x,a,a,a)') 'Dataset: '&
                    ,spectrumcountspath(icountpaths)(1:len_trim(spectrumcountspath(icountpaths)))&
                    ,' not found'
                    nper_return=max(nper_return,1)
                endif
            end do
        else
            call open_data_file(fname,nchan,ndet,nuse,errcode)
            if(errcode.ne.0) then
                write(6,*) 'Error code',errcode &
                ,' when opening file',fname(1:lendirect+nchar)
                stop
            end if
            call getpari(fname,'NPER',idumaks,1,nuse,errcode)
            nper_return=idumaks(1)
        end if
        
    end subroutine get_periods
    
    subroutine get_data_D4C(fname)
!***********************************************************************************
!*
!*        getdat_D4C.for
!*
!*        Copyright A K Soper, March 2005
!*         Written to read the raw counts from a D4C data file
!*
!***********************************************************************************
        USE run_par
        use reallocation_routines
        use inputfilestrings
        USE local_data
        character(len=256), intent(in)      :: fname         !Filenames to be inspected
        character(len=256)                  :: sometext
        real                                :: angarm,tempmon,temptime
        integer, dimension(:), allocatable  :: tempcounts
        integer                             :: i,j,iflag,imod,errcode,iflag1,iflag2
        integer                             :: ifirst,ilast,iref,iref1,iref2,is,isref,icref
        integer                             :: foundformat,nskip
        logical                             :: test

!c open the raw file

        open(10,file=fname,status='old',iostat=errcode)
        if(errcode.ne.0) then
            write(6,100) errcode,fname(1:50)
100         format(/'get_data_D4C> Error code ',i3,' when opening file: ',a)
            stop
        endif

!c Search for either the word 'MonitorCnts'

        sometext=' '
        foundformat=0
        do while(foundformat.eq.0.and.errcode.eq.0)
            read(10,'(a)',iostat=errcode) sometext
            foundformat=index(sometext,'MonitorCnts')
        end do
        if(errcode.ne.0) then
            write(6,101) errcode,fname(1:len_trim(fname))
101         format(/'get_data_D4C> Error code ',i3,' when reading file: ',a)
            stop
        endif
        nskip=7
!c Skip nskip lines and read the monitor count
        do j=1,nskip
            read(10,*)
        end do
        read(10,*) tempmon,temptime
!c Skip 1 more line and read the angular position of the detector arm
        read(10,*)
        read(10,*) angarm,angarm
!c Read 9 blocks of 64
        nproc=2+9*64
        iref1=1
        call reallocate1d_i(tempcounts,nproc)
        tempcounts(iref1)=tempmon
        iref1=iref1+1
! Store the time in the second channel of this module. Because this is stored as an
! integer it is multiplied by 1000 to give 1/1000ths of a second
        tempcounts(iref1)=nint(temptime*1000.0)
        iref1=iref1+1
        do i=1,9
! Now search for the sequence 'IIIIIIIIIIIIIIII' to indicate the beginning of a block of data
            foundformat=0
            do while (foundformat.eq.0.and.errcode.eq.0)
                read(10,'(a)',iostat=errcode) sometext
                foundformat=index(sometext,'IIIIIIIIIIIIIIII')
            end do
            if(errcode.ne.0) then
                write(6,101) errcode,fname(1:len_trim(fname))
                stop
            endif
! It is assumed there are 64 detectors in a block, so we don't read this number
            read(10,*)
            iref2=iref1+63
            read(10,*) (tempcounts(iref),iref=iref1,iref2)
            iref1=iref2+1
        end do

        !c Step through the existing modules and find the module which these data correspond to
        imod=1
        iflag=0
        do while(imod.lt.nmod.and.iflag.eq.0)
            imod=imod+1
            if(abs(moduleangle(imod)-angarm).lt.0.0001) then
                   iflag=imod
            endif
        end do
        if(iflag.gt.0) then
!c The monitor is stored in the module 2 less than this one. The time is stored in the module 1 less than this one.
! Set up the range of module numbers corresponding to these data read in.
            iflag1=iflag-2
            iflag2=iflag1+9
!c           write(6,*) iflag1,iflag2,fname(1:50)
!c Now determine the first and last spectrum numbers corresponding to this module
            i=0
            ifirst=0
            ilast=0
            test=.true.
            do while(i.lt.nspec.and.test)
                i=i+1
                if(imodule(detno(i)).eq.iflag1.and.ifirst.eq.0) ifirst=i
                if(imodule(detno(i)).eq.iflag2) ilast=i
                test=imodule(detno(i)).le.iflag2
            end do
            write(6,*) 'get_data_D4C> ',nspecf,nspecl,ifirst,ilast
!c Set up the counts array depending on the values of ifirst and ilast
            is=nspecf-1
            do while(is.lt.nspecl)
                is=is+1
                icref=1
!c For D4C data only 2 channels per spectrum. The last contains nothing as for other histogram formats
!c The first contains data if the spectrum number requested is within the range read in
                if(is.ge.ifirst.and.is.le.ilast) then
                   isref=is-ifirst+1
                   counts(icref,is,1)=tempcounts(isref)
                else
                   counts(icref,is,1)=0
                endif
                icref=icref+1
                counts(icref,is,1)=0
            end do
        endif
        return
    end subroutine get_data_D4C
    
    subroutine get_data_nxs(filename,nperrq)
        !Gets the specified data starting at spectrum nf and going to spectrum nl. Output in the
        !counts array. Only 3-D datasets can be processed at the present time
        
        USE run_par
        USE local_data
        USE datasetlist
        character(len=256), intent(in)              :: filename
        integer, intent(in)                         :: nperrq
!        integer                                     :: retrieve_dims  !Used to find the data type and dimensions of a specified dataset
        integer                                     :: i,j,idset,icountsref,ibufferref
        logical                                     :: notfound,found
        integer                                     :: nspecfirst,nspeclast,nspecret,nspecrettotal
        integer                                     :: dimen3,first3d,last3d,range3d  !Range for third dimension
        integer                                     :: dimen2     !Range for second dimension
        real(kind=4), DIMENSION(:), ALLOCATABLE, TARGET     :: rdataset1d
        real(kind=8), DIMENSION(:), ALLOCATABLE, TARGET     :: rrdataset1d

        i=1
        found=.false.
        nspecfirst=nspecf
        nspecrettotal=0
        icountsref=0
        !It is assumed the counts array to store the data has been allocated externally
        do while (nspecrettotal.lt.nspecrq.and.i.le.nspectrumpaths)
            found=nspecfirst.le.lastspectrumthispath(i)
            if (found) then
                !Retrieve the identity and dimensions of this spectrum path
                idset=retrieve_dims(filename,spectrumcountspath(i))
                !Get the last spectrum that can be retrieved from this dataset
                nspeclast=min(nspecl,lastspectrumthispath(i))
                !get the coordinates of the first spectrum in the current dataset
                do while (nspecrettotal.lt.nspecrq.and.nspecfirst.le.nspeclast)
                    if(ndims.gt.1) then
                        dimen3=int(nspecfirst-firstspectrumthispath(i)-1)/dims(2)+1
                        if(periodspresent) dimen3=nperrq
                        dimen2=nspecfirst-firstspectrumthispath(i)-dims(2)*(dimen3-1)
                        nspecret=min(dimen2+nspecrq-nspecrettotal-1,dims(2))-dimen2+1
                    else 
                        nspecret=1
                    end if
                    nspecrettotal=nspecrettotal+nspecret
!                    write(6,'(8(1x,i5))') nspecf,nspecl,nspecrq,dimen2,nspecret &
!                    ,nspecrettotal,dimen3
                    !Retrieve the data!
                    if(ndims.eq.3) then
                        countsoffsets(1)=0
                        countsoffsets(2)=nspecfirst-nspecf
                        countsoffsets(3)=0
                        call retrieve_dataslab_i(filename,spectrumcountspath(i) &
                        ,dimen2,nspecret,dimen3,counts,countsdims,countsoffsets)
                    else if(ndims.eq.1.and.nspecret.eq.1) then !If the data is stored in 1D array
                        call retrieve_data_r(filename,spectrumcountspath(i) &
                        ,rdataset1d,rrdataset1d)
                        if(precision.eq.32) then
                            do j=1,nchan
                                counts(j,nspecrettotal,1)=nint(rdataset1d(j))
                            end do
                        else
                            do j=1,nchan
                                counts(j,nspecrettotal,1)=idnint(rrdataset1d(j))
                            end do
                        end if
                        counts(nchanb,nspecrettotal,1)=0
                    end if
                    nspecfirst=nspecfirst+nspecret
                end do
            else
                i=i+1
            end if
        end do
        
    end subroutine get_data_nxs

END MODULE get_data_routines
