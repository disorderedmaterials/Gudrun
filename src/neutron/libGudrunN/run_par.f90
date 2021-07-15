
!     
! File:   run_par.f90
! Author: aks45
!
! Created on 27 April 2012, 10:11
!

MODULE run_par
!
! global run variables common to all routines
!
    implicit none
      
    character(len=80)  :: info        !full title of run
    character(len=80)  :: headinfo        !extra header info
    character(len=80)  :: uname            !username
    character(len=80)  :: sttime            !start time and date
    character(len=80)  :: inst        !instrument name
    character(len=10)  :: ext            !raw file extension
    character(len=3)   :: pathseparator    !string used to described the character(s) between folder names
! folder and file names
    character(len=256) :: directinp        !directory where Gudrun is being run
    character(len=256) :: direct        !name of directory containing raw files
    character(len=256) :: gudrunstartupfolder   !Defines the folder where GudrunGUI was started
    character(len=256) :: gudrunstartupfilefolder    !Defines the folder containing the startup file and instrument specific files
    character(len=256) :: neutronparfilename    !name of file containing neutron scattering cross sections and scattering lengths

    integer :: nspec,mspec            !no. of spectra in 1 period (not including spectrum 0)        
    integer :: nspecb        !no. of spectra in 1 period including spectrum 0 if present
    integer :: nproc            !no. of spectra to process at one time
    integer :: nchan,nchanv,nchans    !no. of time channels (TC's)
    integer :: nchanb,nchanbv,nchanbs !no. of TC boundaries ( = nchan+1)
    integer :: nspch,nspchv,nspchs     !product of no. spec. times no. of tcb's
    integer :: nper            !no. of periods in this file
    integer :: ndet,mdet            !no. of detectors in this file
    integer :: nmod,mmod            !number of detector modules (DAE inputs)
    integer, dimension(:), allocatable :: specno        !spectrum numbers of detectors
    integer, dimension(:), allocatable :: detno        !detector numbers for spectra
    integer, dimension(:), allocatable :: udet        !labels for detectors
    integer, dimension(:), allocatable :: crateno        !crate numbers of modules
    integer, dimension(:), allocatable :: slotno        !slot numbers of modules
    integer, dimension(:), allocatable :: ndetpermodule    !no. of detectors in a module
    integer, dimension(:), allocatable :: crate        !crate number for each detector
    integer, dimension(:), allocatable :: slot        !slot number of each detector
    integer, dimension(:), allocatable :: imodule        !module number of each detector
    integer lenpathseparator    !number of characters in path separator
    integer lendirect,lendirectinp    !length in characters of direct and directinp
    integer ilengudrunstartupfolder   !length of gudrunstartupfoldername
    integer ilengudrunstartupfilefolder   !length of gudrunstartupfilefoldername

    real, dimension(:), allocatable     :: tcb,tcbv,tcbs !time channel boundaries (musec)
    real, dimension(:), allocatable     :: tcw,tcwv,tcws !time channel widths (musec)
    real, dimension(:), allocatable     :: moduleangle    !module scattering angles
! D4C efficiency values. effd4c will be assigned to each spectrum, while effd4cinput are read in
! if the file 'effd4c.eff' is available. Otherwise the efficiencies are assumed to be unity.
    real, dimension(:), allocatable     :: effd4c        
    real, dimension(:,:), allocatable   :: effd4cinput
    real                                :: sampledependentbackgroundfactor
    real                            :: incidentfp        !Incident flight path from input file
    real                            :: incidentwl,zeroangleoffset        !Incident wavelength and zero angle offset (only used for reactor instruments)

    logical subtractsampledependentbackground,haved4ceffs     
    
! Paths to extract data from nexus files

    character(len=256)            :: titlepath     ! Path to run title
    character(len=256)            :: usernamepath     ! Path to run title
    character(len=256)            :: starttimepath     ! Path to run title
    character(len=256)            :: tcbpath       ! Path to time channel boundaries
    character(len=256)            :: goodframespath       ! Path to number of goodframes
    character(len=256)            :: muamphrspath       ! Path to number of micro-amp-hrs
    logical                       :: periodspresent        ! '.true.' or '.false.' depending on whether periods are present
! If true it is assumed the period will be stored in the last rank of the monitor and detector arrays
! The number of periods is obtained from the same value
    character(len=256), dimension(:), allocatable   :: spectrumcountspath     ! Paths to spectrum counts data
    integer                                         :: nspectrumpaths,mspectrumpaths
    character(len=256), dimension(:), allocatable   :: spectrummodulevalues
    integer, dimension(:), allocatable              :: firstspectrumthispath,lastspectrumthispath,dsetlineref
CONTAINS

!***********************************************************************************
!*
!*    get_run_par.for
!*
!*    Copyright A K Soper, January 2003
!*     derived from the ISIS GET routines by Freddie Akeroyd
!*
!*    This program reads the run parameters stored in the header 
!* block of a raw file and sets up the time channel boundaries
!*
!***********************************************************************************
    subroutine get_run_par(runno, nrunno)

        USE inputfilestrings
        USE reallocation_routines
        USE bad_detectors
!
! internal variables
!c
        character(len=256), dimension(:), intent(in) :: runno            !filename to obtain calibration
        integer, intent(in)                          :: nrunno
        character(len=256) :: fname        !name of file to read
        INTEGER :: i,j,id,is,iflag
        integer :: mmodnew            !no. of characters in directory name
        integer :: nuse,errcode !dummy variables
        integer, dimension(1)                   :: idumaks         !dummy variables
        logical :: found
!c
!c determine number of characters in filename
!c
        nchar=len_trim(runno(1))
!c
!c generate filename to read
!c
        fname=direct(1:lendirect)//runno(1)
        write(6,*) 'get_run_par> ',fname(1:len_trim(fname))

!Remove module arrays

        if(allocated(crateno)) deallocate(crateno)
        if(allocated(slotno)) deallocate(slotno)

!c If instrument name is D4C we need an alternative program to
!c to setup the run parameters. In particular we will need to read ALL
!c the specified data files to get the complete detector array.

        if(inst(1:len_trim(inst)).eq.'D4C') then

            call get_run_par_D4C(runno,nrunno)    
!c
!c set up detector numbers for each spectrum. Where more than one detector 
!c refers to the same spectrum, the last detector referred to will be the one 
!c used
!c
            call reallocate1d_i(specno,ndet)
            call reallocate1d_i(detno,nspec)
            call reallocate1d_i(imodule,ndet)
            do id=1,ndet
                is=specno(id)
                specno(id)=is
                detno(is)=id
            end do
!c
!c now set up the modules, namely all the spectra associated with a particular DAE
!c slot are treated as coming from 1 module
!c
            nmod=0
            mmod=0
            do i=1,ndet
!c
!c get the spectrum number for this detector
!c
                iflag=0
                j=0
                do while(iflag.eq.0.and.j.lt.nmod)
                    j=j+1
                    if(crateno(j).eq.crate(i).and.slotno(j).eq.slot(i)) then
                        imodule(i)=j
                        ndetpermodule(j)=ndetpermodule(j)+1
                        iflag=j
                    endif
                end do
                if(iflag.eq.0) then
                    nmod=nmod+1
                    do while (nmod.gt.mmod)
                        mmod=mmod+11
                        call reallocate1d_i(ndetpermodule,mmod)
                        call reallocate1d_i(crateno,mmod)
                        call reallocate1d_i(slotno,mmod)
                    end do
                    imodule(i)=nmod
                    crateno(nmod)=crate(i)
                    slotno(nmod)=slot(i)
                    ndetpermodule(nmod)=1
                end if
            end do
 
        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

            call get_run_par_nxs(fname)
            !No of spectra in file equals number of spectra
            nspecb=nspec
            
        else
!c
!c open the ISIS raw file
!c
            call open_data_file(fname,nchan,ndet,nuse,errcode)
            if(errcode.ne.0) then
                write(6,*) 'Error code',errcode &
                ,' when opening file',fname(1:lendirect+nchar)
                stop
            end if
            call getparc(fname,'HDR ',info,1,nuse,errcode)
            call getpari(fname,'NPER',idumaks,1,nuse,errcode)
            nper=idumaks(1)
            call getpari(fname,'NSP1',idumaks,1,nuse,errcode)
            nspec=idumaks(1)
            call getpari(fname,'NDET',idumaks,1,nuse,errcode)
            ndet=idumaks(1)

!c For calculation of the starting spectrum number for datasets with more than
!c one period we need to take account of the invisible spectrum zero

            nspecb=nspec+1

            ! Allocate the spectrum and detector arrays
            call reallocate1d_i(specno,ndet)
            call reallocate1d_i(detno,nspec)
            call reallocate1d_i(crate,ndet)
            call reallocate1d_i(slot,ndet)
            call reallocate1d_i(imodule,ndet)
            call reallocate1d_i(udet,ndet)
!c
!c spectrum numbers for detectors
!c
            call getpari(fname,'SPEC',specno,ndet,nuse,errcode)
!c
!c get the crate numbers for each detector
!c
            call getpari(fname,'CRAT',crate,ndet,nuse,errcode)
!c
!c get the slot numbers for each detector. For the purposes of deadtime correction all
!c detectors going into a DAE slot will be treated as having been multiplexed with
!c a common deadtime, mdeadtime.
!c
            call getpari(fname,'MODN',slot,ndet,nuse,errcode)
            call getpari(fname,'UDET',udet,ndet,nuse,errcode)


            call close_data_file()
            write(6,*) 'get_run_par> Closed data file'
!c
!c set up detector numbers for each spectrum. Where more than one detector 
!c refers to the same spectrum, the last detector referred to will be the one 
!c used to get the calibration.
!c
            do id=1,ndet
                is=specno(id)
                if(is.lt.1.or.is.gt.nspec) then
                    write(6,*) is,nspec
                else
                    detno(is)=id
                end if
                specno(id)=is
            end do
            write(6,*) 'get_run_par> Set up detector and spectrum numbers'
            !c
!c now set up the modules, namely all the spectra associated with a particular DAE
!c slot are treated as coming from 1 module
!c
            nmod=0
            mmod=0
!            allocate(crateno(mmod),slotno(mmod),ndetpermodule(mmod))
            do i=1,ndet
!c
!c get the spectrum number for this detector
!c
                found=.false.
                j=0
                do while(.not.found.and.j.lt.nmod)
                    j=j+1
!                    write(6,*) 'get_run_par> ',i,j,nmod
                    found=(crateno(j).eq.crate(i).and.slotno(j).eq.slot(i))
                    if(found) then
                        imodule(i)=j
                        ndetpermodule(j)=ndetpermodule(j)+1
                    endif
                end do
                if(.not.found) then
                    !Create a new module
                    nmod=nmod+1
!                    write(6,*) 'get_run_par> ',i,nmod,mmod,size(imodule),size(crateno),size(slotno)
                    do while (nmod.gt.mmod)
                        mmod=mmod+10
                        call reallocate1d_i(ndetpermodule,mmod)
                        call reallocate1d_i(crateno,mmod)
                        call reallocate1d_i(slotno,mmod)
                    end do
                    imodule(i)=nmod
                    crateno(nmod)=crate(i)
                    slotno(nmod)=slot(i)
                    ndetpermodule(nmod)=1
                end if
            end do
        endif
!c
!c write out the parameters read from header block as a check
!c
        open(10,file='gudrun_run_par.dat',status='unknown')
        WRITE(10,'(A)') fname
        WRITE(10,*)' No. of periods = ', nper
        WRITE(10,*)' No. of spectra = ', nspec
        WRITE(10,*)' No. of detectors = ', ndet
        write(10,*)' No. of modules = ', nmod
        write(10,*)' No. of detectors per module'
        write(10,'(10(1x,i5))') (ndetpermodule(i),i=1,nmod)
        write(10,*) ' Detector spectrum numbers:- '
        write(10,'(10(1x,i5))') (specno(i),i=1,ndet)
        write(10,*) ' Spectrum detector numbers:- '
        write(10,'(10(1x,i5))') (detno(i),i=1,nspec)
        write(10,*)' Crate number of each detector'
        write(10,'(10(1x,i5))') (crate(i),i=1,ndet)
        write(10,*)' Slot number of each detector'
        write(10,'(10(1x,i5))') (slot(i),i=1,ndet)
        write(10,*)' Module number of each detector'
        write(10,'(10(1x,i5))') (imodule(i),i=1,ndet)
        close (10)
        return

    end subroutine get_run_par
    
    subroutine get_run_par_d4c(runno,nrunno)

        USE inputfilestrings
        USE reallocation_routines
        USE bad_detectors

        !c Program to read the D4C raw files and set up equivalent parameters for
!c the input to Gudrun 

        character(len=256),dimension(:), intent(in) :: runno !Filenames to be inspected
        integer, intent(in) :: nrunno
        character(len=256) sometext,fname
        integer :: i,j,ierr,ispec,ibadread,imod,idet,iflag,nmon,nspecmon
        real :: angarm,effd4ctemp
        integer :: foundformat,nskip,nentries


!c For D4C we need to read all the files to determine the correct number 
!c of spectra for this run. For fixed wavelength data each spectrum will
!c have precisely 1 time-of-flight channel, and 2 time channel boundaries.

        nper=1

!c For D4C we will process exactly 578 spectra at a time, since this is the number
!c contained in each raw file. 
! The first spectrum is the monitor reading for each set of modules, and the
! second spectrum in the length of the run in milliseconds

        nproc=2+9*64

!        write(6,*) 'get_run_par_D4C> 1'
        call reallocate2d_r(effd4cinput,64,9)
! Set the initial efficiencies in case the file 'effd4c.eff' cannot be opened
! The first 3 and last 3 spectra of each block of 64 will be disregarded due
! to edge effects.
        imod=0
        do while (imod.lt.9)
            imod=imod+1
            idet=0
            do while (idet.lt.64)
                idet=idet+1
                if(idet.lt.4.or.idet.gt.61) then
                    effd4cinput(idet,imod)=-1.0
                else
                    effd4cinput(idet,imod)=1.0
                endif
            end do
        end do
! First try to read the file 'effd4c.eff'
        open(10,file='effd4c.eff',status='old',iostat=ierr)
        do while(ierr.eq.0)
            read(10,'(a)',iostat=ierr) line
!Ignore lines containing a #
            if(ierr.eq.0.and.index(line,'#').eq.0) then
!Parse the line. It must contain at least three words to be useful, module number, detector number and efficiency
                call parse()
                if(nwords.gt.2) then
                    read(line,*,iostat=ierr) imod,idet,effd4ctemp
                    if(ierr.eq.0) then
                        if(effd4ctemp.gt.0.0) then
                            effd4cinput(idet,imod)=effd4ctemp
                        else
                            effd4cinput(idet,imod)=0.0
                        end if
                    end if
                end if
            end if
        end do
        close(10)
        do imod=1,9
            do idet=1,64
                write(6,*) imod,idet,effd4cinput(idet,imod)
            end do
        end do
!        write(6,*) 'get_run_par_D4C> 2'

!c Input module counter

        nmod=0
        mmod=0

!c Monitor counter 

        nmon=0
        nspecmon=1
        ndet=0
        mdet=0

!c Open each file in turn and read the relevant data

        do i=1,nrunno

!c Determine number of characters in filename

            nchar=len_trim(runno(i))

!c generate filename to read

            fname=direct(1:lendirect)//runno(i)(1:nchar)

!c open the raw file

            open(10,file=fname,status='old',iostat=ierr)
            if(ierr.ne.0) write(6,'(a,1x,i6,1x,a)') 'get_run_par_D4C> Error code', ierr &
     ,' when opening file',fname(1:len_trim(fname))

!c Search for the word 'MonitorCnts'
            sometext=' '
            foundformat=0
            do while(foundformat.eq.0.and.ierr.eq.0)
               read(10,'(a)',iostat=ierr) sometext
               foundformat=index(sometext,'MonitorCnts')
            end do
            if(ierr.ne.0) then
               write(6,102) ierr,fname(1:len_trim(fname))
102   format(/'get_run_par_D4C> Error code ',i6,' when reading file: ',a)
               stop
            endif
            nskip=9
!c Skip nskip lines
            do j=1,nskip
               read(10,*)
            end do
!c Read the angular position of the detector arm
            read(10,*) angarm,angarm

            write(6,*) i,' ',fname(1:lendirect+nchar),' ',angarm

!c Step through the existing modules if any and check that this angle
!c does not already exist. If not then generate a new set of modules

            imod=0
            iflag=0
                
            do while(imod.lt.nmod.and.iflag.eq.0)
               imod=imod+1
               if(abs(moduleangle(imod)-angarm).lt.0.0001) then
                   iflag=imod
               endif
            end do
            if(iflag.eq.0) then

!c A new set of detectors have been found. Generate the modules for these detectors

                ! Create 11 new modules
                
                mmod=nmod+11
                call reallocate1d_r(moduleangle,mmod)
                mdet=ndet+nproc
                call reallocate1d_i(specno,mdet)
                call reallocate1d_i(crate,mdet)
                call reallocate1d_i(slot,mdet)
                call reallocate1d_r(effd4c,mdet)
                call reallocate1d_i(ibad,mdet)

!c Monitors are at 180 deg.

                nmod=nmod+1
                moduleangle(nmod)=180.0
!c Spectrum numbers for the incident monitor
                ndet=ndet+1
                specno(ndet)=ndet
                crate(ndet)=nspecmon
                slot(ndet)=1
                effd4c(ndet)=1.0
                ibad(ndet)=-99
! Length of run will be stored in the second module of each set of 11. Set an angle of zero so this module won't
! be detected as a scattering detector.
                nmod=nmod+1
                moduleangle(nmod)=360 !Since we don't use this detector for scattering put in a crazy angle that is unlikely to be found
! by the real scattering detectors
!c Spectrum numbers for the length of run in milliseconds
                ndet=ndet+1
                specno(ndet)=ndet
                crate(ndet)=nspecmon
                slot(ndet)=2
                effd4c(ndet)=1.0
                ibad(ndet)=-1
                imod=0
                do while (imod.lt.9)
                    imod=imod+1
                    nmod=nmod+1
                    moduleangle(nmod)=angarm
                    idet=0
                    do while(idet.lt.64)
                        idet=idet+1
                        ndet=ndet+1
                        specno(ndet)=ndet
                        !We ignore the first 3 and last 3 detectors in each block
                        if(idet.lt.4.or.idet.gt.61) ibad(specno(ndet))=-1
!c The crate number will be used to signify which monitor corresponds to this detector
                        crate(ndet)=nspecmon
!c Each module is assumed to go into a different slot for the specified crate
                        slot(ndet)=imod+2
!c Assign the D4C efficiency factor and set the bad detector status
                        effd4c(ndet)=effd4cinput(idet,imod)
                        if(effd4c(ndet).gt.0.0) then
                            ibad(ndet)=0
                        else
                            ibad(ndet)=-1
                        end if
                    end do
                    angarm=angarm+15.0
                end do
                nspecmon=nspecmon+nproc
            endif
            close(10)
        end do
        do imod=1,nmod
            write(6,*) imod,moduleangle(imod)
        end do
        do idet=1,ndet
            write(6,*) specno(idet),crate(idet),slot(idet),ibad(idet)
        end do
        nspec=ndet
        nspecb=nspec

!c Create and write the bad detector file as a record

        open(10,file='spec.bad',status='unknown',iostat=ierr)
        do ispec=1,nspec
           write(10,*) ispec,ibad(ispec)
        end do
        close(10)

        return
        
    end subroutine get_run_par_D4C
        
    subroutine get_run_par_nxs(filename)

        USE inputfilestrings
        USE datasetlist
        USE reallocation_routines
        USE bad_detectors

        character(len=256), intent(in)      :: filename    ! Name of data file
        integer                             :: i,j,itype
        integer                             :: iprod,iprodmax,icountpaths
        integer                             :: imref,imrefstart,ispecref,iref,itimes,ntimes
        logical                             :: found
        integer                             :: nmodulevalues,lastspectrum
        integer, dimension(:), allocatable  :: modulevalues

        call retrieve_list(filename)
      
        if(error.eq.0) then
          
            write(6,'(a,1x,i4)') 'Number of datasets in file = ',nmembers
            write(6,'(a,a)') 'Dataset number, name, class (0, 1 or 3)'&
            ,'), precision (no. of characters for class = 3),ndims and dims.' 
            do i=1,nmembers
                write(6,'(1x,i4,1x,a)') i,dsetline(i)(1:len_trim(dsetline(i)))
            end do
            
        else
            
            stop
            
        end if
      
!c Get the number of spectra. With nexus files this will be identical to the number of detectors
!c since in general the grouping of detectors into a particular spectrum will not be available outside
!!c of ISIS. In any case it is the measured spectra that count towards the data analysis, not the
!c detectors.

!c Spectra are counted in the order in which spectra are entered in the <inst>_nexus.txt file.

!c As the spectra are read in, so the module number is set according to the value in spectrummodulevalue 
!c for each data block

        nspec=0
        mspec=0
        nmod=0
        mmod=0
        icountpaths=0
        lastspectrum=0
        call reallocate1d_i(firstspectrumthispath,nspectrumpaths)
        call reallocate1d_i(lastspectrumthispath,nspectrumpaths)
        do while (icountpaths.lt.nspectrumpaths)
            icountpaths=icountpaths+1
            found=.false.
            !Find this spectrum path in the dataset list
            found=retrieve_dims(filename,spectrumcountspath(icountpaths)).gt.0
            if(found) then
                if((periodspresent.and.ndims.gt.1).or.(.not.periodspresent.and.ndims.gt.0)) then

!c The number of spectra will the product of dimensions for ranks 2 or more. If only rank 1 is
!c present then this path contributes just one spectrum. If periods are present, the largest rank
!c is the period number and so should not be used to assess the number of spectra.

                    iprod=1
                    iprodmax=ndims
                    if(periodspresent) then
                        iprodmax=ndims-1
                        nper=dims(ndims)
                    else
                        nper=1
                    endif
!   
!c First dimension is the counts versus time channel array
!c Subsequent dimensions represent the spatial pixellation within a detector element, and the detector elements themselves.
!
                    i=1
                    do while(i.lt.iprodmax)!
                        i=i+1
                        iprod=iprod*dims(i)!
                    end do
                    firstspectrumthispath(icountpaths)=lastspectrum
                    lastspectrumthispath(icountpaths)=lastspectrum+iprod
                    lastspectrum=lastspectrumthispath(icountpaths)
                    write(6,'(a,a,2(1x,i6))') 'First and last spectra for path: ',&
                    spectrumcountspath(icountpaths)(1:len_trim(spectrumcountspath(icountpaths))) &
                    ,firstspectrumthispath(icountpaths),lastspectrumthispath(icountpaths)

!c iprod is the number of spectra in this counts array. These need to be grouped in modules 
!c as specified in the instrument file. If no module grouping is specified in the instrument file
!c then all these spectra are assumed to belong to one module

                    call parse(spectrummodulevalues(icountpaths))
                    nmodulevalues=nwords
!                    write(6,'(i4,1x,a)') nmodulevalues,spectrummodulevalues(icountpaths)(ncf(1):ncl(nwords))
                    if(allocated(modulevalues)) deallocate(modulevalues)

                    if(nmodulevalues.gt.0) then

                        allocate(modulevalues(nmodulevalues))
                        do j=1,nmodulevalues
                            read(spectrummodulevalues(icountpaths)(ncf(j):ncl(j)),*) modulevalues(j)
                        end do
                        imref=1
                        ntimes=modulevalues(imref)
                        itimes=1
                        imref=imref+1
                        imrefstart=imref
                        nmod=nmod+1
                        do while(nmod.gt.mmod)
                            ! Need to increment array sizes
                            mmod=mmod+10
                            call reallocate1d_i(crateno,mmod)
                            call reallocate1d_i(slotno,mmod)
                            call reallocate1d_i(ndetpermodule,mmod)
                        end do
!c Increment the module number and save the number of detectors in this module
                        ndetpermodule(nmod)=modulevalues(imref)
                        !The next two values are assigned for completeness, but I don't think they are used
                        !elsewhere in the program.
                        crateno(nmod)=nmod
                        slotno(nmod)=nmod
                        i=0
                        iref=0
                        do while (i.lt.iprod)
                            if(iref.ge.ndetpermodule(nmod)) then

                                iref=0
                                imref=imref+1
!c If the series is terminated incorrectly (not ending on a negative number) then the last module is repeated indefinitely
                                if(imref.gt.nmodulevalues) then
                                    imref=imref-1
                                endif

!c If the next spectrummodulevalue is -ve, we step back the corresponding number of modules
!c unless itimes is equal to ntimes, in which case we move on to the next set of values

                                if(modulevalues(imref).lt.0) then
                                    if(itimes.lt.ntimes.or.ntimes.eq.0) then
                                        itimes=itimes+1
                                        imref=imrefstart
                                    else
                                        imref=imref+1
                                        ntimes=modulevalues(imref)
                                        itimes=1
!c Since the first value of the next batch is the number of times this sequence is repeated, we need to step to the
!c next value
                                        imref=imref+1
                                        imrefstart=imref
                                    endif
                                endif

!c Increment module number
                                nmod=nmod+1
                                do while (nmod.gt.mmod)
                                    ! Need to increment array sizes
                                    mmod=mmod+10
                                    call reallocate1d_i(crateno,mmod)
                                    call reallocate1d_i(slotno,mmod)
                                    call reallocate1d_i(ndetpermodule,mmod)
                                end do
                                ndetpermodule(nmod)=modulevalues(imref)
                                crateno(nmod)=nmod
                                slotno(nmod)=nmod

                            endif
                            i=i+1
                            iref=iref+1
                            ispecref=nspec+i
                            do while (ispecref.gt.mspec)
                                call reallocate1d_i(imodule,ispecref+10)
                                mspec=ispecref+10
                            end do
                            imodule(ispecref)=nmod

!      write(6,*) ispecref,i,iref,imref,itimes,ntimes,nmod

                        end do

                    else

!c All spectra are put into one module

                        nmod=nmod+1
                        do while(nmod.gt.mmod)
                            ! Need to increment array sizes
                            mmod=mmod+10
                            call reallocate1d_i(crateno,mmod)
                            call reallocate1d_i(slotno,mmod)
                            call reallocate1d_i(ndetpermodule,mmod)
                        end do
                        ndetpermodule(nmod)=iprod
                        do while(i.lt.iprod)
                            i=i+1
                            ispecref=nspec+i
                            do while (ispecref.gt.mspec)
                                call reallocate1d_i(imodule,ispecref+10)
                                mspec=ispecref+10
                            end do
                            imodule(ispecref)=nmod
                        end do

                    end if
                    nspec=nspec+iprod
                    
                else
 
                    write(6,'(1x,a,a,a)') 'Dataset: '&
                    ,spectrumcountspath(icountpaths)(1:len_trim(spectrumcountspath(icountpaths)))&
                    ,' has insufficient dimensions to allocate spectra'
                
                endif

            else

                write(6,'(1x,a,a,a)') 'Dataset: '&
                ,spectrumcountspath(icountpaths)(1:len_trim(spectrumcountspath(icountpaths)))&
                ,' not found'
                
            endif

        end do

        ndet=nspec
        call reallocate1d_i(specno,ndet)
        call reallocate1d_i(detno,ndet)
        call reallocate1d_i(crate,ndet)
        call reallocate1d_i(slot,ndet)
!c
!c set up detector numbers for each spectrum. We retain the original ISIS form so the program remains
!c compatible with the old ISIS raw files.
!c
        do i=1,nspec
                specno(i)=i
                detno(i)=i
                crate(i)=crateno(imodule(i))
                slot(i)=slotno(imodule(i))
        end do

        return

    end subroutine get_run_par_nxs    
    
!***********************************************************************************
!*
!*	get_tcbv.for
!*
!*	A K Soper, January 2003
!* 	derived from ISIS GET routines by Freddie Akeroyd	
!*
!*	This program reads the time channel boundaries for a specified run no.
!*
!***********************************************************************************
    subroutine get_tcbv(runno)

        USE reallocation_routines
!c
!c internal variables
!c
        character(len=256), intent(in)          :: runno                !filename to obtain calibration
        character(len=256)                      :: fname                !name of file to read
        integer                                 :: i,ierr
        integer, dimension(:), allocatable      :: cpb                !clock pulse boundaries (clock pulses)
        integer                                 :: pre                        !clock prescaler
        integer, dimension(64)                  :: daep                !used to get dae offset
        integer                                 :: offset                !DAE offset in clock pulses
        integer                                 :: nchar                        !no. of characters in directory name
        integer                                 :: nuse,errcode !dummy variables
        integer, dimension(1)                   :: idumaks         !dummy variables
        real                                    :: div,mid                        !used by routine for TCB generation

!c determine number of characters in filename

	nchar=index(runno,' ')-1
!c
! generate filename to read
!c
	fname=direct(1:lendirect)//runno(1:nchar)

!c For reactor instruments we can generate the calibration from data in the header 
!c blocks

	if(inst(1:index(inst,' ')-1).eq.'D4C') then

!c Only 1 time channel

            nchan=1
            nchanb=nchan+1
            call reallocate1d_i(cpb,nchanb)
            call reallocate1d_r(tcb,nchanb)
            call reallocate1d_r(tcw,nchan)
            offset=0.0
            div=1.0e-2
	    mid=incidentfp*incidentwl/0.0039554
!            mid=10.0/0.0039554
            cpb(1)=int((mid-0.5)/div)
            cpb(2)=int((mid+0.5)/div)
            write(6,*) 'get_tcbv> ',mid,cpb(1),cpb(2)
!c
!c calculate time channel boundaries and widths
!
            i=0
            do while(i.lt.nchanb)
                i=i+1
		tcb(i)=real(cpb(i)+offset)*div
            end do


	else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

            call get_tcb_nxs(fname)

!c Save the time channel boundaries and time channel widths

	else 
!c
!c open the raw file
!c
            call open_data_file(fname,nchan,ndet,nuse,errcode)
            if(errcode.ne.0) write(6,*) 'Error code',errcode,' when opening file',fname
            call getpari(fname,'NTC1',idumaks,1,nuse,errcode)
            nchan=idumaks(1)
            nchanb=nchan+1
            call reallocate1d_i(cpb,nchanb)
            call reallocate1d_r(tcb,nchanb)
            call reallocate1d_r(tcw,nchan)
            call getpari(fname,'TCB1',cpb,nchanb,nuse,errcode)
            call getpari(fname,'PRE1',idumaks,1,nuse,errcode)
            pre=idumaks(1)
            call getpari(fname,'DAEP',daep,64,nuse,errcode)
            offset=daep(24)*4*32		!offset in clock pulses (assumes 32MHz)
            div=real(pre)/real(32)
!c
!c close the raw file
!c
            call close_data_file()
!c
!c calculate time channel boundaries and widths
!c
            i=0
            do while(i.lt.nchanb)
                i=i+1
		tcb(i)=real(cpb(i)+offset)*div
            end do

	endif

	if(nchanv.eq.0) nchanv=nchan
	if(nchan.ne.nchanv) then
            write(6,*) 'No. of time channel boundaries for run ',runno,' is not equal to ',nchanv
            stop
	else
            if(.not.allocated(tcbv)) then
                allocate(tcbv(nchanb),tcwv(nchan))
                i=0
                do while(i.lt.nchanb)
                    i=i+1
                    tcbv(i)=tcb(i)
                end do
                i=0
                do while(i.lt.nchan)
                    i=i+1
                    tcwv(i)=tcb(i+1)-tcb(i)
                end do
                nchanv=nchan
                nchanbv=nchanv+1
            end if
	endif
        nspchv=nspec*nchanbv
	open(10,file='gudrun_van_tcb.dat',status='unknown')
	WRITE(10,'(A)') fname
	WRITE(10,*)' No. of time channels = ', nchanv
	WRITE(10,*)' No. of time channel boundaries = ', nchanbv
	WRITE(10,*)' Time channel boundaries (musec) :- '
	write(10,'(5(1x,e13.6))') (tcbv(I),I=1,nchanbv)
	WRITE(10,*)' Time channel widths (musec) :- '
	write(10,'(5(1x,e13.6))') (tcwv(I),I=1,nchanv)
	close (10)
        if(allocated(cpb)) deallocate(cpb)
        if(allocated(tcb)) deallocate(tcb)
        if(allocated(tcw)) deallocate(tcw)
        
	return
    end
        
    subroutine get_tcbs(runno)
!***********************************************************************************
!*
!*        get_tcbs.for
!*
!*        A K Soper, January 2003
!*        Uses FAA stand alone get routines
!*
!*        This program reads the time channel boundaries for a specified run no.
!*
!***********************************************************************************

        USE reallocation_routines
!       include 'calibration.inc'
!c
!c internal variables
!c
        character(len=256), intent(in)          :: runno                !filename to obtain calibration
        character(len=256)                      :: fname                !name of file to read
        integer                                 :: i,ierr
        integer, dimension(:), allocatable      :: cpb                !clock pulse boundaries (clock pulses)
        integer                                 :: pre                        !clock prescaler
        integer, dimension(64)                  :: daep                !used to get dae offset
        integer                                 :: offset                !DAE offset in clock pulses
        integer                                 :: nchar                        !no. of characters in directory name
        integer                                 :: nuse,errcode !dummy variables
        integer, dimension(1)                   :: idumaks         !dummy variables
        real                                    :: div,mid                        !used by routine for TCB generation

!c determine number of characters in filename

        nchar=len_trim(runno)
!c
!c generate filename to read
!c
        fname=direct(1:lendirect)//runno(1:nchar)

!c For reactor instruments we can generate the calibration from data in the header 
!c blocks

        if(inst(1:index(inst,' ')-1).eq.'D4C') then

!c Only 1 time channel

            nchan=1
            nchanb=nchan+1
            call reallocate1d_i(cpb,nchanb)
            call reallocate1d_r(tcb,nchanb)
            call reallocate1d_r(tcw,nchan)
            offset=0.0
            div=1.0e-2
	    mid=incidentfp*incidentwl/0.0039554
            cpb(1)=int((mid-0.5)/div)
            cpb(2)=int((mid+0.5)/div)
!c
!c calculate time channel boundaries and widths
!c
            i=0
            do while(i.lt.nchanb)
                i=i+1
                tcb(i)=real(cpb(i)+offset)*div
            end do

        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

            call get_tcb_nxs(fname)

        else
!c
!c open the raw file
!c
            call open_data_file(fname,nchan,ndet,nuse,ierr)
            if(ierr.ne.0) write(6,*) 'Error code',ierr,' when opening file',fname
            call getpari(fname,'NTC1',idumaks,1,nuse,errcode)
            nchan=idumaks(1)
            nchanb=nchan+1
!            write(6,*) 'get_tcbs> ',nchan,nchanb
            call reallocate1d_i(cpb,nchanb)
            call reallocate1d_r(tcb,nchanb)
            call reallocate1d_r(tcw,nchan)
            call getpari(fname,'TCB1',cpb,nchanb,nuse,errcode)
            call getpari(fname,'PRE1',idumaks,1,nuse,errcode)
            pre=idumaks(1)
            call getpari(fname,'DAEP',daep,64,nuse,errcode)
            offset=daep(24)*4*32                !offset in clock pulses (assumes 32MHz)
            div=real(pre)/real(32)
            call close_data_file()
!            write(6,*) 'get_tcbs> Closed data file'
!c
!c calculate time channel boundaries and widths
!c
            i=0
            do while(i.lt.nchanb)
                i=i+1
                tcb(i)=real(cpb(i)+offset)*div
            end do
!            write(6,*) 'get_tcbs> Got time channel boundaries and widths'
            
        endif

	if(nchans.eq.0) nchans=nchan
	if(nchan.ne.nchans) then
            write(6,*) 'No. of time channel boundaries for run ',runno,' is not equal to ',nchans
            stop
	else
            if(.not.allocated(tcbs)) then
                allocate(tcbs(nchanb),tcws(nchan))
                i=0
                do while(i.lt.nchanb)
                    i=i+1
                    tcbs(i)=tcb(i)
                end do
                i=0
                do while(i.lt.nchan)
                    i=i+1
                    tcws(i)=tcb(i+1)-tcb(i)
                end do
                nchans=nchan
                nchanbs=nchanb
            end if
	endif
!        write(6,*) 'get_tcbs> Checked time channel boundaries'
        nspchs=nspec*nchanbs
        open(10,file='gudrun_sam_tcb.dat',status='unknown')
        WRITE(10,'(A)') fname
        WRITE(10,*)' No. of time channels = ', nchans
        WRITE(10,*)' No. of time channel boundaries = ', nchanbs
        WRITE(10,*)' Time channel boundaries (musec) :- '
        write(10,'(5(1x,e13.6))') (tcbs(i),i=1,nchanbs)
        WRITE(10,*)' Time channel widths (musec) :- '
        write(10,'(5(1x,e13.6))') (tcws(i),i=1,nchans)
        close (10)
        if(allocated(cpb)) deallocate(cpb)
        if(allocated(tcb)) deallocate(tcb)
        if(allocated(tcw)) deallocate(tcw)
        return
        
    end subroutine get_tcbs
    
    subroutine get_tcb_nxs(filename)

        USE reallocation_routines
        USE datasetlist

        character(len=256), intent(in)                  :: filename    ! Name of data file
        integer                                         :: i,j,myndims
        logical                                         :: found
        real(kind=4), DIMENSION(:), ALLOCATABLE, TARGET         :: rdataset1d
        real(kind=8), DIMENSION(:), ALLOCATABLE, TARGET         :: rrdataset1d
        integer, dimension(:), allocatable              :: mydims

        if(nmembers.eq.0) call retrieve_list(filename)
      
!        if(error.eq.0) then
          
!            write(6,'(a,1x,i4)') 'Number of datasets in file = ',nmembers
!            write(6,'(a,a)') 'Dataset number, name, class (0, 1 or 3)'&
!            ,'), precision (no. of characters for class = 3),ndims and dims.' 
!            do i=1,nmembers
!                write(6,'(1x,i4,1x,a)') i,dsetline(i)(1:len_trim(dsetline(i)))
!            end do
!        end if
      
!c Get the number of spectra. With nexus files this will be identical to the number of detectors
!c since in general the grouping of detectors into a particular spectrum will not be available outside
!!c of ISIS. In any case it is the measured spectra that count towards the data analysis, not the
!c detectors.

!c Spectra are counted in the order in which spectra are entered in the <inst>_nexus.txt file.

!c As the spectra are read in, so the module number is set according to the value in spectrummodulevalue 
!c for each data block

        found=.false.
        !Find the time of flight path in the dataset list
        i=0
        do while (i.lt.nmembers.and..not.found)
            i=i+1
            call parse(dsetline(i))
            found=tcbpath(1:len_trim(tcbpath)).eq.dsetline(i)(ncf(1):ncl(1))
        end do
!       write(6,*) found,nwords
        if(found) then
!           write(6,'(a)') dsetline(i)(ncf(4):ncl(nwords))
            read(dsetline(i)(ncf(4):ncl(4)),*) myndims
            !Should be a 1 dimensional dataset
            if(myndims.ne.1) then
                    
                write(6,'(a)') 'Incorrect number of dimensions for time-of-flight array. Should be 1 dimensional.'
                stop
                    
            end if
            !Can't use reallocate_routines here, because dims has kind defined by HDF5
            call reallocate1d_i(mydims,myndims)
            read(dsetline(i)(ncf(5):ncl(nwords)),*) (mydims(j),j=1,myndims)
            nchanb=mydims(1)
            nchan=nchanb-1
            !Get this block of data
            call retrieve_data_r(filename,tcbpath,rdataset1d,rrdataset1d)
            call reallocate1d_r(tcb,nchanb)
            if(inst(1:len_trim(inst)).eq.'NOVA') then
                do i=1,nchanb
                    tcb(i)=real(rrdataset1d(i))
                end do
            else
                do i=1,nchanb
                    tcb(i)=rdataset1d(i)
                end do
            endif
        else
            write(6,'(a)') 'Path to time-of-flight array not found.'
            stop
        end if
        return
        
    end subroutine get_tcb_nxs
    
    function convertfilename(text)

        character(len=256), intent(in)      :: text
        character(len=256)                  :: convertfilename

        convertfilename = gudrunstartupfolder(1:ilengudrunstartupfolder) &
        //pathseparator(1:lenpathseparator)//text

!      write(6,100) text(1:80)
!      write(6,100) convertfilename(1:80)
100   format(1x,a)
        return
    end function convertfilename

END MODULE run_par

