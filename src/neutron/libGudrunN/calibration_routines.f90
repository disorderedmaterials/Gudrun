!     
! File:   calibration.f90
! Author: aks45
!
! Created on 27 April 2012, 12:56
!

MODULE calibration_routines
    
    implicit none
!c
!c global calibration variables common to all routines
!c
    integer                         :: phiutno        !User table column containing phi values
    real                            :: lenin            !incident flight path
    real, dimension(:), allocatable :: deltat        !time offset in musec for each detector
    real, dimension(:), allocatable :: lendet        !scattered flight path
    real, dimension(:), allocatable :: ttheta        !detector scattering angle
    real, dimension(:), allocatable :: phi        !detector azimuthal angles
    character(len=256)              :: fnamed !name of detector calibration file
     
    
    CONTAINS

    subroutine get_calibration(runno,nrunno)
!***********************************************************************************
!*
!*    get_calibration.for
!*
!*    Copyright A K Soper, January 2003
!*     derived from the ISIS GET routines by Freddie Akeroyd
!*
!*    This program reads the calibration information stored in the header 
!* block of a raw file and sets up the time channel boundaries
!*
!***********************************************************************************

        USE reallocation_routines
        use inputfilestrings
        USE run_par

        character(len=256), dimension(:), intent(in)    :: runno
        integer, intent(in)                             :: nrunno !Number of run files supplied
   
!    internal variables

        character(len=256)                              :: fname        !name of file to read
        character(len=256)                              :: extension          !Filename extension
        character(len=4)                                :: utn            !label for user table containing phi values
        integer                                         :: i,j,is,id,it,iflag,imod,idet,incrementref,im,ierr
        integer                                         :: nchard        !no. of characters in detector filename
        integer                                         :: ndetd            !no. of detectors in det. calib.  file
        integer                                         :: nutno            !no. of user table entries in detector.dat
        integer                                         :: nuse,errcode        !dummy variables
        integer, dimension(1)                           :: idumaks         !dummy variables
        real, dimension(64)                             :: ivpbf        !used to obtain ivpb(23)
        real                                            :: a,b,c,d,e            !dummy variables
        real, dimension(11)                             :: increment     !incremental angles for each module of detectors (D4C)
        character(len=1)                                :: ans
        logical test
!c
!c set up label for user table to be read for phi values
!c
        if(phiutno.gt.0.and.phiutno.lt.10) then
            write(utn,102) phiutno
102         format('UT',i1,' ')
        else if(phiutno.gt.9.and.phiutno.lt.100) then
            write(utn,103) phiutno
103         format('UT',i2)
        else
            phiutno=0
            write(utn,102) phiutno
        endif
        write(6,*) phiutno,utn
        !Allocate the arrays to be used to store the calibration
        call reallocate1d_r(deltat,ndet)
        call reallocate1d_r(lendet,ndet)
        call reallocate1d_r(ttheta,ndet)
        call reallocate1d_r(phi,ndet)
!c
!c determine number of characters in calibration file - if specified
!c
        nchard=index(fnamed,'.')

!c For D4C we can generate the calibration from data in the header 
!c blocks

        if(inst(1:index(inst,' ')-1).eq.'D4C') then
    
!c Set up some numbers specific to D4C

!c Read the angular increments for the nine detector modules on D4C

            do i=1,11
                increment(i)=0.0
            end do
            open(10,file='dec.dec',status='old',iostat=ierr)
            do while(ierr.eq.0)
                read(10,'(a)',iostat=ierr) line
                if(ierr.eq.0.and.index(line,'#').eq.0) then !Ignore lines containing #
                    call parse()
                    if(nwords.gt.1) then
                        read(line,*) im,increment(im+2) !First two values are
! monitor and run time and so are left at zero.
                        increment(im+2)=increment(im+2)+zeroangleoffset !Adjust for zero angle offset.
                    end if
                end if
            end do
            close(10)

            if(incidentfp.gt.0.0) then
                lenin=incidentfp
            else
                write(6,698)
698    format(/'get_calibration> Cannot have incident flight path = 0' &
        ,' for reactor instrument')
                lenin=10.0
            endif

!c Step through the modules and detectors and set up the scattering angles, phi values
!c secondary flight path (0) and time offset.

            do i=1,ndet

                deltat(i)=0.0
                lendet(i)=0.0
                phi(i)=0.0

            end do

            idet=0
            incrementref=0

!c Monitors are at 180deg. 

            do imod=1,nmod
                incrementref=incrementref+1
                if(incrementref.gt.11) incrementref=1
                write(6,*) imod,incrementref,ndetpermodule(imod),moduleangle(imod),increment(incrementref)
                do id=1,ndetpermodule(imod)
                    idet=idet+1
                    ttheta(idet)=moduleangle(imod)+increment(incrementref)+real(id-1)*0.125
                end do
            end do

        else if(ext.eq.'nxs'.or.ext.eq.'NXS') then

!c For NeXus files, the calibration MUST be read from a .calib file

            test=nchard.gt.0
            if(test) then

!c Check the extension of the detector filename

                test=fnamed(nchard+1:nchard+5).eq.'CALIB'.or.fnamed(nchard+1:nchard+5).eq.'calib'

            endif
            if(.not.test) then

                write(6,'(a)') 'NeXus files must have a .calib calibration file'

                stop

            endif
            !Check incident flight path has been defined
            if(incidentfp.gt.0.0) then
                lenin=incidentfp
            else
                write(6,*) 'get_calibration> For NeXus files incident flight path has to be defined!'
                stop
            end if

        else

!c determine number of characters in filename

            nchar=index(runno(1),' ')-1

!c generate filename to read

            fname=direct(1:lendirect)//runno(1)(1:nchar)
            CALL open_data_file(fname,nchan,ndet,nuse,errcode)

!c data for incident flight path from header block

            call getparr(fname,'RVPB',ivpbf,64,nuse,errcode)

            if(incidentfp.le.0.0) then
                LENIN=ivpbf(23)
            else
                lenin=incidentfp
            endif

!c Initially read the rest of calibration from the header

            call getparr(fname,'DELT',deltat,ndet,nuse,errcode)
            call getparr(fname,'LEN2',lendet,ndet,nuse,errcode)
            call getparr(fname,'TTHE',ttheta,ndet,nuse,errcode)

!c get the number of User Table entries

            call getpari(fname,'NUSE',idumaks,1,nuse,errcode)
            nutno=idumaks(1)
            write(6,*) nutno,phiutno

!c get the phi values if phiutno is not zero, otherwise all phis are
!c set to zero

            if(phiutno.gt.0.and.phiutno.le.nutno) then
                call getparr(fname,utn,phi,ndet,nuse,errcode)
            else
                do id=1,ndet
                    phi(id)=0.0
                end do
            endif

!c close RAW file

            call close_data_file()
        endif

!c if nchard.gt.0 read the calibration parameters from the specified file

        if(nchard.gt.0) then
    
!c open detector calibration file and read the detector calibration

            open(10,file=fnamed,status='old',iostat=ierr)
            if(ierr.ne.0) then
                fnamed=convertfilename(fnamed)
                open(10,file=fnamed,status='old',iostat=ierr)
            endif
            
            write(6,*) 'get_calibration> Reading calibration from: ',fnamed(1:len_trim(fnamed))
            write(6,*) nchard,len_trim(fnamed)
            
            if(ierr.eq.0) then
                    write(6,*) 'get_calibration> 1'
                read(10,*)
                
!c read no detectors in this file, and number of User Table entries

                read(10,*) ndetd,nutno
                read(10,*)
                write(6,*) ndetd,nutno

!c check whether this is a .calib file. If so then treat it as a
!c spectrum calibration file, otherwise it is treated as a standard
!c detector calibration file.

                extension=get_ext(fnamed)
                if(extension.ne.'CALIB'.and.extension.ne.'calib') then
 
!c if this number is not equal to ndet then the detector.dat may not be reliable

                    if(ndetd.ne.ndet) then
                        write(6,*) 'ERROR! No. of detectors in calibration file' &
                    ,' is not equal to number in data file',ndetd,ndet
                        stop
                    endif

!c if phiutno is zero, or greater than number of User Table column (nuse)
!c then set phi values to zero

                    do id=1,ndetd
                        if(phiutno.le.0.or.phiutno.gt.nutno) then
                            read(10,*) it,a,b,c,d
                            e=0.0
                        else
                            read(10,*) it,a,b,c,d,(e,j=1,phiutno)
                        endif

!c check that the detector code read in corresponds to that in the header block

                        iflag=0
                        if(it.eq.udet(id)) then
                            deltat(id)=a
                            lendet(id)=b
                            ttheta(id)=d
                            phi(id)=e
                            iflag=1
                        else

!c otherwise find the correct detector number and overwrite the calibration for that
!c detector

                            i=1
                            do while(iflag.eq.0.and.i.le.ndet)
                                if(it.eq.udet(i)) then
                                    deltat(i)=a
                                    lendet(i)=b
                                    ttheta(i)=d
                                    phi(i)=e
                                    iflag=1
                                endif
                                i=i+1
                            end do
                        endif
                        if(iflag.ne.1) write(6,*) &
     'WARNING! Detector ',it,' not found in the header block'
                    end do
                else

!c if number of entries is not equal to nspec, then stop, since the calibration file
!c may not be reliable.

                    if(ndetd.ne.nspec) then
                        write(6,*) 'ERROR: No. of spectra in calibration file is not' &
     ,' equal to number in data file',ndetd,nspec
                        stop
                    endif

!c if phiutno is zero then set phi values to zero

                    do is=1,ndetd
                        if(phiutno.le.0.or.phiutno.gt.nutno) then
                            read(10,*) it,a,b,c,d
                            e=0.0
                        else
                            read(10,*) it,a,b,c,d,(e,j=1,phiutno)
                        endif

!c remember to add 1 to spectrum number since Gudrun spectrum numbers start at 1 not 0 (6/12/2010 Not sure about this
!c anymore - I think it has been changed.)

!c           it=it+1
                        if(is.gt.nspec) then
                            write(6,*) 'Spectrum ',is,' does not exist in raw file'
                            stop
                        endif

!c step through all the detectors and if a detector spectrum number
!c matches the current spectrum value, then update the calibration for this
!c detector

                        id=0
                        test=.false.
                        do while (id.lt.ndet.and..not.test)
                            id=id+1
                            test=is.eq.specno(id)
                        end do
                        if(test) then
                            deltat(id)=a
                            lendet(id)=b
                            ttheta(id)=d
                            phi(id)=e
                        endif
                    end do
                endif
                close(10)
            endif
        endif

        open(10,file='gudrun_calib.dat',status='unknown')
        WRITE(10,*)' No. of detectors = ', ndet
        write(10,*) ' Incident flight path = ',LENIN
        write(10,*) ' Detector flight paths:- '
        write(10,'(5(1x,e13.6))') (lendet(i),i=1,ndet)
        write(10,*) ' Detector scattering angles:- '
        write(10,'(5(1x,e13.6))') (ttheta(i),i=1,ndet)
        write(10,*) ' Detector azimuthal scattering angles'
        write(10,'(5(1x,e13.6))') (phi(i),i=1,ndet)
        close (10)
        return
        
    end subroutine get_calibration

END MODULE calibration_routines
