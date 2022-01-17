!     
! File:   groups.f90
! Author: aks45
!
! Created on 27 April 2012, 12:59
!

MODULE groups_routines
    
!c
!c arrays to define the spectra grouping
!c
    implicit none
    
    character(len=256)                      :: fname_grps        !name of groups file
    integer                                 :: ngroup,ngroup_in        !number of groups of spectra 
    integer, dimension(:), allocatable      :: ndetgr    !no. of spectra in each group
    integer, dimension(:,:), allocatable    :: indgr    !spectrum number for each group
    integer, dimension(:), allocatable      :: igrp        !group number for each spectrum
    integer                                 :: normalisationtype   !0 = nothing, 1 = <f>^2, 2 = <f^2>
    integer                                 :: ilenfname_grps
    integer, dimension(:), allocatable      :: setbackground_factor !signals whether a group background factor has been set
    integer                                 :: isoftgroupedges     !Determines if soft group edges are to be used
    real, dimension(:), allocatable         :: lengrp        !average secondary flight path for group
    real, dimension(:), allocatable         :: tthgrp        !average scattering angle for group
    real, dimension(:), allocatable         :: phigrp       !average azimuthal angle for group
    real                                    :: rmsfac,stdfac    !number of rms deviations to allow and range of standard devs.
    real, dimension(:), allocatable         :: qmingrp,qmaxgrp    !Min and max Q value for each group
    real, dimension(:), allocatable         :: qmingrp_in,qmaxgrp_in    !Min and max Q value for each group
    real, dimension(:), allocatable         :: tthmingrp,tthmaxgrp !Min and max scattering angle for each group.
    real, dimension(:), allocatable         :: background_factor
    real, dimension(:), allocatable         :: background_factor_in
    real, dimension(:), allocatable         :: amutshield
    real                                    :: qmax_for_groups      !Default qmax for groups
    logical                                 :: inagroup           !Determines if a spectrum belongs to a group
    logical                                 :: softgroupedges
    integer, dimension(2)                   :: olddims,newdims
    integer                                 :: oldsize,newsize
    integer, dimension(:), allocatable      :: igroup_in
    integer                                 :: fname_start !start position to extract calibration filename from

    CONTAINS
    
    subroutine get_groups(ignore,wavemin,wavemax,outputunitstype)
!c
!c routine to read the bad spectra from file SPEC.BAD
!c
!c number of spikes from SPIKE.DAT
!c
!c and current spectra grouping from file GROUPS_DEF.DAT
!c
!c get_groups_test.for - this will also read the detector calibration parameters
!c
        use run_par
        use bad_detectors
        use calibration_routines
        use reallocation_routines
        use beam_routines
    
        integer, intent(in)                     :: ignore            !=1 to ignore existing spec.bad, else 0
        real, intent(in)                        :: wavemin,wavemax
        integer, intent(in)                     :: outputunitstype
        integer, dimension(:), allocatable      :: indgrt      !temporary store for spectrum numbers
        integer, dimension(:), allocatable      :: iread            !dummy array
        integer                                 :: nspike            !number of spectra in spike.dat
        integer                                 :: nwrt,nwrti,nsumd,jcount,ipair
        integer                                 :: i,j,jref,jref1,jref2,jref3      !temporary values
        integer                                 :: is,ispec,nextra,ierr,ibadread,igroup,igroup1            !temporary values
        real                                    :: pi,piconv
        real                                    :: qmintest,qmaxtest,anglemin,anglemax,anglestep
        character(len=256)                      :: fname1      !filename for default groups
        logical                                 :: found
        
        pi=4.0*atan(1.0)
!c To convert 2-theta in degrees to theta in radians
        piconv=pi/360.0
!c
!c set spectrum arrays to zero - 
!c
!c a spectrum is assumed bad unless it appears as zero in the spec.bad file
!c
!c the dummy array iread attempts to check that each spectrum number only 
!c occurs once in the groups table.
!c
        if(allocated(ibad)) deallocate(ibad) !This is in case D4C read routines have been used
        allocate(iread(nspec),igrp(nspec),ibad(nspec),spike(nspec))
        do i=1,nspec
            iread(i)=0
            igrp(i)=0
!c           ibad(i)=-1
            ibad(i)=0
            spike(i)=0
        end do
!c
!c read existing bad detector file if it exists
!c
        if(ignore.eq.0) then
            open(10,file='spec.bad',status='old',iostat=ierr)
            nbread=0
            do while(ierr.eq.0)
                ispec=0
                read(10,*,iostat=ierr) ispec,ibadread
                if(ispec.gt.0.and.ispec.le.nspec) ibad(ispec)=ibadread
                if(ispec.gt.0) nbread=nbread+1
            end do
            close(10)
        else
            do i=1,nspec
                ibad(i)=0
            end do
            nbread=nspec
        endif
        if(nbread.ne.nspec) write(6,74) nbread,nspec
74      format(1x,i5,' spectra from spec.bad DIFFERENT from ',i5,' spectra expected')
        if(ignore.eq.0) then
            open(10,file='spike.dat',status='unknown',iostat=ierr)
            if(ierr.eq.0) read(10,*,iostat=ierr) 
            if(ierr.eq.0) read(10,*,iostat=ierr) 
            if(ierr.eq.0) read(10,*,iostat=ierr) 
            if(ierr.eq.0) read(10,*,iostat=ierr) 
            nspike=0
            ierr=0
            do while (ierr.eq.0)
                read(10,*,iostat=ierr) is,ibadread
                if(ierr.eq.0) then
                    nspike=nspike+1
                    if(is.gt.0.and.is.le.nspec) spike(is)=ibadread
                end if
            end do
            close(10)
        else
            do i=1,nspec
                spike(i)=0
            end do
            nspike=nspec
        endif
        if(nspike.ne.nspec) write(6,73)      nspike,nspec
73      format(1x,i5,' spectra from spike.dat DIFFERENT from ' &
  ,i5,' spectra expected')
!c
!c read the default groups file
!c
        if(index(fname_grps,'*').gt.0) then
            fname1='groups_def.dat'
        else
            fname1=fname_grps
        endif
        if(inst(1:index(inst,' ')-1).eq.'D4C'.and.index(fname_grps,'*').gt.0) then
!For D4C the default groups will simply be to put all valid detectors into a single group
            ngroup=1
            newdims(1)=nspec
            newdims(2)=ngroup
            call reallocate1d_i(ndetgr,newdims(2))
            call reallocate2d_i(indgr,newdims(1),newdims(2))
            ndetgr(1)=0
            do i=1,nspec
                if(ibad(i).eq.0) then
                    ndetgr(1)=ndetgr(1)+1
                    indgr(ndetgr(1),1)=i
                    igrp(i)=1
                else
                    igrp(i)=0
                end if
            end do
        else if (index(fname1,'anglegroups').gt.0) then
!Use detector angles to define groups. In this case we simply use first, last and step in scattering angle to determine how detectors are divided amongst groups.
            write(6,604) fname1(1:len_trim(fname1))
            open(10,file=fname1,status='old',iostat=ierr)
            if(ierr.ne.0) then
                fname1=convertfilename(fname1)
                open(10,file=fname1,status='old',iostat=ierr)
            endif
            if(ierr.eq.0) then
               read(10,*) anglemin,anglemax,anglestep
               close(10)
! Assign initial value of number of groups
               ngroup = max(1,int((anglemax-anglemin)/anglestep))
               newdims(2)=ngroup
               call reallocate1d_i(ndetgr,newdims(2))
               ndetgr=0
               newdims(1)=10
               call reallocate2d_i(indgr,newdims(1),newdims(2))
               olddims(1)=newdims(1)
               olddims(2)=newdims(2)
! Step through the detectors and assign each to a group.
               do i=1,nspec
! get the detector number of this spectrum
                  jref=detno(i)
                  if(ttheta(jref).gt.anglemin.and.ttheta(jref).lt.anglemax.and.ibad(i).eq.0) then
                     igroup=int((ttheta(jref)-anglemin)/anglestep)+1
                     ndetgr(igroup)=ndetgr(igroup)+1
                     if(ndetgr(igroup).gt.olddims(1)) then
                        newdims(1)=olddims(1)+10
                        call reallocate2d_i(indgr,newdims(1),newdims(2))
                        olddims(1)=newdims(1)
                     end if
                     indgr(ndetgr(igroup),igroup)=i
                  end if
               end do
               do i=1,ngroup
                  write(6,*) i,ndetgr(i)!,(indgr(j,i),j=1,ndetgr(i))
               end do
! Now remove any groups with zero detectors
               i=1
               do while (i.le.ngroup)
                  if(ndetgr(i).gt.0) then
                     i=i+1
                  else
                     igroup=i
                     do while (igroup.lt.ngroup)
                        igroup1=igroup+1
                        ndetgr(igroup)=ndetgr(igroup1)
                        do j=1,ndetgr(igroup)
                           indgr(j,igroup)=indgr(j,igroup1)
                        end do
                        igroup=igroup1
                     end do
                     ngroup=ngroup-1
                  end if
               end do
               do i=1,ngroup
                  do j=1,ndetgr(i)
                     igrp(indgr(j,i)) = i
                  end do
                  write(6,*) i,ndetgr(i)!,(indgr(j,i),j=1,ndetgr(i))
               end do
            end if
        else
            write(6,604) fname1(1:len_trim(fname1))
604         format(1x,'Name of groups file being read from:',1x,a)
            open(10,file=fname1,status='old',iostat=ierr)
            if(ierr.ne.0) then
                fname1=convertfilename(fname1)
                open(10,file=fname1,status='old',iostat=ierr)
            endif
            if(ierr.eq.0) then
                read(10,*) ngroup
                write(6,*) ngroup
                !allocate the space for ndetgrp and indgr
                newdims(1)=10
                newdims(2)=ngroup
                call reallocate1d_i(ndetgr,newdims(2))
                call reallocate1d_i(indgrt,newdims(1))
                call reallocate2d_i(indgr,newdims(1),newdims(2))
                olddims(1)=newdims(1)
                olddims(2)=newdims(2)
                nsumd=0
                do i=1,ngroup
                    read(10,*) nwrt
                    nwrti=iabs(nwrt)
                    write(6,*) nwrt,nwrti
                    backspace(10)
!c
!c read the first line
!c
                    if(nwrti.le.10) then
                        read(10,*) nwrt,(indgrt(j),j=1,nwrti)
                    else
!c
!c if nwrti.gt.10 the first line is different to others
!
!c number of extra lines beyond the first
!
                        nextra=(nwrti-11)/10+1
!c
!c read the first
!c
                        read(10,*) nwrt,(indgrt(j),j=1,10)
!c
!c read the remainder
!c
                        jref1=11
                        do is=1,nextra
                            jref2=jref1+9
                            if(jref2.gt.nwrti) jref2=nwrti
                            do while (jref2.gt.olddims(1))
                                newdims(1)=olddims(1)+10
                                call reallocate1d_i(indgrt,newdims(1))
                                call reallocate2d_i(indgr,newdims(1),newdims(2))
                                olddims(1)=newdims(1)
                            end do
                            read(10,*) (indgrt(j),j=jref1,jref2)
                            jref1=jref2+1
                        end do
                    endif
                    if(nwrt.lt.0) then
                        ndetgr(i)=nwrti
                        nsumd=nsumd+ndetgr(i)
                        do j=1,ndetgr(i)
                            do while (j.gt.olddims(1))
                                newdims(1)=olddims(1)+10
                                call reallocate1d_i(indgrt,newdims(1))
                                call reallocate2d_i(indgr,newdims(1),newdims(2))
                                olddims(1)=newdims(1)
                            end do
                            indgr(j,i)=indgrt(j)
                            jref=indgr(j,i)
                            igrp(jref)=i
                            if(iread(jref).ne.0) then
                                write(6,609) jref
609   format(1x,'Spectrum ',i4,' found more than once in groups file')
                            end if
                            iread(jref)=iread(jref)+1
                        end do
                    else
                        jcount=0
                        do ipair=1,nwrt,2
                            jref2=2*ipair
                            jref1=jref2-1
                            do jref=indgrt(jref1),indgrt(jref2)
                                jcount=jcount+1
                                do while (jcount.gt.olddims(1))
                                    newdims(1)=olddims(1)+10
                                    call reallocate1d_i(indgrt,newdims(1))
                                    call reallocate2d_i(indgr,newdims(1),newdims(2))
                                    olddims(1)=newdims(1)
                                end do
                                jref3=jref
                                indgr(jcount,i)=jref3
                                igrp(jref3)=i
                                if(iread(jref3).ne.0) then
                                    write(6,609) jref3
                                endif
                                iread(jref3)=iread(jref3)+1
                            end do
                        end do
                        ndetgr(i)=jcount
                        nsumd=nsumd+ndetgr(i)
                    endif
                end do
                close(10)
                write(6,607) nsumd,ngroup,fname1(1:len_trim(fname1))
            else
                write(6,200) fname1(1:len_trim(fname1))
200             format(/'get_groups> Specified groups file: ',a,' does not exist')
                stop
            end if
        end if
607     format(1x,i5,' spectra in ',i2,' groups from file ',a)
!c
!c calculate the average secondary flight path, scattering angle, and azimuthal
!c angle for each group
!c
!c first zero the accumulators
!
        call reallocate1d_r(lengrp,ngroup)
        call reallocate1d_r(tthgrp,ngroup)
        call reallocate1d_r(phigrp,ngroup)
        call reallocate1d_i(iread,ngroup)
        call reallocate1d_r(tthmingrp,ngroup)
        call reallocate1d_r(tthmaxgrp,ngroup)
        call reallocate1d_r(qmingrp,ngroup)
        call reallocate1d_r(qmaxgrp,ngroup)
        call reallocate1d_r(background_factor,ngroup)
        call reallocate1d_i(setbackground_factor,ngroup)
        do j=1,ngroup
            lengrp(j)=0.0
            tthgrp(j)=0.0
            phigrp(j)=0.0
            iread(j)=0
            tthmingrp(j)=180.0
            tthmaxgrp(j)=0.0
            qmingrp(j)=1000.0
            qmaxgrp(j)=0.0
            background_factor(j)=bakfac
        end do
!c
!c step through the spectra and get the appropriate detector parameters if igrp
!c is non-zero.
!c
        do i=1,nspec
            j=igrp(i)
            if(j.gt.0.and.ibad(i).eq.0) then
!
!c get the detector number of this spectrum
!c
                jref=detno(i)
                iread(j)=iread(j)+1
                lengrp(j)=lengrp(j)+lendet(jref)
                tthgrp(j)=tthgrp(j)+ttheta(jref)
                phigrp(j)=phigrp(j)+phi(jref)
! Save the maximum and minimum for this group
                if(ttheta(jref).lt.tthmingrp(j)) tthmingrp(j)=ttheta(jref)
                if(ttheta(jref).gt.tthmaxgrp(j)) tthmaxgrp(j)=ttheta(jref)
            else
!c
!c if a spectrum has not been assigned to a group it is assumed to be bad
!c
!c                  ibad(i)=-1
            endif
        end do
!        write(6,*) 'get_groups> Done sum values', ngroup_in,size(qmingrp_in),size(qmaxgrp_in),size(background_factor_in)
!c
!c Set q-ranges of groups, and normalize total values to the number of values
!c
        do j=1,ngroup
            !Check through any values read in from the input file
            i=0
            found=.false.
            do while(i.lt.ngroup_in.and..not.found)
                i=i+1
                found=igroup_in(i).eq.j
            end do
!            write(6,*) j,i,found,qmingrp_in(i)
!               write(6,*) qmaxgrp_in(i)
!               write(6,*) background_factor_in(i)
            if(found) then
                qmingrp(j)=qmingrp_in(i)
                qmaxgrp(j)=qmaxgrp_in(i)
                background_factor(j)=background_factor_in(i)
                setbackground_factor(j)=1
            else    
                !First set the q-ranges and background factors to default values
                qmingrp(j)=0.0
                qmaxgrp(j)=qmax_for_groups
                background_factor(j)=bakfac
                setbackground_factor(j)=1
            end if
            write(6,*) 'get_groups> Checked input values for group ', j
            if(iread(j).gt.0) then
                lengrp(j)=lengrp(j)/float(iread(j))
                tthgrp(j)=tthgrp(j)/float(iread(j))
                phigrp(j)=phigrp(j)/float(iread(j))
!c The Q range for each group will correspond to the largest possible Q value at the smallest angle and
!c smallest possible Q value at the largest angle. This will not work for reactor data.
                if(index(inst,'D4C').eq.0.and..not.softgroupedges) then
                     if(outputunitstype.eq.1) then
                         qmintest=4.0*pi*sin(tthmaxgrp(j)*piconv)/wavemax
                         qmaxtest=4.0*pi*sin(tthmingrp(j)*piconv)/wavemin
                         if(qmintest.gt.qmingrp(j)) qmingrp(j)=qmintest
                         if(qmaxtest.lt.qmaxgrp(j)) qmaxgrp(j)=qmaxtest
                     else if(outputunitstype.eq.2) then            
                          qmaxtest=wavemax/(2.0*sin(tthmaxgrp(j)*piconv))
                          qmintest=wavemin/(2.0*sin(tthmingrp(j)*piconv))
                          if(qmintest.gt.qmingrp(j)) qmingrp(j)=qmintest
                          if(qmaxtest.lt.qmaxgrp(j)) qmaxgrp(j)=qmaxtest
                     endif
                endif
            end if
        end do

!c
!c write out the specified Q-range for each group
!c
        write(6,*)
        write(6,*) 'Q-ranges for individual groups'
        write(6,*)
        do i=1,ngroup
            write(6,*) i,qmingrp(i),qmaxgrp(i)
        end do

        open(10,file='gudrun_grp.dat',status='unknown')
        write(10,*) ' Group flight paths:- '
        write(10,'(5(1x,e13.6))') (lengrp(i),i=1,ngroup)
        write(10,*) ' Group scattering angles:- '
        write(10,'(5(1x,e13.6))') (tthgrp(i),i=1,ngroup)
        write(10,*) ' Group azimuthal angles:- '
        write(10,'(5(1x,e13.6))') (phigrp(i),i=1,ngroup)
        write(10,*) ' Group number of each spectrum:-'
        do i=1,nspec,5
            jref=min(nspec,i+4)
            write(10,'(5(1x,i5,1x,i2,1x))') (j,igrp(j),j=i,jref)
        end do
        close (10)
!c Write a file to store group calibration suitable to be read by Gudrun_GUI
            fname1=fname1(1:len_trim(fname1))//'.cal'
#ifdef GUDPY_COMPATIBILITY
            fname_start=INDEX(fname1, "/", BACK=.TRUE.)
            IF (fname_start == 0 ) THEN
                fname_start=INDEX(fname1, "\", BACK=.TRUE.)
            END IF
            IF (fname_start > 0) THEN
                fname_start=fname_start+1
            END IF
            fname1=fname1(fname_start:len_trim(fname1))
#endif
        write(6,*) fname1
        open(10,file=fname1,status='unknown')
        do i=1,ngroup
            write(10,'(i5,1x,3(1x,e13.6))') i,lengrp(i),tthgrp(i),phigrp(i)
        end do
        close(10)
        return

    end subroutine get_groups
    
    subroutine write_groups(runs)

        use run_par
        use bad_detectors
    
!c
!c routine to write the bad spectra to file SPEC.BAD and 'INST'nnnnn.BAD
!c
!c the number of spikes to file SPIKE.DAT
!c
!c and current spectra grouping to 'INST'nnnnn.GRP
!c
        character(len=256), intent(in)      :: runs            !current run number
        integer                             :: i,is,j,idet,ndets
        integer                             :: k,k1,nwrt,nsumd    !temporary values
        character(len=256)                  :: fname        !temporary filename
!c
!c write the bad spectrum file
!c
        write(fname,502) runs(1:index(runs,'.')-1)
502     format(a,'.bad')
        open(10,file='spec.bad',status='unknown')
        open(11,file=fname,status='unknown')
        do i=1,nspec
            write(10,111) i,ibad(i)
            write(11,111) i,ibad(i)
111     format(1x,i5,1x,i8)
        end do
        close(10)
!c
!c write the spike.dat file
!c
        open(10,file='spike.dat',status='unknown')
        write(10,*) 'spike.dat'
        write(10,*)
        write(10,*)
        write(10,*)
        do is=1,nspec
            write(10,*) is,spike(is),0.0
        end do
        close(10)
!c
!c write the new groups file.
!c
        write(fname,504) runs(1:index(runs,'.')-1)
504     format(a,'.grp')
!c
!c First step through each group and remove bad spectra
!c
        nsumd=0
        do i=1,ngroup
            ndets=ndetgr(i)
            j=1
            do while (j.le.ndetgr(i))
                idet=indgr(j,i)
!c no. of detectors in a group cannot be less than 1
                if((ibad(idet).ne.0.or.spike(idet).ne.0).and.ndetgr(i).gt.1) then
                    ndetgr(i)=ndetgr(i)-1
                    k=j
                    do while(k.le.ndetgr(i))
                        k1=k+1
                        indgr(k,i)=indgr(k1,i)
                        k=k+1
                    end do
                else
                    j=j+1
                endif
            end do
            nsumd=nsumd+ndetgr(i)
        end do
!c
!c now write the new groups file
!c
        open(10,file=fname,status='unknown')
            write(10,302) ngroup
302     format(1x,i5)
            do i=1,ngroup
                nwrt=-ndetgr(i)
                write(10,303) nwrt,((indgr(j,i)),j=1,ndetgr(i))
303     format(1x,i6,((10(1x,i5))))
            end do
        close(10)
        write(6,608) nsumd,ngroup,fname
608     format(1x,i5,' spectra in ',i2,' groups to file ',a50)
        return
    end subroutine write_groups

END MODULE groups_routines
