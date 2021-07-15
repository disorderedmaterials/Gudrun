!     
! File:   normalisation_routines.f90
! Author: aks45
!
! Created on 15 November 2013, 16:30
!

MODULE normalisation_routines
    
    implicit none
    
    CONTAINS
    
    subroutine vdcs_read()

    use reallocation_routines
    use run_par
    use groups_routines
    use van_par
    use merge_routines
    use local_data
    use write_routines

        integer data_unit,ierr
      parameter(data_unit=20)
      character*256 sometext,fname
      character*256 string
        character*32 formatin
      character*1 hash
        
        integer :: mbraggv,iq,ig,nblock,ngrpf,iv
        real    :: pi,vlevel,bavsq,qstep,qmax,qstep2,e
        
        pi=4.0*atan(1.0)
!c ideal level of differential c/s for calibration
      vlevel=vscatav/(4.0*pi)
!c set the minimum level for the interference differential c/s
      bavsq=0.0
      do iv=1,nvelement
            bavsq=bavsq+vfrac(iv)*vscatlen(iv)*vscatlen(iv)
      end do
!c open file
      open(data_unit,file=vdcsname,status='old',iostat=ierr)
        if(ierr.ne.0) then
            vdcsname=convertfilename(vdcsname)
            open(data_unit,file=vdcsname,status='old',iostat=ierr)
        endif
!      write(6,*) 'vdcs_read> ',vdcsname(1:len_trim(vdcsname))
        if(ierr.eq.0) then
!c If this file contains Bragg peak positions, heights and widths, then
! read them now. The corresponding dcs values will be defined once the q-scale is defined
            if(index(vdcsname,'.bragg').gt.0) then
!c Get the Bragg data
!c Ignore lines with # in them
                sometext='#'
                do while(index(sometext,'#').gt.0)
                    read(data_unit,*,iostat=ierr) sometext
                end do
!c Backspace 1 line
                backspace(data_unit)
!c Now read the data. It is assumed that the last line containing data is terminated
!c with a new line character. Otherwise the last line will be lost.
                ierr=0
                nbraggv=0
                mbraggv=0
!c Read the overall factor
              read(data_unit,*,iostat=ierr) braggf
                do while(ierr.eq.0)
                    nbraggv=nbraggv+1
                    do while (nbraggv.gt.mbraggv)
                        mbraggv=mbraggv+10
                        call reallocate1d_r(braggq,mbraggv)
                        call reallocate1d_r(braggh,mbraggv)
                        call reallocate1d_r(braggw,mbraggv)
                    end do
                    read(data_unit,*,iostat=ierr) braggq(nbraggv),braggh(nbraggv),braggw(nbraggv)
                end do
!c Throw away the last line in case it contains garbage
                nbraggv=nbraggv-1
                write(6,102) nbraggv
102      format(/' vdcs_read> ',i3,' Bragg peaks to be defined')
                return
         else if(index(vdcsname,'*').le.0) then
!c Standard Gudrun mint file
!c Read 2 lines
                read(data_unit,*,iostat=ierr)
                read(data_unit,*,iostat=ierr) 
                read(data_unit,'(1x,a80)',iostat=ierr) string
                read(string,*) ngrpf,nqv
                call reallocate1d_r(vqdcs,nqv)
                call reallocate2d_r(vdcs,nqv,ngrpf)
                call check_and_allocate_local_data(nqv)
!c
!c ngrpf is the number of groups in the file
!c If this is different from the number of groups in the current groups file, it
!c is assumed an error has occurred and further analysis will stop
!c
!c Now read lines, ignoring all those beginning with a #

                read(data_unit,*,iostat=ierr) string
                do while(index(string,'#').gt.0.and.ierr.eq.0)
                    read(data_unit,*,iostat=ierr) string
                end do
!c Last line read would not contain the # so backup 1 line
                backspace(data_unit)
!c
!c setup input format specification
!c
                nblock=2*ngrpf+1
                write(formatin,332) '(a1,',nblock,'(1x,e14.7))'
332   format(a4,i3.3,a11)
!c
!c now get data
!c
                iq=0
                do while(ierr.eq.0)
                    iq=iq+1
                    read(data_unit,formatin,iostat=ierr) hash,vqdcs(iq),(vdcs(iq,ig),e,ig=1,ngrpf)
!c Save a copy of the vdcs.
                    detcount(iq)=vdcs(iq,1)
!c The dcs is not allowed to go below its official lower limit
!c
!c                   if(detcount(iq).lt.-bavsq) detcount(iq)=-bavsq
                end do
                nqv=iq-1
                write(6,200) ngrpf,nqv,vdcsname(1:len_trim(vdcsname))
200      format(1x,'vdcs_read> Read ',i3,' groups and ',i5,' Q values from ',a)
                close(data_unit)
!c
!c It is assumed the data read in are interference differential c/s 
!c and oscillate about 0.
!c
!c For the programme we need to divide by expected self scattering level to make result dimensionless
                do iq=1,nqv
                    errcount(iq)=0.0
                    do ig=1,ngrpf
                        vdcs(iq,ig)=vdcs(iq,ig)/vlevel
                    end do
                end do
!c
!c write the result in a file suitable for checking
!c
                fname='vanadium.soq'
                call w_diag_file(fname,nqv,vqdcs,detcount,errcount)
                return
            endif
        else
!c
!c error condition occurred - print a message and stop
!c
            write(6,*) 'Problems with specified VDCS file ',vdcsname
            stop
        endif

    end subroutine vdcs_read

    subroutine vdcs_bragg()

    use reallocation_routines
    use run_par
    use groups_routines
    use van_par
    use merge_routines
    use local_data
    use write_routines

    character*256 fname,string
    integer i,ic,iq,mqv,istart
    real pi,vlevel,qstep,qstep2,qmax,qbin,q1,q2,sum,rat,qdif,qdif2,anorm,e,sqrt2pi

    if(index(vdcsname,'.bragg').gt.0) then
       pi=4.0*atan(1.0)
       sqrt2pi=sqrt(2.0*pi)
!c ideal level of differential c/s for calibration
       vlevel=vscatav/(4.0*pi)
!c The vanadium DCS is calculated in the Q range specified in the input
!c If outputunitstype is 1, then we use the Q-scale specified on input, otherwise we use a 
!c Q scale 0 - 20 in steps of 0.01
       if(outputunitstype.eq.1) then
!Check that qbound is not zero, otherwise start from later value
           istart=1
           do while(istart.lt.nq.and.qbound(istart).le.0.0)
              istart=istart+1
           end do
           nqv=nq-istart+1
           call reallocate1d_r(vqdcs,nqv)
           do ic=1,nqv
              vqdcs(ic)=qbound(ic+istart-1)
           end do
       else
           mqv=100
           qstep=0.01
           qstep2=0.5*qstep
           qmax=20.0
           q1=qstep2
           call reallocate1d_r(vqdcs,mqv)
           vqdcs(1)=q1
           nqv=1
           do while(q1.lt.qmax)
               nqv=nqv+1
               do while (nqv.gt.mqv) 
                   mqv=mqv+100
                   call reallocate1d_r(vqdcs,mqv)
               end do
               ic=nqv
               q2=(2*ic-1)*qstep2
               vqdcs(nqv)=q2
               q1=q2
           end do
       endif
       write(6,*) 'vdcs_bragg> ',nqv,' Q values to be generated in vdcs_bragg'
       call check_and_allocate_local_data(nqv)
!c Step through the bins and generate Gaussians, bearing in mind that these
!c are bin boundaries.
       q1=vqdcs(1)
       do iq=1,nqv-1
          q2=vqdcs(iq+1)
          qbin=0.5*(q1+q2)
!c Step through Bragg data and generate a Gaussian for each peak
          sum=0.0
!          write(6,*) 'vdcs_bragg> ',iq,qbin,nbraggv,braggq(nbraggv),braggh(nbraggv),braggw(nbraggv)
          do i=1,nbraggv
!c Apply Porod's law for Q=0 peak
             if(braggq(i).eq.0.0.and.qbin.gt.0.0) then
                rat=braggw(i)/qbin
                if(rat.gt.0.001) then
                   sum = sum+braggh(i)*rat*rat*rat*rat
                endif
             else
                qdif=(qbin-braggq(i))/braggw(i)
                anorm=1.0/(braggw(i)*sqrt2pi)
                qdif2=0.5*qdif*qdif
                if(qdif2.lt.20.0) sum=sum+braggh(i)*anorm*exp(-qdif2)
             endif
          end do
!         write(6,*) 'vdcs_bragg> ',iq,qbin,sum,braggf
          detcount(iq)=sum*braggf
          q1=q2
       end do
       detcount(nqv)=0.0
       call reallocate2d_r(vdcs,nqv,1)
!c
!c It is assumed the data read in are interference differential c/s 
!c and oscillate about 0.
!c
!c For the programme we need to divide by expected level to make result dimensionless
       do iq=1,nqv
          errcount=0.0
          vdcs(iq,1)=detcount(iq)/vlevel
       end do
!c
!c write the result in a file suitable for checking
!c
       fname='vanadium.soq'
       call w_diag_file(fname,nqv,vqdcs,detcount,errcount)
       return
    else if(index(vdcsname,'*').gt.0) then
!c Set up the vanadium dcs in the event that no data for this have been read in
!c If outputunitstype is 1, then we use the Q-scale specified on input, otherwise we use a 
!c Q scale 0 - 20 in steps of 0.01
       if(outputunitstype.eq.1) then
          nqv=nq
          call reallocate1d_r(vqdcs,nqv)
          call reallocate2d_r(vdcs,nqv,1)
          do iq=1,nq
             vqdcs(iq)=qbound(iq)
             vdcs(iq,1)=0.0
          end do
       else
          qstep=0.01
          qstep2=0.5*qstep
          qmax=20.0
          q1=0
          nqv=1
          mqv=100
          call reallocate1d_r(vqdcs,mqv)
          call reallocate2d_r(vdcs,mqv,1)
          vqdcs(nqv)=q1
          vdcs(nqv,1)=0.0
          do while(q1.lt.qmax)
             q2=(2*nqv-1)*qstep2
             nqv=nqv+1
             do while (nqv.gt.mqv)
                mqv=mqv+100
                call reallocate1d_r(vqdcs,mqv)
                call reallocate2d_r(vdcs,mqv,1)
             end do
             vqdcs(nqv)=q2
             vdcs(nqv,1)=0.0
!          detcount(nqv)=0.0
             q1=q2
          end do
       endif
       call check_and_allocate_local_data(nqv)
       do iq=1,nqv
          detcount(iq)=vdcs(iq,1)
          errcount(iq)=0.0
       end do
!c
!c write the result in a file suitable for checking
!c
       fname='vanadium.soq'
       call w_diag_file(fname,nqv,vqdcs,detcount,errcount)
       return
    endif
    end subroutine vdcs_bragg
    
    subroutine get_nsmoov(smoovlimit)

    use reallocation_routines
    use run_par
    use groups_routines
    use van_par
    use merge_routines
    use local_data
    use write_routines

!c
!c Gets a list of factors for specified runs. Run numbers which are not specified are automatically
!c given a factor of 1.0. If the file runfactor_list.dat does not exist, then all the runs 
!c automatically will have a factor of 1.0.
!c
        use spec_van
        
        integer is,nsmoovalue,ierr,ispec
        real smoovlimit
!c Initialise the array in case no file exists
        nsmoovalue=0
        if(smoovlimit.gt.0.0) nsmoovalue=2
        call reallocate1d_i(nsmoov,nspec)
        do is=1,nspec
            nsmoov(is)=nsmoovalue
        end do
!c Only revise these values if smoovlimit > 0
        if(smoovlimit.gt.0.0) then
            open(10,file='nsmoov.dat',status='unknown',iostat=ierr)
            do while (ierr.eq.0)
                read(10,*,iostat=ierr) ispec,nsmoovalue
                if(ierr.eq.0.and.ispec.gt.0.and.ispec.le.nspec) nsmoov(ispec)=nsmoovalue
            end do
            close(10)
        endif
        return

    end subroutine get_nsmoov

    subroutine save_nsmoov(smoovlimit)

        use reallocation_routines
        use run_par
        use groups_routines
        use van_par
        use merge_routines
        use local_data
        use write_routines
        use spec_van

        real smoovlimit
        integer is,ierr
!c Only save the values if smoovlimit > 0
        if(smoovlimit.gt.0.0) then
            open(10,file='nsmoov.dat',status='unknown',iostat=ierr)
            is=0
            do while (is.lt.nspec)
                is=is+1
                write(10,*,iostat=ierr) is,nsmoov(is)
            end do
            close(10)
        endif
        return
        
    end subroutine save_nsmoov
    
    subroutine dcs_read(fnamedcs,ilenfnamedcs)

        use sam_par
        use merge_routines
        use local_data
        use corrections_routines
        use reallocation_routines
        use run_par
        use groups_routines
        use van_par
        use write_routines
        
        integer data_unit
        integer ilenfnamedcs
        integer nfirst,nlast,iq,iqref,j,ierr,nblock,ngrpf,nsmoodcs

        real pi,qstep,qmax,qstep2,ratio

        character*256 fnamedcs
        character*256 string
        character*32 formatin
        character*1 hash

        parameter(data_unit=20)
        
        pi=4.0*atan(1.0)
!c open file
        ilenfnamedcs=len_trim(fnamedcs)
        open(data_unit,file=fnamedcs,status='old',iostat=ierr)
        if(ierr.ne.0) then
!c error condition occurred - print a message and stop
!c
            write(6,*) 'dcs_read> Problems with specified DCS file ',fnamedcs(1:ilenfnamedcs)
            stop

        else
!c Standard Gudrun dcs file
!c Read 2 lines
            read(data_unit,*,iostat=ierr)
            read(data_unit,*,iostat=ierr) 
            read(data_unit,'(1x,a80)',iostat=ierr) string
            read(string,*) ngrpf,nxdcs
!c
!c ngrpf is the number of groups in the file
!c If this is greater than the expected number, it
!c is assumed an error has occurred and further analysis will stop
!c
            if(ngrpf.gt.ngroup) then
                write(6,104) fnamedcs(1:ilenfnamedcs),ngroup
104             format(/'dcs_read> Number of groups in dcs file ',a,' > ',i5)
                stop
            endif
            call reallocate1d_r(sxdcs,nxdcs)
            call reallocate2d_r(sdcs,nxdcs,ngrpf)
            call reallocate2d_r(esdcs,nxdcs,ngrpf)
!c
!c Now read lines, ignoring all those beginning with a #
            read(data_unit,*,iostat=ierr) string
            do while(index(string,'#').gt.0.and.ierr.eq.0)
                read(data_unit,*,iostat=ierr) string
            end do
!c Last line read would not contain the # so backup 1 line
            backspace(data_unit)
!c
!c setup input format specification
!c
            nblock=2*ngrpf+1
            write(formatin,332) '(a1,',nblock,'(1x,e14.7))'
332   format(a4,i3.3,a11)
!c
!c now get data
!c
            iq=0
            do while(iq.lt.nxdcs.and.ierr.eq.0)
                iq=iq+1
                read(data_unit,formatin,iostat=ierr) hash,sxdcs(iq),(sdcs(iq,j),esdcs(iq,j),j=1,ngrpf)
            end do
            nxdcs=iq-1
            ngrpdcs=ngrpf
            write(6,200) nxdcs,ngrpdcs,fnamedcs(1:ilenfnamedcs)
200         format(/'dcs_read> Read ',i5,' values and ',i5,' groups from ',a)
            !Convert xboundaries to xbins
!            do iq=1,nxdcs-1
!                sxdcs(iq)=0.5*(sxdcs(iq)+sxdcs(iq+1))
!            end do
!            nxdcs=nxdcs-1
            ! Perform a simple smoothing on data read in to remove noise in these data
!            call reallocate1d_r(detcount,nxdcs)
!            call reallocate1d_r(dettemp,nxdcs)
!            call reallocate1d_r(errcount,nxdcs)
!            do j=1,ngrpf
!                do iq=1,nxdcs-1
!                    detcount(iq)=sdcs(iq,j)
!                    errcount(iq)=esdcs(iq,j)**2
!                end do
!                detcount(nxdcs)=0.0
!                errcount(nxdcs)=0.0
!                !Trim the zeros off these values
!                call trim_zeros(1,nxdcs-1,nfirst,nlast,detcount)
!                nsmoodcs=100
!                call smooe(nfirst,nlast,nxdcs,nsmoodcs,detcount,dettemp,errcount,ratio,1.0)
!                write(6,*) 'dcs_read> First and last ',nfirst,nlast,' with ',nsmoodcs,' smoothings on dcs group ',j,' ratio ',ratio
!                do iq=1,nxdcs
!                    sdcs(iq,j)=dettemp(iq)
!                end do
!            end do
        endif
        close(data_unit)
        return

    end subroutine dcs_read
    
!***********************************************************************************
!*
!*      form_smoo_van_b.FOR
!*
!*      A K Soper, January 2003
!*
!*      forms the smoothed ratio detector/monitor for vanadium, using a square or
!*       other broadening function.
!*
!***********************************************************************************
    subroutine form_smoo_van(smoolimit,nspecwrt,alimit,nchanout)

        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines
        use monitor_routines
        use local_data
        use spec_van
        use van_par
        use groups_routines
        use corrections_routines
        use interpolation_routines
        use write_routines
        use bad_detectors

!c
!c internal variables
!c
        character*256 fname            !name of file to write data to
        integer i,ic,ic1,id,is,ib      !internal indices
        integer j,jf,jl,jref            !internal indices
        character*256 run                  !file to be smoothed
        integer naccept,nchanout            !1 if spectrum fit is successful, else 0
        integer nsmoo,nsmoomax            !no. of iterations of vanadium fitting
        integer nfirsts,nlasts      !range of non-zero values
        integer nfirst,nlast            !range of values to fit.
        integer ispec                  !spectrum counter
        integer nspecwrt            !spectrum number to write diagnostic file
        integer nbragg,nl1,nl2
        integer nup,ndown,nbroad
        integer npos,nneg,irefmin
        real lentot                  !total flight path for detector
        real rat,err,v2,b2            !temporary real values
        real wave,wfac            !temporary wavelengths
        real pi      
        real lbragg1,lbragg2            !range of wavelengths for Bragg peak
        real wave1,wave2,val1,val2,grad
        real alimit,chisqratio,smoolimit,vmin

        integer ierr                  !Signals whether an attempt to smooth the vanadium was successful or not

        pi=4.0*atan(1.0)

!c
!c set up the time channel boundaries for the vanadium
!c
        nchan=nchanv
        nchanb=nchanbv
        call reallocate1d_r(tcb,nchanb)
        do ic=1,nchanb
            tcb(ic)=tcbv(ic)
        end do
        call reallocate1d_r(detcount,nchanb)
        call reallocate1d_r(errcount,nchanb)
        call reallocate1d_r(wavebound,nchanb)
        call reallocate1d_r(corrtemp,nchanb)
!
!c Set the maximum number of vanadium smoothings

        nsmoomax=0
        if(smoolimit.gt.0.0) nsmoomax=nchan
        nsmoo=1
!c
!c run number for vanadium
!c
        fname=runv(1)
!c
!c for each detector....
!c
        nspecprocgood=0
        call reallocate1d_i(specprocgood,nspecproc)! Array specprocgood can never be larger than nspecproc
        do is=1,nspecproc
!c
!c Gudrun spectrum number
!c
            ispec=specproc(is)
            nfirst=1
            nlast=nchan
         
            ndettested=ndettested+1
!c
!c temporary store of normalised vanadium data
!c
            do ic=1,nchan
                detcount(ic)=smovandet(ic,is)
                errcount(ic)=errvandet(ic,is)
            end do

            nsmoo=nsmoov(ispec)
            if(inst(1:len_trim(inst)).eq.'D4C'.and.nsmoo.eq.0) then
                !Use simply the relative detector efficiencies to normalise the data
                do ic=1,nchan
                    corrtemp(ic)=effd4c(ispec)
                end do
                if(effd4c(ispec).le.0) then
                    ierr=-1
                else
                    ierr=0
                end if
                chisqratio=0.0
            else
                call smoov_limit(nsmoo,nsmoomax,nchan,detcount,errcount &
                ,corrtemp,chisqratio,smoolimit,nfirst,nlast,alimit,ierr,nchanout)
            end if

!c If the smooth was unsuccessful (went below alimit) signal this as a bad detector 
!c (unless instrument is D4C)

!            if(inst(1:len_trim(inst)).ne.'D4C') then
            if(ierr.ne.0.or.ibad(ispec).ne.0) then
                if(ierr.eq.1) then
                    write(nchanout,103) igrp(ispec),ispec,nsmoo,chisqratio
103                 format(1x,'form_smoo_van> Group: ',i5, ', spectrum: ',i5 &
                    ,' rejected, with ',i4,' smooths and chisq ratio ',f10.5)
                    ibad(ispec) = -1
                else if(ierr.eq.2) then
                    write(nchanout,105) igrp(ispec),ispec,alimit,nsmoo,chisqratio
105                 format(1x,'form_smoo_van> Group: ',i5, ', spectrum: ',i5 &
                    ,' rejected, below limit: ',e10.3,' with ',i4 &
                    ,' smooths and chisq ratio ',f10.5)
                    ibad(ispec) = -2
                else
                    write(nchanout,106) igrp(ispec),ispec
106                 format(1x,'form_smoo_van> Group: ',i5, ', spectrum: ',i5 &
                    ,' rejected before smooth')
                end if
                ndetrejected=ndetrejected+1
            else
                write(nchanout,104) igrp(ispec),ispec,nsmoo,chisqratio
                write(6,104) igrp(ispec),ispec,nsmoo,chisqratio
104             format(1x,'form_smoo_van> Group: ',i5, ' spectrum: ',i5 &
                ' with ',i4,' smooths and chisq ratio ',f10.5)
                nspecprocgood=nspecprocgood+1
                specprocgood(nspecprocgood)=ispec
!If the smooth was successful, we set the smoothed error bar to zero
                if(nsmoo.gt.nchan) then
                    do ic=1,nchan
                        errcount(ic)=0.0
                    end do
                end if
!c
!c save it
!c
                do ic=1,nchan
                    smovandet(ic,nspecprocgood)=corrtemp(ic)
                    errvandet(ic,nspecprocgood)=errcount(ic)
                end do
                smovandet(nchanb,nspecprocgood)=0.0
                errvandet(nchanb,nspecprocgood)=0.0
                corrtemp(nchanb)=0.0
            endif
            nsmoov(ispec)=nsmoo
!            endif
            if(ispec.eq.nspecwrt) then
                call change_ext(fname,'smovan')
!c
!c get the spectrum number of the last spectrum in the last group
!c
                id=detno(ispec)
                lentot=lenin+lendet(id)
                wfac=0.0039554/lentot
                do ic=1,nchanb
                    wavebound(ic)=wfac*(tcb(ic)+deltat(id))
                end do
                call w_diag_file(fname,nchanb,wavebound,corrtemp,errcount)
            endif
        end do
        return
        
    END subroutine form_smoo_van
    
    
!***********************************************************************************
!*
!*      rebin_van_b.FOR
!*
!*      A K Soper, Nov. 2013
!*
!*      Rebins the vanadium onto the same wavelength scale as the corresponding sample detector
!*      Also checks that the vanadium spectrum numbers are the same as the sample. 
!*      (They can change as a result of the vanadium smooth.) Any unwanted spectra in smovandet
!*      are eliminated.!!!!
!*
!***********************************************************************************
    subroutine rebin_van(nspecwrt)

        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines
        use monitor_routines
        use local_data
        use spec_van
        use van_par
        use groups_routines
        use interpolation_routines
        use write_routines
!c
!c internal variables
!c
        character*256 fname            !name of file to write data to
        integer i,ic,ic1,id,is,ib      !internal indices
        integer j,jf,jl,jref            !internal indices
        integer ispec                  !spectrum counter
        integer nspecwrt            !spectrum number to write diagnostic file
        real lentot                  !total flight path for detector
        real rat,err,v2,b2            !temporary real values
        real wave,wfac            !temporary wavelengths
        real pi
        real, dimension(:,:), allocatable       :: smovantemp,errvantemp!Temporary stores

        logical ierr                  !Signals whether an attempt to smooth the vanadium was successful or not

        pi=4.0*atan(1.0)
        
        nchan=nchanv
        nchanb=nchanbv
        call reallocate1d_r(tcb,nchanb)
        do ic=1,nchanb
            tcb(ic)=tcbv(ic)
        end do
!First put the existing data into a temporary store
        call reallocate2d_r(smovantemp,nchanbv,nspecproc)
        call reallocate2d_r(errvantemp,nchanbv,nspecproc)
        smovantemp=smovandet
        errvantemp=errvandet
        call reallocate2d_r(smovandet,nchanbs,nspecproc)
        call reallocate2d_r(errvandet,nchanbs,nspecproc)
!Allocate temporary arrays to the relevant size
        call reallocate1d_r(wavebound,nchanbv)
        call reallocate1d_r(detcount,nchanbv)
        call reallocate1d_r(errcount,nchanbv)
        call reallocate1d_r(wavetemp,nchanbs)
        call reallocate1d_r(dettemp,nchanbs)
        call reallocate1d_r(errtemp,nchanbs)
!c
!c run number for vanadium
!c
        fname=runv(1)
!c
!c for each detector....
!c
        do is=1,nspecproc
!c
!c Gudrun spectrum number
!c
            ispec=specproc(is)
!c
!c set up the vanadium wavelength bins
!c
            id=detno(ispec)
            lentot=lenin+lendet(id)
            wfac=0.0039554/lentot
            do ic=1,nchanb
                wavebound(ic)=wfac*(tcb(ic)+deltat(id))
            end do
!c
!c temporary store of normalised vanadium data
!c
            do ic=1,nchanv
                detcount(ic)=smovantemp(ic,is)
                errcount(ic)=errvantemp(ic,is)
            end do
!c
!c set up the wavelength values for the sample and containers
!c
            do ic=1,nchanbs
                wavetemp(ic)=wfac*(tcbs(ic)+deltat(id))
            end do
!c
!c rebin the vanadium data onto the same wavelength scale as the
!c sample
!c
            call rebinq(nchanv,wavebound,detcount,nchans,wavetemp,dettemp,2,0)
            call rebinq(nchanv,wavebound,errcount,nchans,wavetemp,errtemp,3,0)
            if(ispec.eq.nspecwrt) then
                call change_ext(fname,'smovanbin')
                call w_diag_file(fname,nchanbs,wavetemp,dettemp,errtemp)
            endif
!c
!c save it
!c
            do ic=1,nchans
                smovandet(ic,is)=dettemp(ic)
                errvandet(ic,is)=errtemp(ic)
            end do
            smovandet(nchanbs,is)=0.0
            errvandet(nchanbs,is)=0.0

        end do
        deallocate(smovantemp,errvantemp)
        return
        
    END subroutine rebin_van        

END MODULE normalisation_routines

