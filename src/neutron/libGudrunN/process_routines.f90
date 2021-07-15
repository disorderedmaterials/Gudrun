!     
! File:   process_routines.f90
! Author: aks45
!
! Created on 20 November 2013, 16:43
!

MODULE process_routines
      
    !Routines which process the data detector by detector
    
    implicit none
    
    CONTAINS
    
!***********************************************************************************
!*
!*    divide_by_mon.FOR
!*
!*    A K Soper, November 1999
!*
!*    forms the normalised ratio detector/smoothed monitor
!*
!***********************************************************************************
    subroutine divide_by_mon(ntype,nspecwrt)

        USE run_par
        use inputfilestrings
        use reallocation_routines
        USE calibration_routines
        USE groups_routines
        USE beam_routines
        use local_data
        use monitor_routines
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use write_routines
        use interpolation_routines
        use summed_data_routines
!c
!c internal variables
!c
    character(len=256) fname    !name of file to write data to
    integer i,ic,ic1,id,is,imfind    !internal indices
    integer j,jf,jl,jref,ind    !internal indices
    integer naccept        !1 if spectrum fit is successful, else 0
    integer nsmoo            !no. of smoothings on detector counts
    integer nfirst,nlast        !range of non-zero values
    integer ntype            !type of data 1=v,2=bv,3=s,4=bs
    integer nspecwrt        !spectrum no. for diagnostic file
    integer isref            !Gudrun spectrum number
        integer mchan
        integer nspecmon                !Used to get the correct monitor spectrum for D4C data
    real, dimension(:), allocatable         :: wavemon,smomon,smoscount!(mchan)        !monitor wavelength boundaries
    real lentot            !total flight path for detector
    real rat,err            !temporary real values
    real wave1,wave2        !temporary wavelengths
    real wavemin,wavemax        !minimum and maximum wavelengths to use
    real monsum            !temporary monitor sum

!c
!c Set up the time channel boundaries
!c
        if(ntype.lt.3) then
            nchan=nchanv
            nchanb=nchanbv
            call reallocate1d_r(tcb,nchanb)
            do ic=1,nchanb
                  tcb(ic)=tcbv(ic)
            end do
        else
            nchan=nchans
            nchanb=nchanbs
            call reallocate1d_r(tcb,nchanb)
            do ic=1,nchanb
                  tcb(ic)=tcbs(ic)
            end do
        endif
        mchan=nchanb !For compatibility with fixed array versions
        !Allocate the required arrays
        call reallocate1d_r(wavemon,nchanb)
        call reallocate1d_r(smomon,nchanb)
        call reallocate1d_r(smoscount,nchanb)
        call reallocate1d_r(wavebound,nchanb)
        call reallocate1d_r(detcount,nchanb)
        call reallocate1d_r(errcount,nchanb)

!c run number for vanadium and monitor values for time-of-flight experiment

    if(ntype.eq.1) then
            fname=runv(1)
            do ic=1,nchanb
        wavemon(ic)=waveimon(ic)
        smomon(ic)=smovanmon(ic)
            end do
    else if(ntype.eq.2) then
            fname=runbv(1)
            do ic=1,nchanb
                wavemon(ic)=waveimon(ic)
        smomon(ic)=smobvanmon(ic)
            end do
    else if(ntype.eq.3) then
            fname=runbs(1)
            do ic=1,nchanb
        wavemon(ic)=waveimon(ic)
        smomon(ic)=smobsammon(ic)
            end do
    else 
            ind=ntype-3
            write(6,*) 'divide_by_mon ',nchan,nchanb,ind,size(smosammon,1),size(smosammon,2)
            fname=runs(1,ind)
            do ic=1,nchanb
        wavemon(ic)=waveimon(ic)
        smomon(ic)=smosammon(ic,ind)
            end do
    end if

!c
!c integrate the monitor between specified wavelength range
!c
    monsum=0.0
    imfind=0
    do ic=1,nchan
       if(wavemon(ic).ge.monwave1.and.wavemon(ic+1).le.monwave2) then
          imfind=imfind+1
          monsum=monsum+smomon(ic)
       endif
    end do
    if(imfind.gt.0) then
       monsum=monsum/imfind
!c
!c set the monitor count
!c
       do ic=1,nchanb
          smoscount(ic)=monsum
       end do
    endif
!c
!c for each spectrum....
!c
    do is=1,nspecproc

        isref=specproc(is)
!c
!c detector number for this spectrum
!c
        id=detno(isref)
!c
!c total flight path
!c
        lentot=lenin+lendet(id)
        rat=0.0039554/lentot

!c
!c set up wavelength boundaries and width of each channel
!c
        do ic=1,nchanb
                    wavebound(ic)=rat*(tcb(ic)+deltat(id))
        end do
!c
!c get counts for this spectrum
!c
        call get_spec(ntype,is,detcount,errcount)
                !For D4C the relevant monitor is stored in the particular module block associated with this spectrum
                if(index(inst,'D4C').gt.0) then
                    nspecmon=crate(detno(isref))
                    call get_summed_data_spec(ntype,nspecmon,smomon)
                    if(isref.eq.nspecwrt) write(6,*) 'divide_by_mon> ',isref,nspecmon,smomon(1),detcount(1),errcount(1)
                end if

!c
!c if monitor wavelength range is zero then
!c
        if(imfind.le.0) then
!c
!c rebin monitors onto same wavelength as detector
!c
!           if(isref.eq.nspecwrt) then
!              call rebinq(nchanb,wavemon,smomon,nchanb,wavebound,smoscount,1,1)
!           else
              call rebinq(nchan,wavemon,smomon,nchan,wavebound,smoscount,1,0)
                      if(isref.eq.nspecwrt) write(6,*) 'divide_by_mon> ',smomon(1),smoscount(1)
!           endif
        endif
!c
!c divide detector by monitor
!c
        nfirst=0
        nlast=0
        do ic=1,nchanb
                    rat=smoscount(ic)
                    if(rat.gt.0.0) then
            if(nfirst.eq.0) nfirst=ic
            nlast=ic
            detcount(ic)=detcount(ic)/rat
            errcount(ic)=errcount(ic)/rat/rat
                    else
            detcount(ic)=0.0
            errcount(ic)=0.0
                    endif
        end do
!c
!c save results
!c        
        call put_spec(ntype,is,detcount,errcount)
!c
!c write diagnostic file if required
!c
                if(isref.eq.nspecwrt) then
 !                   write(6,*) 'divide_by_mon> ',fname(1:len_trim(fname))
                    call change_ext(fname,newext='normmon')
 !                   write(6,*) 'divide_by_mon> ',fname(1:len_trim(fname))
                    call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
        endif
    end do
    return    
    END subroutine divide_by_mon
    
!***********************************************************************************
!*
!*      subtract_bak.FOR
!*
!*      A K Soper, November 1999
!*
!*      subtracts background from vanadium or sample
!*
!***********************************************************************************
    subroutine subtract_bak(ntype,wavemin,wavemax,nspecwrt)

        USE run_par
        use inputfilestrings
        use reallocation_routines
        USE calibration_routines
        USE groups_routines
        USE beam_routines
        use local_data
        use monitor_routines
        use van_par
        use sam_par
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use bad_detectors
        use write_routines
        use interpolation_routines
!c
!c internal variables
!c
        character*256 fname            !name of file to write data to
        integer i,ic,ic1,id,is      !internal indices
        integer j,jf,jl,jref,ind      !internal indices
        integer ntype                  !type of data 1=v,2=bv,3=bs,4=s,5=c1,6=c2,...
        character*256 run                  !filename of data
        integer nspecwrt            !spectrum no. for diagnostic file
        integer isref                  !Gudrun spectrum number
        integer igrpref                  !Group number of specified spectrum
        integer nwav,nwavtrans,nfirst,nlast
        real lentot,wave,wfac            !total flight path for detector
        real rat,err                  !temporary real values
        real bakfact,prod,prod1                  !temporary background factors
        real wav1,wav2
        real wavemin,wavemax,backaccept
        real, dimension(:), allocatable :: templa,temptrans

        logical transdependentbackground   !decides whether sample dependent background is to be subtracted

!Set a factor to determine whether data is acceptable relative to the background
        backaccept=vanbackfraction
        ind=ntype-3
!c
!c Set up the time channel boundaries
!c
        if(ntype.lt.3) then
            nchan=nchanv
            nchanb=nchanbv
            nwav=nlambv
            nwavtrans=nlambtransv
            call reallocate1d_r(tcb,nchanb)
            do ic=1,nchanb
                tcb(ic)=tcbv(ic)
            end do
        else
            nchan=nchans
            nchanb=nchanbs
            nwav=nlambs
            nwavtrans=nlambtranss
            call reallocate1d_r(tcb,nchanb)
            do ic=1,nchanb
                tcb(ic)=tcbs(ic)
            end do
        endif
        !Allocate the required arrays
        call reallocate1d_r(wavebound,nchanb)
        call reallocate1d_r(wavebin,nchanb)
        call reallocate1d_r(detcount,nchanb)
        call reallocate1d_r(errcount,nchanb)
        call reallocate1d_r(corrbin,nchanb)
        call reallocate1d_r(wavecorr,nwav)
        call reallocate1d_r(corrtemp,nwav)
        call reallocate1d_r(templa,nwavtrans)
        call reallocate1d_r(temptrans,nwavtrans)
!
!c get the wavelength boundary values which are needed for the sample transmission values.
!c
        if(ntype.eq.1) then
            fname=runv(1)
            do ic=1,nwavtrans
                templa(ic)=vlambtrans(ic)
                temptrans(ic)=vantransmission(ic)
            end do
            do ic=1,nwav
                wavecorr(ic)=vlamb(ic)
                corrtemp(ic)=ainter(templa,temptrans,wavecorr(ic),nwavtrans) !Interpolate measured transmission onto wavelength scale for correctiosn.
            end do
        else
            ind=ntype-3
            fname=runs(1,ind)
            do ic=1,nwavtrans
                templa(ic)=slambtrans(ic)
                temptrans(ic)=samtransmission(ic,ind)
            end do
            do ic=1,nwav
                wavecorr(ic)=slamb(ic)
                corrtemp(ic)=ainter(templa,temptrans,wavecorr(ic),nwavtrans)
            end do
        endif
!c
!c for each spectrum
!c
        do is=1,nspecproc
!c
!c get the group number and background factor for this spectrum
!c
            isref=specproc(is)
            igrpref=igrp(isref)
!            write(6,*) 'subtract_bak> ',is,isref,igrpref
!c
!c background factor for this group.
!c
            bakfact=background_factor(igrpref)
            if(bakfact.lt.0.0) then

                  transdependentbackground = .true.

            else

                  transdependentbackground = .false.

            endif
            bakfact=abs(bakfact)

            id=detno(isref)
            lentot=lenin+lendet(id)
            wfac=0.0039554/lentot
            do ic=1,nchanb
                wavebound(ic)=wfac*(tcb(ic)+deltat(id))
            end do
            wav1=wavebound(1)
            do ic=1,nchan
                wav2=wavebound(ic+1)
                wavebin(ic)=0.5*(wav1+wav2)
                wav1=wav2
            end do

!c Now interpolate the current sample transmission onto the wavelength scale for this spectrum

!            write(6,*) 'subtract_bak> 1 ',is,isref,id,lentot,wavebin(1)
            call inter(nchan,wavebin,corrbin,nwav,wavecorr,corrtemp)
!            write(6,*) 'subtract_bak> 2 ',is,isref,id,nchan,size(wavebin),size(corrbin),nwav,size(wavecorr),size(corrtemp)

!c Write this to a diagnostic file

            if(isref.eq.nspecwrt) then
                call change_ext(fname,'trans')
                call w_diag_file(fname,nchan,wavebin,corrbin,errcount)
            endif
!c
!c get the data to be background corrected
!c
            if(ntype.eq.1) then
!c
!c run number
!c
                run=runv(1)
!c
!c subtract background - errors are not saved
!c
                do ic=1,nchan

                    !Select the first wavelength bin so that the upper limit is within specifed wavelength range
                    !and lower limit is within specified wavlength range
                    if(wavebound(ic+1).gt.wavemin.and.wavebound(ic).lt.wavemax) then 
                        prod=bakfact
                        prod1=0.0
                        if(transdependentbackground) then
                            prod1=1.0-prod
                            prod = prod*corrbin(ic)
                        endif
                        prod=prod+prod1
                        smovandet(ic,is)=smovandet(ic,is)-prod*normbakdet(ic,is)
                        errvandet(ic,is)=errvandet(ic,is)+prod*prod*errbakdet(ic,is)
                    else
                        smovandet(ic,is)=0.0
                        errvandet(ic,is)=0.0
                    endif
                end do
                !Now check that background subtracted data is at least backaccept times bigger than background
                !First find the first wavelength for which the background acceptance criterion is satisfied
                nfirst=0
                nlast=0
                do ic=1,nchan
                    if(smovandet(ic,is).gt.backaccept*normbakdet(ic,is)) then
                        nlast=ic
                        if(nfirst.eq.0) nfirst=ic
                    end if
                end do
                if(nfirst.gt.0.and.nlast.gt.0) then
                    !Now set it
                    do ic=1,nchan
                        if(ic.lt.nfirst.or.ic.gt.nlast) then
                            smovandet(ic,is)=0.0
                            errvandet(ic,is)=0.0
                        end if
                    end do
                else
                    !Otherwise signal a bad detector
                    ibad(isref)=-1
                end if

            else if(ind.gt.0) then
!c
!c subtract background and calculate std.s 
!c
                do ic=1,nchan
                    if(wavebound(ic+1).gt.wavemin.and.wavebound(ic).lt.wavemax) then
                        prod=bakfact
                        prod1=0.0
                        if(transdependentbackground) then
                           prod1=1.0-prod
                           prod = prod*corrbin(ic)
                        endif
                        prod=prod+prod1
                        normsamdet(ic,is,ind)=normsamdet(ic,is,ind)-prod*normbakdet(ic,is)
                        errsamdet(ic,is,ind)=errsamdet(ic,is,ind)+prod*prod*errbakdet(ic,is)
                    else
                        normsamdet(ic,is,ind)=0.0
                        errsamdet(ic,is,ind)=0.0
                    endif
                end do
            endif
!c
!c write out diagnostic file if needed
!c
            if(isref.eq.nspecwrt) then
                call change_ext(fname,'subbak')
                if(ntype.eq.1) then
                    do ic=1,nchanb
                        detcount(ic)=smovandet(ic,is)
                        errcount(ic)=errvandet(ic,is)
                    end do
                else if(ind.gt.0) then
                    do ic=1,nchanb
                        detcount(ic)=normsamdet(ic,is,ind)
                        errcount(ic)=errsamdet(ic,is,ind)
                    end do
                endif
                call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
            endif
        end do
        return
        
    end subroutine subtract_bak 
    
!***********************************************************************************
!*
!*      do_mul_corr.for
!*
!*      A K Soper, December 1999
!*
!*      performs the multiple scattering correction. 
!*
!*      For vanadium it first multiplies the single scattering by 1+P, 
!*      adds the multiple scattering, then divides the data by this 
!*      total scattering
!*
!***********************************************************************************
    subroutine do_mul_corr(ntype,nspecwrt)
        
        USE run_par
        use inputfilestrings
        use reallocation_routines
        USE calibration_routines
        USE groups_routines
        USE beam_routines
        use local_data
        use monitor_routines
        use van_par
        use sam_par
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use mul_corr_azi
        use van_placzek
        use van_mul_azi
        use sam_mul_azi
        use corrections_routines
        use bad_detectors
        use write_routines
        use interpolation_routines

!c
!c internal variables
!c
        character*256 fname            !name of file to write data to
        integer i,ic,ic1,ic2,id,is,ib      !internal indices
        integer j,jf,jl,jref,ind,ig,k      !internal indices
        integer ntype                  !vanadium, sample, can, furnace
        integer nfirst,nlast            !indices for range of channels to use
        integer npfirst,nplast            !indices for range of channels to use
        integer nfirst0,nlast0            !indices for range of channels to use
        integer npfirst0,nplast0            !indices for range of channels to use
        integer ncfirst,nclast      !range of non-zero values in data
        character*256 run                  !temporary filename
        integer nwavpla            !no.of wavelengths for Placzek correction
        integer nspecwrt            !spectrum no. for diagnostic file
        integer isref                  !Gudrun spectrum number
        integer ngeom                  !sample geometry for correction angles
        real, dimension(:), allocatable         :: onebin,totbin,plabin,vdcsbin!(mchan)            !interpolated single scattering
        real, dimension(:), allocatable         :: vwavdcs,vdcscor!(mq)
        real lentot                  !total flight path for detector
        real theta,phid             !Scattering direction for a detector
        real theta0,phid0            !scattering direction for a detector at 0,0
        real rat,err,ratp            !temporary real values
        real rat0,ratp0                  !temporary real values
        real cordif1,cordif2,q1,q2
        real onecor1,onecor2            !temporary values
        real mulcor1,mulcor2            !temporary values
        real wave1,wave2            !temporary wavelengths
        real sang                        !sample rotation angle for f.p. corrections
        real pi,pi4,tconv,qconv
        real chansum,chansume             !sum of channels
        real amushield,amufact,amuterm    !Used to correct multiple scattering for attenuation in shielding

!c
!c define pi
!c
        pi=4.0*atan(1.0)
        pi4=4.0*pi
        tconv=pi/360.0
!c
!c run number for sample - used for output
!c
        if(ntype.eq.1) then
            fname=runv(1)
            ngeom=vgeom
            sang=vrot
            ind=1
            nchan=nchanv
            nchanb=nchanbv
            call reallocate1d_r(tcb,nchanb)
            do ic=1,nchanb
                tcb(ic)=tcbv(ic)
            end do
        else if(ntype.gt.2) then
            ind=ntype-3
            fname=runs(1,ind)
            ngeom=sgeom
            sang=srot(1)
            nchan=nchans
            nchanb=nchanbs
            call reallocate1d_r(tcb,nchanb)
            do ic=1,nchanb
                tcb(ic)=tcbs(ic)
            end do
        else
            return
        end if
!c
!c get the corrections in temporary array
!c
        if(ntype.eq.1) then
            nwavmul=nvwavmul
            nangmul=nvangmul
            nazimul=nvazimul
            nwavpla=nvwavpla
            call reallocate1d_r(wavmul,nwavmul)
            call reallocate1d_r(angmul,nangmul)
            call reallocate1d_r(azimul,nazimul)
            call reallocate3d_r(onescat,nwavmul,nangmul,nazimul)
            call reallocate3d_r(mulscat,nwavmul,nangmul,nazimul)
            call reallocate3d_r(totscat,nwavmul,nangmul,nazimul)
            do ic=1,nwavmul
                  wavmul(ic)=vwavmul(ic)
            end do
            do ib=1,nangmul
                angmul(ib)=vangmul(ib)
                do ic=1,nwavmul
                    do k=1,nazimul
                        onescat(ic,ib,k)=vonescat(ic,ib,k)
                        mulscat(ic,ib,k)=vmulscat(ic,ib,k)
!c
!c Use the single scattering to estimate the total scattering by the sample in all directions.
!c This is because the multiple scattering may contain oscillations specific to each particular direction, and
!c so introduce additional structure not present in the data. These m.s. oscillations will be subtracted with
!c the regular multiple scattering
!c
!                     totscat(ic,ib,k)=onescat(ic,ib,k)+mulscat(ic,ib,k)
                        totscat(ic,ib,k)=onescat(ic,ib,k)
                    end do
                end do
            end do
            do ic=1,nazimul
                  azimul(ic)=vazimul(ic)
            end do
            call reallocate1d_r(vwavdcs,nqv)
            call reallocate1d_r(vdcscor,nqv)
            call reallocate1d_r(vdcsbin,nchan)
            call reallocate1d_r(onebin,nchan)
            call reallocate1d_r(totbin,nchan)
            call reallocate1d_r(vancor,nspec)
            call reallocate1d_r(vancore,nspec)
        else
            ind=ntype-3
            nwavmul=nswavmul
            nangmul=nsangmul
            nazimul=nsazimul
            nwavpla=nvwavpla
            call reallocate1d_r(wavmul,nwavmul)
            call reallocate1d_r(angmul,nangmul)
            call reallocate1d_r(azimul,nazimul)
            call reallocate3d_r(onescat,nwavmul,nangmul,nazimul)
            call reallocate3d_r(mulscat,nwavmul,nangmul,nazimul)
            call reallocate3d_r(totscat,nwavmul,nangmul,nazimul)
            do ic=1,nwavmul
                wavmul(ic)=swavmul(ic)
            end do
            do ic=1,nazimul
                azimul(ic)=sazimul(ic)
            end do
            do ib=1,nangmul
                angmul(ib)=sangmul(ib)
                do k=1,nazimul
                    do ic=1,nwavmul
!c Get the total scattering for the sample. This will be used to estimate the sample dependent background
                        onescat(ic,ib,k)=sonescat(ic,ib,k,ind)
                        mulscat(ic,ib,k)=smulscat(ic,ib,k,ind)
                        totscat(ic,ib,k)=onescat(ic,ib,k)+mulscat(ic,ib,k)
                    end do
                end do
            end do
            call reallocate2d_r(samcor,nspec,ind)
            call reallocate2d_r(samcore,nspec,ind)
        endif
        call reallocate1d_r(wavebound,nchanb)
        call reallocate1d_r(wavebin,nchan)
        call reallocate1d_r(detcount,nchanb)
        call reallocate1d_r(errcount,nchanb)
        call reallocate1d_r(corrbin,nchan)
        if(subtractsampledependentbackground) call reallocate1d_r(totbin,nchan)
        call reallocate1d_r(plabin,nchan)
        call reallocate1d_r(corrtemp,nwavpla)
!c
!c Get correction angles for a detector at 0,0
!c
        call get_corr_ang(ngeom,sang,0.0,0.0,theta0,phid0)
!c
!c Get interpolation values the detector at 0,0 (for sample dependent background term)
!c
        call get_interp_values(angmul,nangmul,azimul,nazimul,theta0,phid0,rat0,nfirst0,nlast0,ratp0,npfirst0,nplast0)
!      write(6,*) theta0,phid0,rat0,nfirst0,nlast0,ratp0,npfirst0,nplast0
!c
!c      write(6,*) nangmul,(angmul(ib),ib=1,nangmul)
!c      write(6,*) (ibad(is),is=1,20)
!c      write(6,*) (detno(is),is=1,20)
!c      write(6,*) ndet,(ttheta(id),id=1,20)
!c
!c for each detector....
!c
        do is=1,nspecproc
!c
!c Gudrun spectrum number
!c
            isref=specproc(is)
            if(ibad(isref).eq.0) then

!c Get the shielding attenuation factor for this grp

                amushield=amutshield(igrp(isref))
!c
!c detector number of this spectrum
!c
                id=detno(isref)
!c
!c get the corresponding correction angles for this detector
!c
                call get_corr_ang(ngeom,sang,ttheta(id),phi(id),theta,phid)
!c
!c total flight path
!c
                lentot=lenin+lendet(id)
                rat=0.0039554/lentot
                wave1=rat*(tcb(1)+deltat(id))      
                wavebound(1)=wave1
!c
!c set up wavelength boundaries
!c
                do ic=1,nchan
                    ic1=ic+1
                    wave2=rat*(tcb(ic1)+deltat(id))
                    wavebin(ic)=0.5*(wave1+wave2)
                    wave1=wave2
                    wavebound(ic1)=wave2
                end do
!c
!c get normalised data for this detector
!c
                call get_spec(ntype,is,detcount,errcount)
!c
!c determine the first and last non-zero values of these data
!c
                call trim_zeros(1,nchan,ncfirst,nclast,detcount)
!                write(6,*) 'do_mul_corr> 1',is,isref,nchan,ncfirst,nclast,detcount(nchan) &
!                ,errcount(nchan)
!c
!c First get the interpolation values for the current detector
!c
                call get_interp_values(angmul,nangmul,azimul,nazimul,theta,phid,rat,nfirst,nlast,ratp,npfirst,nplast)
!                write(6,*) 'do_mul_corr> 2 ',id,ttheta(id),phi(id),theta,phid,rat,nfirst,nlast,ratp &
!                ,npfirst,nplast
!c
!c Do interpolation of the multiple scattering term
!c
                call do_interp_corr(nchan,wavebin,corrbin,nwavmul,wavmul,mulscat,rat,nfirst,nlast,ratp,npfirst,nplast)
!                write(6,*) 'do_mul_corr> 3 ',wavebin(1),corrbin(1),size(wavebin),size(corrbin)
!c
!c Also do the interpolation of total scattering at 0,0 if required
!c
                if(subtractsampledependentbackground) then

                    call do_interp_corr(nchan,wavebin,totbin,nwavmul,wavmul,totscat,rat0,nfirst0,nlast0,ratp0,npfirst0,nplast0)
!c
!c Add this to the multiple scattering using the supplied factor
!c
                    do ic=1,nchan
                        if(ic.lt.ncfirst.or.ic.gt.nclast) then
                            corrbin(ic)=0.0
                        else
                            corrbin(ic)=corrbin(ic)+sampledependentbackgroundfactor*totbin(ic)
                        endif
                    end do
                endif

!c Correct multiple scattering for attenuation in the shielding

                if(amushield.gt.0.0) then

                    amufact=amushield*lendet(id)/1.798

                    do ic=ncfirst,nclast
                        amuterm=amufact*wavebin(ic)
                        if(amuterm.lt.30.) then
                            amuterm=exp(-amuterm)
                        else
                            amuterm=0.0
                        endif
                        corrbin(ic)=corrbin(ic)*amuterm
                    end do

                endif
                if(ntype.eq.1) then
!c
!c for vanadium need to do the single atom scattering as well
!c
                    call do_interp_corr(nchan,wavebin,onebin,nwavmul,wavmul,onescat,rat,nfirst,nlast,ratp,npfirst,nplast)
!                    write(6,*) 'do_mul_corr> 4 ',wavebin(1),onebin(1),size(wavebin),size(onebin)
!c
!c now get the Placzek correction for this detector
!c
!c first determine the relevant group number, and get the Placzek correction
!c into a temporary array
!c
                    ig=igrp(isref)
                    do ic=1,nwavpla
                        corrtemp(ic)=vplaczek(ic,ig)
                    end do
                    call inter(nchan,wavebin,plabin,nwavpla,vwavpla,corrtemp)
!                    write(6,*) 'do_mul_corr> 5 ',wavebin(1),plabin(1),size(wavebin),size(plabin)
!c
!c write diagnostic file if spectrum number is correct
!c
                    if(isref.eq.nspecwrt) then
                        call change_ext(fname,'corrtemp')
                        call w_diag_file(fname,nwavpla,vwavpla,corrtemp,errcount)
                        call change_ext(fname,'plabin')
                        call w_diag_file(fname,nchan,wavebin,plabin,errcount)
                    endif
!c
!c convert vanadium dcs q-scale to wavelength for this detector
!c
                    qconv=pi4*sin(tconv*ttheta(id))
                    q1=qconv/vqdcs(nqv)
                    vwavdcs(1)=q1
                    do ic=2,nqv
!c
!c Put the vdcs data into a temporary array, remembering to invert the order of the data, and that q-scale
!c is bin boundaries
!c
                        jref=nqv-ic+1
                        q2=qconv/vqdcs(jref)
                        vwavdcs(ic)=q2
                        vdcscor(ic-1)=vdcs(jref,1)
                        q1=q2

                    end do
                    vdcscor(nqv)=0.0
!                write(6,*) 'do_mul_corr> 5',is,isref,nqv,size(vwavdcs),size(vdcscor)
!c
!c write diagnostic file if spectrum number is correct
!c
                    if(isref.eq.nspecwrt) then
                        call change_ext(fname,'vdcscor')
                        call w_diag_file(fname,nqv-1,vwavdcs,vdcscor,vdcscor)
                    endif
!c
!c interpolate onto same wavelength scale as scattering data
!c
!                write(6,*) 'do_mul_corr> 6',is,isref
                    call rebinq(nqv-1,vwavdcs,vdcscor,nchan,wavebound,vdcsbin,2,0)
!                write(6,*) 'do_mul_corr> 7',is,isref
!c
!c write diagnostic file if spectrum number is correct
!c
                    if(isref.eq.nspecwrt) then
                        call change_ext(fname,'vdcsbin')
                        call w_diag_file(fname,nchan,wavebin,vdcsbin,vdcsbin)
                    endif
!c
!c for vanadium divide data by sum of single and multiple scattering
!c for everthing else subtract the multiple scattering
!c
                    chansum=0.0
                    chansume=0.0
                    do ic=ncfirst,nclast
                        cordif1=onebin(ic)*(1.0+plabin(ic)+vdcsbin(ic))+corrbin(ic)
                        corrbin(ic)=cordif1
                        if(cordif1.gt.0.0) then
                            detcount(ic)=detcount(ic)/cordif1
                            errcount(ic)=errcount(ic)/(cordif1*cordif1)
                            chansum=chansum+detcount(ic)
                            chansume=chansume+errcount(ic)
                        else
                            detcount(ic)=0.0
                            errcount(ic)=0.0
                        endif
                    end do
                    detcount(nchanb)=0.0
                    errcount(nchanb)=0.0
!c
!c write diagnostic file if spectrum number is correct
!c
                    if(isref.eq.nspecwrt) then
                        call change_ext(fname,'vancor')
                        call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
                        call change_ext(fname,'corrbin')
                        call w_diag_file(fname,nchan,wavebin,corrbin,errcount)
                        call change_ext(fname,'onebin')
                        call w_diag_file(fname,nchan,wavebin,onebin,errcount)
                    endif
                    vancor(isref)=chansum/(nclast-ncfirst+1)
                    vancore(isref)=chansume/(nclast-ncfirst+1)
                else
                    chansum=0.0
                    chansume=0.0
                    do ic=ncfirst,nclast
!c
!c only modify values with non-zero error bars
!c
                        detcount(ic)=detcount(ic)-corrbin(ic)
                        chansum=chansum+detcount(ic)
                        chansume=chansume+errcount(ic)
                    end do
                    detcount(nchanb)=0.0
!c
!c write diagnostic file if spectrum number is correct
!c
                    if(isref.eq.nspecwrt) then
                        call change_ext(fname,'mulcor')
                        call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
                        call change_ext(fname,'corrbin')
                        call w_diag_file(fname,nchan,wavebin,corrbin,errcount)
                    endif
                    samcor(isref,ntype-3)=chansum/(nclast-ncfirst+1)
                    samcore(isref,ntype-3)=chansume/(nclast-ncfirst+1)
                endif
            else
                do ic=1,nchanb
                    detcount(ic)=0.0
                    errcount(ic)=0.0
                end do
            end if
!c
!c save the corrected data in the appropriate array
!c
            call put_spec(ntype,is,detcount,errcount)
        end do
        return      
    
    END subroutine do_mul_corr
    
!***********************************************************************************
!*
!*    divide_by_van.FOR
!*
!*    A K Soper, November 1999
!*
!*    forms the normalised ratio detector/vanadium
!*
!***********************************************************************************
    subroutine divide_by_van(ntype,nspecwrt,factornorm)

        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines
        use local_data
        use spec_van
        use spec_sam
        use spec_baks
        use groups_routines
        use write_routines
        use interpolation_routines
!c
!c internal variables
!c
        character*256 fname        !name of file to write data to
        integer i,ic,ic1,id,is,ind    !internal indices
        integer ntype            !type of data 1=v,2=bv,3=s,4=bs
        integer nspecwrt        !spectrum no. for diagnostic file
        integer isref            !Gudrun spectrum number
        real rat,err,wfac,lentot,wave    !temporary real values
        real factornorm             !Factor to multiply the normalisation by prior to dividing into sample.
        real samratio,errratio  !Used to calculate the errors
        ind=ntype-3
!c
!c run number for sample
!c
        if(ind.gt.0) then
            nchan=nchans
            nchanb=nchanbs
            call reallocate1d_r(tcb,nchanb)
            do ic=1,nchanb
                tcb(ic)=tcbs(ic)
            end do
            fname=runs(1,ind)
!c
!c for each detector....
!c
            do is=1,nspecproc
                isref=specproc(is)
!c
!c spectrum number for this detector
!c
                id=detno(isref)
!c
!c calculate counts and std.s per bin width
!c
                do ic=1,nchan
                    rat=smovandet(ic,is)*factornorm
                    if(rat.gt.1.0e-10) then
                        samratio=normsamdet(ic,is,ind)/rat
                        normsamdet(ic,is,ind)=samratio
                        errsamdet(ic,is,ind)=errsamdet(ic,is,ind)/rat/rat
                        if(errvandet(ic,is).gt.0) then
                            errratio=samratio*samratio*errvandet(ic,is)/rat/rat
                            errsamdet(ic,is,ind)=errsamdet(ic,is,ind)+errratio
                        end if
                    else
                        normsamdet(ic,is,ind)=0.0
                        errsamdet(ic,is,ind)=0.0
                    endif
                end do
!c
!c write out this detector if needed for diagnostic purposes
!c
                if(isref.eq.nspecwrt) then
                    id=detno(nspecwrt)
                    call change_ext(fname,'normvan')
                    lentot=lenin+lendet(id)
                    wfac=0.0039554/lentot
                    call reallocate1d_r(wavebound,nchanb)
                    call reallocate1d_r(detcount,nchanb)
                    call reallocate1d_r(errcount,nchanb)
                    do ic=1,nchanb
                        wavebound(ic)=wfac*(tcb(ic)+deltat(id))
                        detcount(ic)=normsamdet(ic,is,ind)
                        errcount(ic)=errsamdet(ic,is,ind)
                    end do
                    call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
                endif
            end do
        end if
    return
        
    END subroutine divide_by_van   

!***********************************************************************************
!*
!*    do_abs_corr.for
!*
!*    A K Soper, December 1999
!*    P Buchanan
!*    performs the absorption correction.
!*
!***********************************************************************************
    subroutine do_abs_corr(nspecwrt)

        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines
        use local_data
        use spec_sam
        use abs_corr_azi
        use groups_routines
        use beam_routines
        use sam_par
        use corrections_routines
        use write_routines
        use interpolation_routines
!c
!c internal variables
!c
    character(len=256) fname        !name of file to write data to
    integer i,ic,ic1,ic2,id,is,iq    !internal indices
    integer j,jf,jl,jref,ind    !internal indices
    integer ind3,ncontac        !internal indices
    integer ntype            !1 = S, 2 = S+C, 3 = S+C+F
    integer nfirst,nlast        !indices for range of channels to use
    integer npfirst,nplast    !indices for range of channels to use
    integer ccont            !container/sample counter
    integer nspecwrt        !spectrum no. for diagnostic file
    integer isref            !Gudrun spectrum number
    real, dimension(:,:), allocatable       :: samcount!(mchan,mcont)    !normalised data for sample
    real, dimension(:,:), allocatable       :: serrcount!(mchan,mcont)    !error array for sample
    real, dimension(:,:), allocatable       :: atten!(mcorrwav,mcont1)    !temporary corrections at current angle
    real, dimension(:,:), allocatable       :: attenbin!(mchan,mcont1)    !temporary corrections at current angle
    real                                    :: lentot            !total flight path for detector
    real rat,ratp,err        !temporary real values
    real angdif,phidif,amutterm
    real cordif1,cordif2
    real atten1,atten2        !temporary values
    real wave1,wave2        !temporary wavelengths
    real samfac            !temporary factors to apply attenuation
        real tweakextra             !additional tweak factor over 1.0
!c
!c run number for sample - used for output
!c
    fname=runs(1,1)
!c
!c number of attenuation corrections is (ncont*(ncont+1))/2
!c
    ncontac=(ncont*(ncont+1))/2
        nchan=nchans
        nchanb=nchanbs
        call reallocate1d_r(tcb,nchanb)
        do ic=1,nchanb
            tcb(ic)=tcbs(ic)
        end do
        call reallocate2d_r(samcount,nchanb,ncont)
        call reallocate2d_r(serrcount,nchanb,ncont)
        call reallocate2d_r(atten,nwavabs,ncontac)
        call reallocate2d_r(attenbin,nchanb,ncontac)
        call reallocate1d_r(wavebound,nchanb)
        call reallocate1d_r(wavebin,nchanb)
        call reallocate1d_r(corrtemp,nwavabs)
        call reallocate1d_r(corrbin,nchanb)
        call reallocate1d_r(detcount,nchanb)
        call reallocate1d_r(errcount,nchanb)
        call reallocate1d_r(dettemp,nchanb)
        call reallocate1d_r(errtemp,nchanb)
!c
!c for each detector....
!c
    do is=1,nspecproc
!c
!c Gudrun spectrum number 
!c
            isref=specproc(is)
!c
!c detector number of this spectrum
!c
            id=detno(isref)
!            write(6,*) 'do_abs_corr> 1 ',is,isref,id,ncont
!c
!c get the corresponding correction angles for this detector
!c
            call get_corr_ang(sgeom,srot(1),ttheta(id),phi(id),theta,phid)
!c
!c total flight path
!c
            lentot=lenin+lendet(id)
            rat=0.0039554/lentot
            wave1=rat*(tcb(1)+deltat(id))    
            wavebound(1)=wave1
!c
!c set up wavelength boundaries
!c
            do ic=1,nchan
        ic1=ic+1
        wave2=rat*(tcb(ic1)+deltat(id))
                wavebound(ic1)=wave2
        wavebin(ic)=0.5*(wave1+wave2)
        wave1=wave2
            end do
            call get_spec(4,is,detcount,errcount)
!c
!c get normalised data for this detector
!c
            do ind=1,ncont
                ind3=ind+3
        call get_spec(ind3,is,detcount,errcount)
        do ic=1,nchan
!c Calculate the sample environment container/sample correction factor if requested
                    if(sampleenvamut(ind).gt.0.0) then
                        amutterm=sampleenvamut(ind)*wavebin(ic)
                        if(amutterm.lt.30.0) then
                            amutterm=sampleenvscatfrac(ind)+(1.0-sampleenvscatfrac(ind))*exp(-amutterm)
                        else
                            amutterm=sampleenvscatfrac(ind)
                        endif
                        samcount(ic,ind)=detcount(ic)/amutterm
                        serrcount(ic,ind)=errcount(ic)/(amutterm*amutterm)
                    else
                        samcount(ic,ind)=detcount(ic)
                        serrcount(ic,ind)=errcount(ic)
                    endif
        end do
            end do
!c
!c determine the nearest pair of correction scattering angles to theta (which
!c are assumed to be in order of increasing angle).
!c
            ic=1
            nfirst=ic
            do while(ic.lt.nangabs.and.theta.gt.angabs(ic))
                nfirst=ic
        ic=ic+1
            end do
            nlast=ic
            angdif=angabs(nlast)-angabs(nfirst)
            if(angdif.gt.0.0) then
                rat=(theta-angabs(nfirst))/angdif
            else
                rat=1.0
            endif
!c
!c determine the nearest pair of phi values, which are assumed to be in order of 
!c increasing angle
!c
            ic=1
            npfirst=1
            do while (ic.lt.naziabs.and.phid.gt.aziabs(ic))
                npfirst=ic
        ic=ic+1
            end do
            nplast=ic
            phidif=aziabs(nplast)-aziabs(npfirst)
            if(phidif.gt.0.0) then
                ratp=(phid-aziabs(npfirst))/phidif
            else
                ratp=1.0
            endif
!c
!c set up (by interpolation or extrapolation) the correction at this angle
!c
            do ic=1,nwavabs
                do ind=1,ncontac
!c
!c interpolate between theta and phi angles
!c
                    cordif1=abscor(ic,nlast,npfirst,ind)-abscor(ic,nfirst,npfirst,ind)
                    cordif2=abscor(ic,nlast,nplast,ind)-abscor(ic,nfirst,nplast,ind)
                    atten1=abscor(ic,nfirst,npfirst,ind)+rat*cordif1
                    atten2=abscor(ic,nfirst,nplast,ind)+rat*cordif2
                    atten(ic,ind)=atten1+ratp*(atten2-atten1)
                end do
            end do
!c
!c interpolate attenuation corrections onto the same wavelength scale as data
!c
            do ind=1,ncontac
                do ic=1,nwavabs
                    corrtemp(ic)=atten(ic,ind)
        end do
        call inter(nchan,wavebin,corrbin,nwavabs,wavabs,corrtemp)
                do ic=1,nchan
                    attenbin(ic,ind)=corrbin(ic)
        end do
            end do
!c
!c Now that the absorption correction coefficients are obtained, they
!c need to be applied to the data
!c
!c Counter to indicate the current can or sample 
!c
            ccont=ncont
!c
!c Index to indicate the correction to be used
!c
            ind=ncontac
            do while (ccont.gt.0) 
!c
!c We start by correcting the outermost container or sample for attenuation
!c
                do ic=1,nchan
                    if(tweak(ccont).ge.0.0.or.ccont.eq.1) then
                        tweakextra=0.0
                        samfac=abs(tweak(ccont))/attenbin(ic,ind)
                    else
!c The next lines were inserted to allow for part of the can scattering to not pass through the sample
                        tweakextra=abs(tweak(ccont))
                        samfac=(1.0-tweakextra)/attenbin(ic,ind)
                    endif
                    samfac=samfac/snorm(ccont)
                    detcount(ic)=samcount(ic,ccont)*samfac
                    errcount(ic)=serrcount(ic,ccont)*(samfac*samfac)
!c Save the uncorrected counts: this is for the case where part of the container scattering does
!c not pass through the inner containers and sample
                    dettemp(ic)=tweakextra*samcount(ic,ccont)
                    errtemp(ic)=tweakextra*tweakextra*serrcount(ic,ccont)
                end do
!c
!c exit the loop if ccont=1
!c
                if(ccont.gt.1) then
!c
!c otherwise subtract the data for sample ccont from the remainder datasets
!c
                    jref=0
                    do j=1,ccont-1
                        jref=jref+ccont-j+1
                        do ic=1,nchan
                            samfac=attenbin(ic,jref)
                            samcount(ic,j)=samcount(ic,j)-detcount(ic)*samfac-dettemp(ic)
                            serrcount(ic,j)=serrcount(ic,j)+errcount(ic)*(samfac*samfac)+errtemp(ic)
                        end do
                    end do
                    ind=jref-1
                endif
                ccont=ccont-1
            end do
!c
!c save the corrected data in the appropriate array
!c
            call put_spec(4,is,detcount,errcount)
!c Also save the corrected data into the array with no subtractions
            do ic=1,nchanb
                normsamdetnosub(ic,is)=detcount(ic)
            end do
!c
!c write the data to a diagnostic file if this is the correct spectrum
!c
            if(isref.eq.nspecwrt) then
                call change_ext(fname,'abscor')
                write(6,*) 'do_abs_corr> 3 ',fname(1:len_trim(fname))
                call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
            endif
!            write(6,*) 'do_abs_corr> 4 ',is,isref,id,nchanb,wavebound(nchanb/2),detcount(nchanb/2),errcount(nchanb/2)
    end do
!        write(6,*) 'do_abs_corr> 5 ',nspecproc,specproc(1),specproc(nspecproc)
    return
        
    END subroutine do_abs_corr
        
!***********************************************************************************
!*
!*    do_Placzek_corr.for
!*
!*    A K Soper, May 2001
!*
!*    performs the Placzek correction. 
!*
!*    it assumes all attenuation corrections have been done.
!*
!*    this routine is not used for the vanadium, where the Placzek correction
!*    is done in do_mul_corr
!*
!***********************************************************************************
    subroutine do_Placzek_corr(nspecwrt)

        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines    
        use groups_routines
        use local_data
        use spec_van
        use spec_sam
        use pla_corr
        use sam_placzek
        use sam_par
        use write_routines
        use interpolation_routines
!c
!c internal variables
!c
    character*256 fname        !name of file to write data to
    integer i,ic,ic1,ic2,id,is,iq    !internal indices
    integer j,jf,jl,jref,ind,ig    !internal indices
    integer ntype            !vanadium, sample, can, furnace
    integer nfirst,nlast        !indices for range of channels to use
    integer nspecwrt        !spectrum no. for diagnostic file
    integer isref            !Gudrun spectrum number
    real lentot            !total flight path for detector
    real theta            !scattering angle for a detector
    real rat,err            !temporary real values
    real angdif,cordif        !temporary values
    real wave1,wave2        !temporary wavelengths
    real pi,pi4,pifac

!c cross section factor for Placzek correction

    pi=4.0*atan(1.0)
    pi4=4.0*pi
    pifac=sscatav(1)/pi4

!c This routine can only be applied to the corrected sample DCS

!c run number for sample - used for output

    fname=runs(1,1)
        nchan=nchans
        nchanb=nchanbs
        call reallocate1d_r(tcb,nchanb)
        do ic=1,nchanb
            tcb(ic)=tcbs(ic)
        end do
        call reallocate1d_r(wavebound,nchanb)
        call reallocate1d_r(wavebin,nchanb)
        call reallocate1d_r(detcount,nchanb)
        call reallocate1d_r(errcount,nchanb)
        call reallocate1d_r(corrtemp,nswavpla)
        call reallocate1d_r(corrbin,nchanb)

!c for each detector....

    do is=1,nspecproc

!c Gudrun spectrum number

            isref=specproc(is)

!c detector number of this spectrum

            id=detno(isref)

!c scattering angle of this detector

            theta=ttheta(id)

!c total flight path

            lentot=lenin+lendet(id)
            rat=0.0039554/lentot
            wave1=rat*(tcb(1)+deltat(id))    
            wavebound(1)=wave1

!c set up wavelength boundaries

            do ic=1,nchan
                ic1=ic+1
        wave2=rat*(tcb(ic1)+deltat(id))
        wavebin(ic)=0.5*(wave1+wave2)
        wave1=wave2
        wavebound(ic1)=wave2
            end do

!c get normalised data for this detector

            call get_spec(4,is,detcount,errcount)

!c Determine range of channels over which to perform correction. We find the first and last non-zero values
!c and use that range

            nfirst=0
            nlast=0
            do ic=1,nchan
                if(errcount(ic).ne.0.0) then
                    nlast=ic
                    if(nfirst.eq.0) nfirst=ic
                endif
            end do

!c Determine the nearest Placzek angle to this detector

            ig=igrp(isref)

!c Interpolate the correction for this angle onto the detector wavelength scale

            do ic=1,nswavpla
                corrtemp(ic)=splaczek(ic,ig)
            end do
            call inter(nchan,wavebin,corrbin,nswavpla,swavpla,corrtemp)

!c subtract the correction

            do ic=1,nchan
                if(ic.ge.nfirst.and.ic.le.nlast) then
                    detcount(ic)=detcount(ic)-corrbin(ic)*pifac
        else
                    detcount(ic)=0.0
                    errcount(ic)=0.0
        endif
            end do

!c save the corrected data in the appropriate array

            call put_spec(4,is,detcount,errcount)

!c write diagnostic file if needed

            if(isref.eq.nspecwrt) then
                call change_ext(fname,'placor')
        call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
            endif
    end do
    return
        
    END subroutine do_Placzek_corr
        
!***********************************************************************************
!*
!*      do_self_corr.for
!*
!*      A K Soper, June 2009
!*
!*      performs the inelasticity correction read in from a file, stored in sdcs. 
!*
!*      it assumes all attenuation corrections have been done.
!*
!************************************************************************************
    subroutine do_self_corr(outputunitstype,nspecwrt)

        use run_par
        use inputfilestrings
        use reallocation_routines
        use calibration_routines
        use groups_routines
        use local_data
        use spec_van
        use spec_sam
        use pla_corr
        use sam_placzek
        use sam_par
        use corrections_routines
        use write_routines
        use interpolation_routines
        
        integer outputunitstype,nspecwrt  !Defines the units type and diagnostic spectrum
!c
!c internal variables
!c
        character*256 fname            !name of file to write data to
        integer i,ic,ic1,ic2,id,is,iq,icref,ncref      !internal indices
        integer j,jf,jl,jref,ind,ig,igold      !internal indices
        integer ntype                  !vanadium, sample, can, furnace
        integer nfirst,nlast,nfirstdata,nlastdata            !indices for range of channels to use
        integer isref,isrefold                  !Gudrun spectrum number
        real lentot                  !total flight path for detector
        real theta                  !scattering angle for a detector
        real rat,err                  !temporary real values
        real angdif,cordif            !temporary values
        real wave1,wave2            !temporary wavelengths
        real pi,pi4,pifac,angconv,stheta

        nchan=nchans
        nchanb=nchanbs
        call reallocate1d_r(tcb,nchanb)
        do ic=1,nchanb
            tcb(ic)=tcbs(ic)
        end do
 
        !c cross section factor for Placzek correction

        pi=4.0*atan(1.0)
        pi4=4.0*pi
        angconv=pi/360.0
        pifac=sscatav(1)/pi4

!c This routine can only be applied to the corrected sample DCS

!c run number for sample - used for output

        fname=runs(1,1)
        call reallocate1d_r(wavecorr,nxdcs)
        call reallocate1d_r(wavetemp,nxdcs)
        call reallocate1d_r(corrtemp,nxdcs)
        call reallocate1d_r(wavebound,nchanb)
        call reallocate1d_r(wavebin,nchanb)
        call reallocate1d_r(detcount,nchanb)
        call reallocate1d_r(errcount,nchanb)

!c Get the inelasticity correction into a temporary array. It is assumed the supplied
!c data are in point format. If the input data are Q values, we have to reverse
!c the order of values in order to get the wavelength values in increasing order.

        isref=specproc(1)
        ig=igrp(isref)
        !Ensure that if the group requested does not exist in the sdcs array, then use the last group
        !This will be used if the merged data have been read in rather than individual groups.
        if(ig.gt.ngrpdcs.or.ig.lt.1) ig=ngrpdcs 
        igold=ig
        nfirst=0
        nlast=0
        do iq=1,nxdcs
            if(outputunitstype.eq.3) then !Assume for wavelength units, the dcs has been read in as a function of Q, so we need to invert the order of Q values
                ic=nxdcs-iq+1
            else
                ic=iq
            endif
            wavecorr(ic)=sxdcs(iq)
            corrtemp(ic)=sdcs(iq,ig)
        end do
        !Need to trim these data to remove any leading or trailing zeros
        call trim_zeros(1,nxdcs,nfirst,nlast,corrtemp)
        if(nfirst.gt.1) then
            ic=0
            icref=nfirst-1
            do while (icref.lt.nlast)
                ic=ic+1
                icref=icref+1
                corrtemp(ic)=corrtemp(icref)
                wavecorr(ic)=wavecorr(icref)
            end do
            ncref=ic
        else
            ncref=nlast
        endif
      
!c for each detector....

        do is=1,nspecproc

!c Gudrun spectrum number

            isref=specproc(is)

!c detector number of this spectrum

            id=detno(isref)

! Get the group number of this detector and generate new self scattering term if different from igold

            ig=igrp(isref)
            if(ig.gt.ngrpdcs.or.ig.lt.1) ig=ngrpdcs !Ensure that if the group requested does not exist in the sdcs array, then use the last group
            if(ig.ne.igold) then !Only get these data if the group number has changed from previously
                igold=ig
                nfirst=0
                nlast=0
                do iq=1,nxdcs
                    if(outputunitstype.eq.3) then
                        ic=nxdcs-iq+1
                    else
                        ic=iq
                    end if
                    wavecorr(ic)=sxdcs(iq)
                    corrtemp(ic)=sdcs(iq,ig)
                end do
                !Need to trim these data to remove any leading or trailing zeros
                call trim_zeros(1,nxdcs,nfirst,nlast,corrtemp)
                if(nfirst.gt.1) then
                    ic=0
                    icref=nfirst-1
                    do while (icref.lt.nlast)
                        ic=ic+1
                        icref=icref+1
                        corrtemp(ic)=corrtemp(icref)
                        wavecorr(ic)=wavecorr(icref)
                    end do
                    ncref=ic
                else
                    ncref=nlast
                end if
            end if

!c scattering angle of this detector

            theta=ttheta(id)
            if(outputunitstype.ne.3) then
                do ic=1,ncref
                    wavetemp(ic)=wavecorr(ic)
                end do
            else
                !For wavelength units we assume input dcs is on a Q scale
                stheta=pi4*sin(theta*angconv)
                do ic=1,ncref
                    wavetemp(ic)=stheta/wavecorr(ic)
                end do
            end if

!c total flight path

            lentot=lenin+lendet(id)
            rat=0.0039554/lentot
            wave1=rat*(tcb(1)+deltat(id))      
            wavebound(1)=wave1

!c set up wavelength boundaries

            do ic=1,nchan
                ic1=ic+1
                wave2=rat*(tcb(ic1)+deltat(id))
                wavebin(ic)=0.5*(wave1+wave2)
                wave1=wave2
                wavebound(ic1)=wave2
            end do

!c get normalised data for this detector

            call get_spec(4,is,detcount,errcount)

!c Determine range of channels over which to perform correction. We find the first and last non-zero values
!c and use that range

            call trim_zeros(1,nchan,nfirstdata,nlastdata,detcount)
            
!c interpolate correction onto the same wavelength scale as data

!            call inter(nchan,wavebin,corrbin,ncref,wavetemp,corrtemp)
            call rebinq(ncref-1,wavetemp,corrtemp,nchan,wavebound,corrbin,2,0)
            call trim_zeros(1,nchan,nfirst,nlast,corrbin)
! Set the minimum and maximum range for the subtraction
            nfirst=max(nfirst,nfirstdata)
            nlast=min(nlast,nlastdata)
!c subtract the correction

            do ic=1,nchan
                if(ic.ge.nfirst.and.ic.le.nlast) then
                    detcount(ic)=detcount(ic)-corrbin(ic)
                else
                    detcount(ic)=0.0
                    errcount(ic)=0.0
                endif
            end do

!c save the corrected data in the appropriate array

            call put_spec(4,is,detcount,errcount)

!c write diagnostic file if needed

            if(isref.eq.nspecwrt) then
                call change_ext(fname,'selfcor')
                call w_diag_file(fname,nchanb,wavebound,detcount,errcount)
                call change_ext(fname,'selfcorbinned')
                call w_diag_file(fname,nchan,wavebin,corrbin,errcount)
!                  call w_diag_file(fname,nxdcs,wavetemp
!     1,corrtemp,errcount)
            endif
        end do
        return
        
    END subroutine do_self_corr         
        
!***********************************************************************************
!*
!*    get_spec.for
!*
!*    A K Soper, December 1999
!*
!*    puts a temporary array back into the spectrum array
!*
!***********************************************************************************
    subroutine get_spec(ntype,is,tarray,earray)

        use run_par
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
!c
!c internal variables
!c
    integer ic,is,index        !internal indices
    integer ntype            !vanadium, sample, can, furnace, etc.
    real, dimension(:)              :: tarray
        real, dimension(:), optional    :: earray!temporary array to store data
!c    write(6,*) nfilebv,runbv(1)
!c
!c get normalised data for this spectrum
!c
    if(ntype.eq.1) then
            do ic=1,nchanv
        tarray(ic)=smovandet(ic,is)
            end do
    else if(ntype.eq.2) then
            do ic=1,nchanv
        tarray(ic)=normbakdet(ic,is)
            end do
    else if(ntype.eq.3) then
            do ic=1,nchans
        tarray(ic)=normbakdet(ic,is)
            end do
        else if(ntype.gt.3) then
            index=ntype-3
            do ic=1,nchans
                tarray(ic)=normsamdet(ic,is,index)
            end do
    end if
        if(present(earray)) then
!c
!c put back corrected errors for this spectrum
!c
            if(ntype.eq.1) then
                do ic=1,nchanv
                    earray(ic)=errvandet(ic,is)
                end do
            else if(ntype.eq.2) then
                do ic=1,nchanv
                    earray(ic)=errbakdet(ic,is)
                end do
            else if(ntype.eq.3) then
                do ic=1,nchanv
                    earray(ic)=errbakdet(ic,is)
                end do
            else if(ntype.gt.3) then
                index=ntype-3
                do ic=1,nchans
                    earray(ic)=errsamdet(ic,is,index)
                end do
            end if
        end if
    return    
    END subroutine get_spec
    
!***********************************************************************************
!*
!*    put_spec.for
!*
!*    A K Soper, December 1999
!*
!*    puts a temporary array back into the spectrum array
!*
!***********************************************************************************
    subroutine put_spec(ntype,is,tarray,earray)

        use run_par
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
!c
!c internal variables
!c
    integer ic,is,index        !internal indices
    integer ntype            !vanadium, sample, can, furnace, etc.
    real, dimension(:)              :: tarray
        real, dimension(:), optional    :: earray!temporary array to store data
!c    write(6,*) nfilebv,runbv(1)
    index=ntype-3
!c
!c get normalised data for this spectrum
!c
    if(ntype.eq.1) then
            do ic=1,nchanv
        smovandet(ic,is)=tarray(ic)
        errvandet(ic,is)=tarray(ic)
            end do
    else if(ntype.eq.2) then
            do ic=1,nchanv
        normbakdet(ic,is)=tarray(ic)
            end do
    else if(ntype.eq.3) then
            do ic=1,nchans
        normbakdet(ic,is)=tarray(ic)
            end do
        else if(ntype.gt.3) then 
            index=ntype-3
            do ic=1,nchans
                normsamdet(ic,is,index)=tarray(ic)
            end do
    end if
        if(present(earray)) then
!c
!c put back corrected errors for this spectrum
!c
            if(ntype.eq.1) then
                do ic=1,nchanv
                    errvandet(ic,is)=earray(ic)
                end do
            else if(ntype.eq.2) then
                do ic=1,nchanv
                    errbakdet(ic,is)=earray(ic)
                end do
            else if(ntype.eq.3) then
                do ic=1,nchanv
                    errbakdet(ic,is)=earray(ic)
                end do
            else if(ntype.gt.3) then
                index=ntype-3
                do ic=1,nchans
                    errsamdet(ic,is,index)=earray(ic)
                end do
            end if
        end if
    return    
    END subroutine put_spec
    
END MODULE process_routines

