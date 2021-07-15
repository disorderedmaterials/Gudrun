!     
! File:   transmission_routines.f90
! Author: aks45
!
! Created on 06 November 2013, 05:53
!

MODULE transmission_routines
    
    implicit none
    
    real, dimension(:), allocatable         :: aln
    integer, parameter                      :: nm=500
    
    CONTAINS
    

!***********************************************************************************
!*
!*    get_trans.FOR
!*
!*    A K Soper, June 2001
!*
!*    calculates total cross section of sample from transmission monitor
!*
!***********************************************************************************
    subroutine get_trans(ntype,wavemin,wavemax,wavestep)
        
        use reallocation_routines
        use inputfilestrings
        use run_par
        use calibration_routines
        use local_data
        use monitor_routines
        use van_par
        use sam_par
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use math_routines
        use interpolation_routines
        use write_routines
        use beam_routines
        use corrections_routines
!c
!c internal variables
!c
        character(len=256)                  :: fname        !dummy file name
        integer                             :: i,is,is1,ic,jc,lenfname        !internal indices
        character(len=256)                  :: run            !filename
        integer                             :: nperrq        !period number of data in raw file
        integer                             :: ntype            !type of file being read
        integer                             :: geom            !geometry of sample
        integer                             :: nan            !no. of containers
        integer                             :: nwav,nwavb,mwav,mchan,mcorrwav,mcorrwavmax      !no. of wavelengths for c/s
        integer                             :: nfirst,nlast,nwavnew
        integer                             :: nsmoo,nord            !no. of smoothings on c/s, order of polynomial to fit transmission ratio
        real, dimension(:), allocatable     :: sammon!(mcorrwav)        !rebinned monitor counts
        real, dimension(:), allocatable     :: samtrans!(mcorrwav)    !rebinned transmission mon counts
        real, dimension(:), allocatable     :: bsammon!(mcorrwav)    !ditto
        real, dimension(:), allocatable     :: bsamtrans!(mcorrwav)    !ditto
        real, dimension(:), allocatable     :: bakmon!(mcorrwav)        !ditto
        real, dimension(:), allocatable     :: baktrans!(mcorrwav)    !ditto
        real, dimension(:), allocatable     :: xfit,yfit!(mcorrwav)  !polynomial fit to transmission
        real                                :: subbi,subbt,subi,subt    !temporary stores
        real, dimension(:), allocatable     :: rad1,rad2!(mcont)    !temporary dimension values
        real, dimension(:), allocatable     :: den!(mcont)        !density of each sample and cont.
        real, dimension(:,:), allocatable   :: sigtl!(mcorrwav,mcont)    !total cross section of sample or cont.
        real                                :: sum,rat,err,x
        real                                :: rot            !rotation angle of sample for flat plate
        real                                :: wavemin,wavemax,wavestep!c
        integer(kind=8) time1,time2,sumtime
!c get the wavelength boundary values which are needed for the transmission.
!c
        nwav=nint((wavemax-wavemin)/wavestep)
        nwavb=nwav+1
        mcorrwavmax=max(nwavb,nchanb)
        write(6,*) 'get_trans> nwav,nwavb,mcorrwavmax ',nwav,nwavb,mcorrwavmax
        call reallocate1d_r(detcount,mcorrwavmax)
        call reallocate1d_r(errcount,mcorrwavmax)
        call reallocate1d_r(wavecorr,mcorrwavmax)
        do ic=1,nwavb
            wavecorr(ic)=wavemin+real(ic-1)*wavestep
        end do

        if(ntype.eq.1) then
            fname=runv(1)
            nperrq=nperv
        else if(ntype.gt.3) then
            is=ntype-3
            fname=runs(1,is)
            nperrq=npers(is)
        endif
        call reallocate1d_r(sammon,nwavb)
        call reallocate1d_r(samtrans,nwavb)
        call reallocate1d_r(bsammon,nwavb)
        call reallocate1d_r(bsamtrans,nwavb)
        call reallocate1d_r(bakmon,nwavb)
        call reallocate1d_r(baktrans,nwavb)
        call reallocate1d_r(xfit,nwavb)
        call reallocate1d_r(yfit,nwavb)
        mchan=nchanb !For compatibility with fixed array versions
        mcorrwav=nwavb !For compatibility with fixed array versions

        if(ntrmon.gt.0) then
            if(ntype.eq.1) then
!c
!c get the vanadium monitor data, and vanadium background data
!c
                call rebinq(nchan,waveimon,smovanmon,nwav,wavecorr,sammon,2,0)
                call rebinq(nchan,wavetmon,smovantrans,nwav,wavecorr,samtrans,2,0)
                call rebinq(nchan,waveimon,smobvanmon,nwav,wavecorr,bsammon,2,0)
                call rebinq(nchan,wavetmon,smobvantrans,nwav,wavecorr,bsamtrans,2,0)
!c
!c set up filename of file to be written
!c
                call change_ext(fname,'rebinq')
                write(6,*) fname
                write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
                do ic=1,nwavb
                    errcount(ic)=0.0
                end do
                call w_diag_file(fname,nwavb,wavecorr,sammon,errcount)
                do ic=1,nwav
                    wavecorr(ic)=0.5*(wavecorr(ic)+wavecorr(ic+1))
                    bakmon(ic)=bsammon(ic)
                    baktrans(ic)=bsamtrans(ic)
                end do

            else if(ntype.gt.3) then

                is=ntype-3
!c get the sample monitor data and get the sample background monitor data
                do ic=1,nchanb
                    detcount(ic)=smosammon(ic,is)
                    errcount(ic)=smosamtrans(ic,is)
                end do
                call rebinq(nchan,waveimon,detcount,nwav,wavecorr,sammon,2,0)
                call rebinq(nchan,wavetmon,errcount,nwav,wavecorr,samtrans,2,0)
                call rebinq(nchan,waveimon,smobsammon,nwav,wavecorr,bsammon,2,0)
                call rebinq(nchan,wavetmon,smobsamtrans,nwav,wavecorr,bsamtrans,2,0)
                do ic=1,nwav
                    wavecorr(ic)=0.5*(wavecorr(ic)+wavecorr(ic+1))
                    errcount(ic)=0.0
                end do
                do ic=1,nwav
                    bakmon(ic)=bsammon(ic)
                    baktrans(ic)=bsamtrans(ic)
                end do
            endif
!c
!c integrate each empty instrument monitor, and calculate time independent 
!c fraction for each 
!c
            subbi=integ(nwav,wavecorr,bakmon,mcorrwav)*fraci
            subbt=integ(nwav,wavecorr,baktrans,mcorrwav)*fract
!c
!c now integrate sample and background monitors, and subtract the corresponding
!c time independent backgrounds
!c
            subi=integ(nwav,wavecorr,sammon,mcorrwav)*fraci
            if(subbi.gt.0.0) then
                subt=(subi*subbt)/subbi
            else
                subt=0.0
            endif
            do ic=1,nwav
                sammon(ic)=sammon(ic)-subi
                samtrans(ic)=samtrans(ic)-subt
            end do
!c
!c repeat this for the sample background
!c
            subi=integ(nwav,wavecorr,bsammon,mcorrwav)*fraci
            if(subbi.gt.0.0) then
                subt=(subi*subbt)/subbi
            else
                subt=0.0
            endif
            do ic=1,nwav
                bsammon(ic)=bsammon(ic)-subi
                bsamtrans(ic)=bsamtrans(ic)-subt
            end do
!c
!c form ratios of transmission monitor to incident monitor, then ratio sample
!c to sample background
!c
            do ic=1,nwav
!c            write(6,*) samtrans(ic),sammon(ic),bsamtrans(ic),bsammon(ic)
                if(sammon(ic).gt.0.0.and.bsammon(ic).gt.0.0) then
                    samtrans(ic)=samtrans(ic)/sammon(ic)
                    bsamtrans(ic)=bsamtrans(ic)/bsammon(ic)
                else
                    samtrans(ic)=0.0

                    bsamtrans(ic)=0.0
                endif
                if(bsamtrans(ic).gt.0.0) then
                    samtrans(ic)=samtrans(ic)/bsamtrans(ic)
                endif
                if(samtrans(ic).le.0.0) samtrans(ic)=0.0
            end do
            !Just in case of binning errors, we need to be sure there are no zeros in these data. If there are then only use the first
            !non-zero region
            call trim_zeros(1,nwav,nfirst,nlast,samtrans)
            ic=0
            do jc=nfirst,nlast
                ic=ic+1
                wavecorr(ic)=wavecorr(jc)
                samtrans(ic)=samtrans(jc)
            end do
            nwav=nlast-nfirst+1
            nwavb=nwav+1
!c
!c set up filename of file to be written
!c
            call change_ext(fname,'transnofit')
            write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
            lenfname=len_trim(fname)
            write(6,*) 'get_trans> ',fname(1:lenfname),ntype
            do ic=1,nwav
                errcount(ic)=0.0
            end do
            call w_diag_file(fname,nwav,wavecorr,samtrans,errcount)
!c set up xfit scale in the range 0 - 1.0
            do ic=1,nwav
!            write(6,*) ic,nwav,size(xfit),size(yfit),size(wavecorr),size(errcount)
                xfit(ic)=(wavecorr(ic)-wavecorr(1))/(wavecorr(nwav)-wavecorr(1))
!            xfit(ic)=0.0
!            write(6,*) ic,xfit(ic)
                yfit(ic)=samtrans(ic)
!            write(6,*) ic,xfit(ic),yfit(ic)
            end do
!c Fit polynomial to data
            nord=int(0.5*wavecorr(nwav))+1
            write(6,*) nwav,nord,ncont
            call system_clock(time1)
            call polyfit_nofix(nwav,xfit,yfit,nord,samtrans,1000000,0.001,0.0001,mcorrwav)
            call system_clock(time2)
            sumtime=(time2-time1)/1000
            write(6,*) 'get_trans> Time for polyfit> ',sumtime
        else
            do ic=1,nwav
                wavecorr(ic)=0.5*(wavecorr(ic)+wavecorr(ic+1))
                samtrans(ic)=1.0
            end do
        endif
!c
!c Save it in the appropriate array
!c
        if(ntype.eq.1) then
            nlambtransv=nwav
            call reallocate1d_r(vlambtrans,nlambtransv)
            call reallocate1d_r(vantransmission,nlambtransv)
            do ic=1,nwav
                vlambtrans(ic)=wavecorr(ic)
                vantransmission(ic)=samtrans(ic)
            end do
            nlambv=nwav
        else
            is=ntype-3
            nlambtranss=nwav
            call reallocate1d_r(slambtrans,nlambtranss)
            call reallocate2d_r(samtransmission,nlambtranss,ncont)
            do ic=1,nwav
                slambtrans(ic)=wavecorr(ic)
                samtransmission(ic,is)=samtrans(ic)
            end do
        endif
!c
!c set up filename of file to be written
!c
        call change_ext(fname,newext='trans')
        write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
        lenfname=len_trim(fname)
        write(6,*) 'get_trans> ',fname(1:lenfname),ntype
        do ic=1,nwav
            errcount(ic)=0.0
        end do
        call w_diag_file(fname,nwav,wavecorr,samtrans,errcount)
        if(allocated(sammon)) deallocate(sammon,samtrans,bsammon,bsamtrans,bakmon,baktrans,xfit,yfit)
    return
        
    end subroutine get_trans
        
!***********************************************************************************
!*
!*    get_tran_cs.FOR
!*
!*    A K Soper, June 2001
!
!*     11 November 2009. Have separated the calculation of the transmission for each
!*     sample from the calculation of the transmission cross section
!*
!*    calculates total cross section of sample from transmission monitor
!*
!***********************************************************************************
    subroutine get_trans_cs(ntype)

        use reallocation_routines
        use inputfilestrings
        use run_par
        use calibration_routines
        use local_data
        use monitor_routines
        use van_par
        use sam_par
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use math_routines
        use interpolation_routines
        use write_routines
        use beam_routines
!c
!c internal variables
!c
    character(len=256)                      :: fname        !dummy file name
    integer                                 :: i,is,is1,ic        !internal indices
    character(len=256)                      :: run            !filename
    integer                                 :: nperrq        !period number of data in raw file
    integer                                 :: ntype            !type of file being read
    integer                                 :: geom            !geometry of sample
    integer                                 :: nan            !no. of containers
    integer                                 :: nwav            !no. of wavelengths for c/s
    integer                                 :: nfirst,nlast,nwavnew
    integer                                 :: nsmoo            !no. of smoothings on c/s
    real, dimension (:), allocatable        :: samtrans !(mcorrwav)    !rebinned transmission mon counts
    real, dimension(:), allocatable         :: rad1,rad2 !(mcont)    !temporary dimension values
    real, dimension(:), allocatable         :: den!(mcont)        !density of each sample and cont.
    real, dimension(:,:), allocatable       :: sigtl!(mcorrwav,mcont)    !total cross section of sample or cont.
    real, dimension(:), allocatable         :: temptrans!(mcorrwav)   !temporary array of sample transmissions
    real                                    :: sum,rat,err,x,xx,rats,captav
    real                                    :: rot            !rotation angle of sample for flat plate
    
        if(ntype.eq.1) then
            run=runv(1)
            nperrq=nperv
            geom=vgeom
            if(geom.eq.2) rot=vrot
!c
!c get the wavelength boundary values which are needed for the total 
!c cross section, and the transmission data
!c
            nwav=nlambtransv
            call reallocate1d_r(samtrans,nwav)
            call reallocate1d_r(wavecorr,nwav)
            do ic=1,nwav
                wavecorr(ic)=vlambtrans(ic)
                samtrans(ic)=vantransmission(ic)
            end do
!c
!c now get the relevant sample data
!c
            nan=1
            call reallocate1d_r(rad1,nan)
            call reallocate1d_r(rad2,nan)
            call reallocate1d_r(den,nan)
            call reallocate2d_r(sigtl,nwav,nan)
            rad1(nan)=vdimen1
            rad2(nan)=vdimen2
            den(nan)=vrho
            do ic=1,nwav
                sigtl(ic,nan)=vtscat(ic)
            end do
        else if(ntype.gt.3) then
            is=ntype-3
            run=runs(1,is)
            nperrq=npers(is)
            geom=sgeom
            if(geom.eq.2) rot=srot(1)
!c
!c get the wavelength boundary values which are needed for the total 
!c cross section, and the sample transmission data
!c
            nwav=nlambtranss
            call reallocate1d_r(samtrans,nwav)
            call reallocate1d_r(wavecorr,nwav)
            do ic=1,nwav
                wavecorr(ic)=slambtrans(ic)
                if(is.eq.ncont) then
                    samtrans(ic)=samtransmission(ic,is)
                else

!c Get ratio of transmissions for all but the outermost containers

                    samtrans(ic)=samtransmission(ic,is)/samtransmission(ic,is+1)
                endif
            end do
!c
!c now get the relevant sample data
!c
            nan=ncont-is+1
            call reallocate1d_r(rad1,nan)
            call reallocate1d_r(rad2,nan)
            call reallocate1d_r(den,nan)
! It is conceivable that either the sample or containers will have cross sections read in from a file, in which case we should use the wavelength
! scale of the lowest order container
            nwav=nlambs
            call reallocate2d_r(sigtl,nwav,nan)
            call reallocate1d_r(temptrans,nwav)
            nan=0
            do i=is,ncont
                nan=nan+1
                rad1(nan)=sdimen1(i)
                rad2(nan)=sdimen2(i)
                den(nan)=srho(i)
                if(i.eq.1) den(nan)=den(nan)/abs(tweak(i))
                do ic=1,nlambs
                    xx=slamb(ic)
                    temptrans(ic)=ainter(wavecorr,samtrans,xx,nlambtranss)
                    sigtl(ic,nan)=stscat(ic,i)
                end do
            end do
        endif
!c
!c o.k. now calculate the cross section
!c
        if(geom.eq.1) then
            call transcyl(nan,nwav,temptrans,rad1,rad2,den,sigtl,nwav,nan)
        else if(geom.eq.2) then
            call transflat(nan,nwav,temptrans,rad1,rad2,rot,den,sigtl,nwav,nan)
        endif
!c
!c step through the values, and use only range for which the cross sections
!c are .gt. 0
!c
        nfirst=1
        nlast=1
        do ic=1,nwav
            if(sigtl(ic,1).gt.0.0) then
                if(nfirst.eq.0) nfirst=ic
                nlast=ic
            endif
        end do
        nwavnew=nlast-nfirst+1
!c
!c apply a little smoothing to the data to correct for fluctuations
!c
!    do ic=1,nwav
!            detcount(ic)=sigtl(ic,1)
!    end do
!    nsmoo=nwavnew/40+1
!    call smood(nfirst,nlast,nsmoo,detcount,mchan)
!    do ic=1,nwav
!        sigtl(ic,1)=detcount(ic)
!    end do
!c
!c write the data
!c
        write(fname,105) run(1:index(run,'.')-1),nperrq
105     format(a,'.mut',i2.2)
        open(10,file=fname,status='unknown')
        write(10,106) '#',nwavnew
106     format(a1,i5)
        if(ntype.eq.1) then
            captav=vcaptav
        else
            captav=scaptav(ntype-3)
        end if
        do ic=nfirst,nlast
            x=wavecorr(ic)
            rat=sigtl(ic,1)
            rats=max(0.0,(rat-x*captav/1.798))
            write(10,101) x,rat,rats
101     format(3(1x,e13.6))
        end do
        close(10)
!c
!c and save it in the appropriate array
!c
        if(ntype.eq.1) then
            nlambv=nwavnew
            do ic=nfirst,nlast
                i=ic-nfirst+1
                vlamb(i)=vlamb(ic)
                vtscat(i)=sigtl(ic,1)
                vsscat(i)=max(0.0,(vtscat(i)-vlamb(i)*vcaptav/1.798)) !Prevent scattering cross section from going below zero!
            end do
        else
            is=ntype-3
            nlambs=nwavnew
            do ic=nfirst,nlast
                i=ic-nfirst+1
                slamb(i)=slamb(ic)
                stscat(i,is)=sigtl(ic,1)
                ssscat(i,is)=max(0.0,(stscat(i,is)-slamb(i)*scaptav(is)/1.798))
            end do
        endif
        deallocate(rad1,rad2,den,sigtl,temptrans)
        return
        end subroutine get_trans_cs
    
        subroutine transcyl(nan,nwav,samtrans,rad1,rad2,den,sigtl,mcorrwav,mcont)

        use reallocation_routines
        use inputfilestrings
        use run_par
        use calibration_routines
        use local_data
        use monitor_routines
        use van_par
        use sam_par
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use math_routines
        use interpolation_routines
        use write_routines
        use beam_routines
        integer                             :: nan,mcont            !no. of containers
        integer                             :: nwav,mcorrwav            !no. of wavelengths for c/s
        real                                :: samtrans(mcorrwav)    !transmission monitor ratio
        real                                :: rad1(mcont),rad2(mcont)    !sample dimension values
        real                                :: den(mcont)        !density of each sample and cont.
        real                                :: sigtl(mcorrwav,mcont)    !total cross section of sample or cont.
        real, dimension(:), allocatable     :: amuc!(mcont)        !temporary attenuation coefficients
        integer                             :: ic,in
        real                                :: tra,amu2,x

        call reallocate1d_r(amuc,nan)
        do ic=1,nwav
            x=wavecorr(ic)
            write(6,*) nan,ic,x
            if(nan.gt.1) then
                do in=2,nan
                    amuc(in)=den(in)*sigtl(ic,in)
                end do
            endif
            call lmomnt(nan,amuc,rad1,rad2)
            tra=samtrans(ic)
            call trans2(tra,amu2)
            write(6,*) ic,x,tra,amu2
            sigtl(ic,1)=amu2/den(1)
        end do
        return
    end subroutine transcyl
        
    subroutine lmomnt(nan,amuc,rad1,rad2)

        use reallocation_routines
        use beam_routines

        integer nan
        real amuc(nan)
        real rad1(nan)
        real rad2(nan)
        integer nx,i,j
        real al0,almu0,xst,xhalf,x,al1,al2,al,sum,cont,cont1,pr
        
        call reallocate1d_r(aln,nm)
        nx=100
!c
!c zero all moments
!c
        do i=1,nm
            aln(i)=0.
        end do
        al0=0.
        almu0=0.
!c
!c step across beam and calculate distance through sample
!c and profile value
!c
        xst=(a-b)/nx
        xhalf=0.5*xst
        prstep=(a-b)/(nprof-1)
        do i=1,nx
            x=a-i*xst+xhalf
            al1=distt(rad1(1),x)
            al2=distt(rad2(1),x)
            al=al2-al1
!c
!c if containers are present, then calculate the attenuation through these
!c 
            sum=0
            if(nan.gt.1) then
                do j=2,nan
                    al1=distt(rad1(j),x)
                    al2=distt(rad2(j),x)
                    sum=sum+amuc(j)*(al2-al1)
                end do
            endif
            pr=probe(x,profil,nprof,prstep,a,b)
!c
!c generate contributions to moments
!c
            cont=pr*xst
            al0=al0+cont
            cont1=cont*exp(-sum)
            almu0=almu0+cont1
            do j=1,nm
                cont1=al*cont1/j
                aln(j)=aln(j)+cont1
            end do
        end do
        do i=1,nm
            aln(i)=aln(i)/almu0
        end do
!c    write(6,101) almu0,(aln(i),i=1,50)
101     format((5(1x,e13.6)))
        return
    end subroutine lmomnt

    function distt(rad,x)
        
        real distt,rad,x,y
        
    if(abs(x).gt.rad) then
            distt=0.0
        else
            y=rad*rad-x*x
            distt=2.*sqrt(y)
        end if
    return
    end function distt
    
    subroutine trans2(tra,amu2)
        
        real tra,amu2
        real acur,amu1,amun,tr1,tr2,diff1,diff2,grad
        integer ntry,nc1,nc2
!c
!c Find the value of amu2 consistent with transmission value tra by successive
!c approximation
!c

        acur=0.001
!c
!c generate trial values
!c
        ntry=0
        amu1=1.0
        tr1=sexp(nc1,amu1)
        amu2=amu1*0.5
        if(tra.lt.tr1) amu2=2.0*amu1
        tr2=sexp(nc2,amu2)
        do while (abs(tr2-tr1).ge.acur.and.ntry.lt.200)
            if(abs(tra-tr2).lt.acur) then
                return
            else
                ntry=ntry+1
!        write(6,*) ntry,nc2,tra,amu1,amu2,tr1,tr2
                !Predict the best likely new value
                grad=(amu2-amu1)/(tr2-tr1)
                amun=grad*(tra-tr1)+amu1
                !Test to see which is better - save the best
                diff1=abs(amun-amu1)
                diff2=abs(amun-amu2)
                if(diff2.lt.diff1) then
                    amu1=amu2
                    tr1=tr2
                    nc1=nc2
                end if
            end if
            !Try the new value
            amu2=amun
            tr2=sexp(nc2,amu2)
        end do
        return
    end subroutine trans2

    function sexp(nmuse,amu)

        real sexp,amu
        integer nmuse,nref,i
        real sum
        
    nmuse=int(amu*aln(1)/0.05)+1
    if(nmuse.lt.10) nmuse=10
    if(nmuse.gt.nm) nmuse=nm
    nref=nm
    sum=0.
    do i=1,nm
            sum=aln(nref)-amu*sum
!c    write(6,999) amu,aln(nref),sum
999    format(5(1x,e13.6))
            nref=nref-1
        end do
    sexp=1.-amu*sum
!c    write(6,999) amu,almu0,sexp
    return
    end function sexp

    subroutine transflat(nan,nwav,samtrans,rad1,rad2,srot,den,sigtl,mcorrwav,mcont)
    
        use reallocation_routines
        use math_routines

        integer                             :: nan,mcont            !no. of containers
        integer                             :: nwav,mcorrwav            !no. of wavelengths for c/s
        real                                :: samtrans(mcorrwav)    !transmission monitor ratio
        real                                :: rad1(mcont),rad2(mcont)    !sample dimension values
        real                                :: den(mcont)        !density of each sample and cont.
        real                                :: sigtl(mcorrwav,mcont)    !total cross section of sample or cont.
        real, dimension(:), allocatable     :: amuc!(mcont)        !temporary attenuation coefficients
        real                                :: srot            !rotation angle of sample (degrees)
        integer                             :: ic
        real                                :: fact,tra,thick,amu2
        
        call reallocate1d_r(amuc,nwav)
        fact=den(1)*cosd(srot)
        do ic=1,nwav
            tra=samtrans(ic)
            thick=rad2(1)+rad1(1)
            call transf(tra,thick,amu2)
!c    write(6,999) x,tra,amu2
!999    format(5(1x,e13.6))
            sigtl(ic,1)=amu2/fact
        end do
        return
        
    end subroutine transflat
    
    subroutine transf(tra,rad,amu2)
        
        real tra,rad,amu2
        
        if(tra.gt.0.0) then
            amu2=-alog(tra)/rad
        else
            amu2=0.0
        endif
        return

    end subroutine transf
    
END MODULE transmission_routines


