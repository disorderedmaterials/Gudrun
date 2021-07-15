!     
! File:   attenuation_correction_routines.f90
! Author: aks45
!
! Created on 19 November 2013, 09:07
!

MODULE abs_corr_azi
    
    implicit none
    
    integer                                     :: nwavabs        !number of wavelengths for a.c.
    integer                                     :: nangabs        !number of angles for a.c.
    integer                                     :: naziabs        !number of azimuthal angles for a.c.
    real, dimension(:), allocatable             :: wavabs!(mcorrwav)        !wavelengths for a.c.
    real, dimension(:), allocatable             :: angabs!(mcorrang)            !angles for a.c.
    real, dimension(:), allocatable             :: aziabs!(mcorrang)
    real, dimension(:,:,:,:), allocatable       :: abscor!(mcorrwav,mcorrang,mcorrang,mcont1) !full corrections array
    real, dimension(:), allocatable             :: rad1,rad2!(mcont)    !temporary dimension values
    real, dimension(:), allocatable             :: den!(mcont)        !density of each sample and cont.
    real, dimension(:), allocatable             :: captcs!(mcont)            !capture cross section of sample or cont
    real, dimension(:,:), allocatable           :: sigsl!(mcorrwav,mcont)    !total cross section of sample or cont.
    real, dimension(:,:), allocatable           :: sigtl,sigtl_out!(mcorrwav,mcont)    !total cross section of sample or cont.
    real, dimension(:,:,:,:), allocatable       :: abscort!(mcorrwav,mcorrang,mcorrang,mcont) !temporary attenuation correction
    real, dimension(:,:), allocatable           :: MUS,MUT,mut_out!(mcorrwav,mcont)
    real, dimension(:), allocatable             :: amus,amut,amut_out
    real, dimension(:), allocatable             :: abstemp
    real                                        :: theta,theta1,phid,azi,pi,piconv,height,astep,areas,areac
    
END MODULE abs_corr_azi
    
MODULE attenuation_correction_routines
    
!    use reallocation_routines
!    use inputfilestrings
!    use run_par
!    use calibration_routines
!    use abs_corr_azi
!    use beam_routines
!    use sam_par
!    use corrections_routines

    implicit none
    
    CONTAINS
    
!***********************************************************************************
!*
!*    get_abs_corr.FOR
!*
!*    A K Soper, December 1999
!*
!*    gets attenuation correction for any run
!*
!***********************************************************************************
    subroutine get_abs_corr(run,nperrq,forcecalc)

    use reallocation_routines
    use inputfilestrings
    use corrections_routines
    use abs_corr_azi
    use sam_par

!c internal variables
!c
    character(len=256) fname        !dummy file name
    integer i,ib,j,k,ia        !internal indices
    character(len=256) run            !raw filename
    integer nperrq        !period number of data in raw file
    integer ncorr            !no. of corrections
    integer ierr,forcecalc        != 1 forces calculation of corrs.
    real aa
        real, dimension(:), allocatable     ::bb!(mcorrang,mcont1)        !temporary real values

        write(6,*) 'Got to Get_ABS_CORR'
!c
!c force the calculation of the corrections if desired
!c
        ierr=1 !Signal that corrections need to be calculated in case corrections file has problems
        if(forcecalc.eq.0) then
!c
!c set up total number of corrections expected to be read
!c
            ncorr=(ncont*(ncont+1))/2
!c
!c set up filename of file to be read
!c
            fname=run
            call change_ext(fname,'abs')
            write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
!c
!c try to open attenuation correction file. If an error occurs programme
!c redirects to calculate the attenuation  correction 
!c
            open(10,file=fname,status='old',iostat=ierr)
            if(ierr.eq.0) then
                read(10,*,iostat=ierr) nwavabs,nangabs,naziabs
                write(6,*) nwavabs,nangabs,naziabs
            endif
            if(ierr.eq.0) then
                !Allocate the space required
                call reallocate1d_r(wavabs,nwavabs)
                call reallocate1d_r(angabs,nangabs)
                call reallocate1d_r(aziabs,naziabs)
                call reallocate4d_r(abscor,nwavabs,nangabs,naziabs,ncorr)
                read(10,*,iostat=ierr) (aziabs(k),k=1,naziabs)
            endif
            ia=0
            do while (ia.lt.nangabs.and.ierr.eq.0)
                ia=ia+1
        read(10,*,iostat=ierr) ib,angabs(ia)
        j=0
                do while (j.lt.nwavabs.and.ierr.eq.0)
                    j=j+1
                    read(10,*,iostat=ierr) wavabs(j),((abscor(j,ib,k,i),i=1,ncorr),k=1,naziabs)
        end do
            end do
            close(10)
            if (ia.lt.nangabs.or.j.lt.nwavabs) then
                write(6,1002) fname(1:len_trim(fname))
1002    format(' get_abs> ERROR ** unexpected end of file in ',a)
                stop
                return
            end if
        endif
        if(ierr.ne.0) then
!
!Calculate attenuation 
!c correction from input sample parameters
!c
            call calc_abs(run,nperrq)
            return
        end if
    end subroutine get_abs_corr
        
!***********************************************************************************
!*
!*    calc_abs.FOR
!*
!*    A K Soper, December 1999
!*
!*    calculates attenuation correction for any run
!*
!***********************************************************************************
    subroutine calc_abs(run,nperrq)
    use reallocation_routines
    use inputfilestrings
    use corrections_routines
    use abs_corr_azi
    use sam_par
!c
!c internal variables
!c
    character(len=256) run        !file name for m.s.corrections
    character(len=256) fname        !file name for m.s.corrections
    character(len=32) formatout
    integer i,ib,j,ind,il,iang,k    !internal indices
    integer nperrq        !period number of data in raw file
    integer nan            !temporary array of number of annuli
    integer ncorr            !no. of absorption corrections per lamb.
    integer nblock        !number of blocks in output file
    integer lrec        !record length 
!c
!c set up total number of corrections expected to be read
!c
        ncorr=(ncont*(ncont+1))/2
!c
!c set up some initial arrays
!c
        nwavabs=nlambs
        call reallocate1d_r(wavabs,nwavabs)
        do i=1,nwavabs
            wavabs(i)=slamb(i)
        end do
        call reallocate1d_r(rad1,ncont)
        call reallocate1d_r(rad2,ncont)
        call reallocate1d_r(den,ncont)
        call reallocate1d_r(captcs,ncont)
        call reallocate2d_r(sigtl,nwavabs,ncont)
!c
!c set up angles for corrections
!c
    call set_corr_ang(sgeom,srot(1),angabs,aziabs,nangabs,naziabs)
        call reallocate4d_r(abscor,nwavabs,nangabs,naziabs,ncorr)
!c
!c calculate the corrections
!c
    ind=0
    do j=1,ncont
            nan=0
            do i=j,ncont
        nan=nan+1
        rad1(nan)=sdimen1(i)
        rad2(nan)=sdimen2(i)
        den(nan)=srho(i)
!c
!c for sample include tweak factor in sample density
!c
        if(i.eq.1) den(nan)=den(nan)/abs(tweak(i))
        captcs(nan)=scaptav(i)
        do il=1,nwavabs
                    sigtl(il,nan)=stscat(il,i)
        end do
            end do
            if(sgeom.eq.1) then
        call cylabstof(nan,sheight)
            else if(sgeom.eq.2) then
        call fltabstof(nan,srot(1))
            endif
!c
!c save the attenuation correction for this set of sample and containers
!c
!c the corrections are always stored in order of increasing radius.
!c
            do i=1,nan
        ind=ind+1
        do iang=1,nangabs
                    do il=1,nwavabs
            do k=1,naziabs
                            abscor(il,iang,k,ind)=abscort(il,iang,k,i)
            end do
                    end do
        end do
            end do
        end do
    ncorr=ind
    write(6,*) ncorr,' attenuation corrections created in calc_abs'
!c
!c set the output format
!c
    nblock=naziabs*ncorr+1
    write(formatout,332) '(',nblock,'(1x,e12.5))'
332    format(a1,i3.3,a11)
    write(6,*) formatout
    lrec=nblock*13
    if(lrec.lt.80) lrec=80
!c
!c set up filename of file to write
!c
        fname=run
        call change_ext(fname,'abs')
        write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
!
!c Open the attenuation corrections file and write it. 
!c
    open(10,file=fname,status='unknown',form='formatted',recl=lrec)
    write(10,*) nwavabs,nangabs,naziabs
    write(10,formatout) (aziabs(k),k=1,naziabs)
    do ib=1,nangabs
            write(10,*) ib,angabs(ib)
            do j=1,nwavabs
                write(10,formatout) wavabs(j),((abscor(j,ib,k,i),i=1,ncorr),k=1,naziabs)
            end do
    end do
    close(10)
    return

    end subroutine calc_abs
    
    subroutine CYLABSTOF(nan,theight)
        
        use reallocation_routines
        use corrections_routines
        use abs_corr_azi
        use sam_par
        use beam_routines

        real                        :: theight,fac
        integer                     :: nan,i,il,iazi,j,ir,k,ms
!c
!c    ORIGINAL from AKS
!c    modified 27/03/2001 to become a subroutine within GUDRUN
!c
!c    modified 5/02/2003 to include corrections for azimuthal detectors
!c
!c    Note that the arrays angabs and aziabs are in sample cylindrical coordinates
!c     (as opposed to detector coordinates) so that angabs is relative a z-axis in the 
!c    vertically upwards direction
!c

        PI=4.0*atan(1.0)
    astep=stepa
    height=theight
        call reallocate2d_r(mus,nwavabs,nan)
        call reallocate2d_r(mut,nwavabs,nan)
        call reallocate2d_r(sigsl,nwavabs,nan)
        call reallocate1d_r(amus,nan)
        call reallocate1d_r(amut,nan)
        call reallocate4d_r(abscort,nwavabs,nangabs,naziabs,nan)
        call reallocate1d_r(abstemp,nan)
        do i=1,nan
!C
!C CALCULATE  SCATTERING C/S AT EACH WAVELENGTH ASSUMING 1/LAMBDA ABSORPTION
!C
            DO IL=1,nwavabs
                FAC=wavabs(IL)/1.7979
                SIGSL(IL,I)=SIGTL(IL,I)-FAC*captcs(I)
                MUS(IL,I)=DEN(I)*SIGSL(IL,I)
                MUT(IL,I)=DEN(I)*SIGTL(IL,I)
            end do
        end do
        MS=(RAD2(1)-RAD1(1)+0.0001)/ASTEP
        IF (MS.LT.1) MS=1
        do iazi=1,naziabs
            azi=aziabs(iazi)*pi/180.0
            DO I=1,nangabs
                THETA=angabs(I)*pi/180.0
                DO J=1,nwavabs
                    DO IR=1,NAN
                        AMUS(IR)=MUS(J,IR)
                        AMUT(IR)=MUT(J,IR)
                    end do
                    CALL ACYL(MS,nan)
                    do k=1,nan
                        abscort(j,i,iazi,k)=abstemp(k)
                    end do
                end do
            end do
        end do
    return
        
    END subroutine cylabstof
      
    SUBROUTINE ACYL(MS1,nan)
        
        use abs_corr_azi
        use beam_routines

        integer                 :: ms,ms1,nan,i
        real                    :: aaaa,aaab,areaa,areab,areas1
!C
!C     DETERMINE MAXIMUM OF THE PROFILE FUNCTION
!C
        ms=ms1
        DO I=1,NAN
            abstemp(i)=1.0
            IF(amus(I).gt.0.0) then
!C
!C  NO STEPS ARE CHOSEN SO THAT STEP WIDTH IS THE SAME FOR ALL ANNULI
!C
                MS=MS1*(RAD2(i)-RAD1(i))/(RAD2(1)-RAD1(1))
                if(ms.lt.1) ms=1 
                CALL SUMROM(nan,AAAA,AREAA,I,A,RAD1(I),RAD2(I),MS)
                CALL SUMROM(nan,AAAB,AREAB,I,B,RAD1(I),RAD2(I),MS)
                AREAS1=SIGN(AREAA,A)-SIGN(AREAB,B)
                if(areas1.gt.0.) abstemp(i)=(SIGN(AAAA,A)-SIGN(AAAB,B))/areas1
            end if
        end do
        RETURN
    END subroutine acyl

    SUBROUTINE SUMROM(nan,AAA,AREA,NRAD,A2,R1,R2,MS)

        use reallocation_routines
        use corrections_routines
        use beam_routines
        use abs_corr_azi

        real aaa,area,a2,r1,r2
        integer nan,nrad,ms
        
        REAL, dimension(:), allocatable       :: LIS,LSS!(5),
        real :: LIST,LISN,LSST,LSSN,a3,detfac,omgadd,width,p1,p2,x,oad,rstep,radd,r,omegst,omegad
        real :: areay,sum1,sum2,arsum,omega,d,prob,o,angle,path,alim
        integer :: m,m1,m2,n,nomeg,i,ii,j
        logical :: secondloop

        call reallocate1d_r(lis,nan)
        call reallocate1d_r(lss,nan)
!c
!c sine of detector azimuthal angle
!c
        detfac=abs(sin(theta))
        A3=A2
        OMGADD=0.
        IF(A3.LT.0.) OMGADD=PI
        A3=ABS(A3)
        IF(A3.GE.R2) ALIM=R2
        IF(A3.LT.R2) ALIM=A
        IF(A3.LE.R1) ALIM=R1
        M1=int(real(MS)*((ALIM-R1)/(R2-R1)+0.0001))
        M2=MS-M1
!        PSUM=0.
!C
!C INTEGRATE THE BEAM PROFILE FROM 0 TO ALIMIT
!C
!        WIDTH=A3
!        IF(A3.GE.R2) WIDTH=R2
!        ASTEP=WIDTH/40
!        P1=0.5*PROBE(0.,PROFIL,NPROF,PRSTEP,A,B)
!        DO I=1,40
!            X=I*ASTEP
!            X=SIGN(X,A2)
!            P2=0.5*PROBE(X,PROFIL,NPROF,PRSTEP,A,B)
!            PSUM=PSUM+P1+P2
!            P1=P2
!        end do
!   10 CONTINUE
!        PSUM=PSUM*ASTEP
        OAD=PI-azi
!        IF(M1.EQ.0) GO TO 41
        AAA=0.
        AREA=0.
        secondloop=.true.!In general will perform two loops, one for M1 and one for M2
        if(m1.gt.0) then
            N=M1
            RSTEP=(ALIM-R1)/M1
            RADD=-0.5*RSTEP+R1
        else if(M2.gt.0) then
            N=M2
            RSTEP=(R2-ALIM)/M2
            RADD=-0.5*RSTEP+ALIM
            secondloop=.false.!Only loop once in this case
        else
            return
        end if
!   31 DO 30 M=1,N
        m=0
        do while (m.lt.n)
            m=m+1
            R=M*RSTEP+RADD
            NOMEG=PI*R/RSTEP
            OMEGST=PI/NOMEG
            OMEGAD=-0.5*OMEGST+OMGADD
            AREAY=R*RSTEP*OMEGST*AMUS(NRAD)
            SUM1=0.
            SUM2=0.
            ARSUM=0.
            I=0
            do while (i.lt.nomeg)
                i=i+1
!   35 IF(I.GT.NOMEG) GO TO 36
                OMEGA=I*OMEGST+OMEGAD
                D=R*SIN(OMEGA)
!                IF(ABS(D).GT.A3) GO TO 101
                IF(ABS(D).le.A3) then
!C
!C DETERMINE A PROFILE VALUE FOR THIS 'D' VALUE
!C
                    PROB=PROBE(D,PROFIL,NPROF,PRSTEP,A,B)
!C
!C CALCULATE DISTANCE INCIDENT NEUTRON PASSES THROUGH EACH ANNULUS
!C
                    DO J=1,NAN
                        LIST=DIST(R,RAD1(J),OMEGA)
                        LISN=DIST(R,RAD2(J),OMEGA)
                        LIS(J)=LISN-LIST
                    end do
!C
!C CALCULATE DISTANCE SCATTERED NEUTRON PASSES THROUGH EACH ANNULUS
!C
                    O=OMEGA+OAD
!c
!c check this element is visible to the detector
!c (Apparently this was not done in previous versions of the programme!)
!c
                    D=R*sin(O)
                    if(D.lt.A1.and.D.gt.B1) then
                        DO J=1,NAN
                            LSST=DIST(R,RAD1(J),O)
                            LSSN=DIST(R,RAD2(J),O)
                            LSS(J)=LSSN-LSST
                            if(detfac.gt.1.0e-5) then
                                lss(j)=lss(j)/detfac
                            end if
                        end do
!C
!C CALCULATE ABSORPTION FOR PATH THROUGH ALL ANNULI
!C
                        PATH=0
                        DO  II=1,NAN
                            PATH=PATH+AMUT(II)*(LIS(II)+LSS(II))
                        end do
                        if(path.lt.30.0) then
                            SUM1=SUM1+EXP(-PATH/detfac)*PROB
                        endif
                        ARSUM=ARSUM+PROB
                    endif
                else
                    I=NOMEG-I+2
                end if
            end do
            AAA=AAA+SUM1*AREAY
            AREA=AREA+ARSUM*AREAY
            if(m.eq.n) then !Repeat for the remainder of the radius values if relevant.
                secondloop=secondloop.and.m2.ne.0
                if(.not.secondloop) return
                N=M2
                RSTEP=(R2-ALIM)/M2
                RADD=-0.5*RSTEP+ALIM
                M=0
                secondloop=.false. !Make sure we exit for sure after second loop
            end if
        end do
        return
!      GO TO 31
    END subroutine sumrom
    
    subroutine fltabstof(ncan,rot1)

        use reallocation_routines
        use corrections_routines
        use abs_corr_azi
        use sam_par
!c
! to calculate flat plate absorption factors
!c
!c      original   by aks
!c      modified   by wsh    10-3-88   for use with coral
!c      modified by aks 16-7-90 to put out a unity absorption factor
!c if secondary angle is close to 90 degrees.
!c
!c      completely revamped by aks 17-05-01 to become a subroutine of the gudrun
!c suite
!c
!c      modified 5/02/2003 to include situation with azimuthal detectors
!c      modified January 2012 to replace inner and outer dimensions with upstream and downstream slabs
!c      also corrected logic which was flawed in some cases. 
!c

        integer ncan            !no. of containers + sample
        integer ic,ic1,ib,il
        real rot,rot1              !angle of flat plate normal to beam
        real, dimension(:), allocatable     :: ts
        real, dimension(:), allocatable     :: pl_in_dn,pl_out_dn,pl_in_up,pl_out_up!(mcont)
        real fs_up,fs_dn,exp_dn,exp_up      !temporary exponential values
        real tsec,sec1,sec2      !secants of angles
        real ts_in_up,ts_out_up,ts_in_dn,ts_out_dn      !in and out flight path through each container
        real footprintfac,diskfootprint
      
        pi=4.0*atan(1.0)
        piconv=pi/180.
        rot=rot1*piconv

!c Assign the 'out' values assuming elastic values

        call reallocate2d_r(sigtl_out,nwavabs,ncan)
        call reallocate1d_r(amut,ncan)
        call reallocate1d_r(amut_out,ncan)
        call reallocate1d_r(ts,ncan)
        call reallocate1d_r(pl_in_dn,ncan)
        call reallocate1d_r(pl_out_dn,ncan)
        call reallocate1d_r(pl_in_up,ncan)
        call reallocate1d_r(pl_out_up,ncan)
        call reallocate4d_r(abscort,nwavabs,nangabs,naziabs,ncan)
        do ic=1,ncan
            do il=1,nwavabs
                sigtl_out(il,ic)=sigtl(il,ic)
            end do
        end do
        do ib=1,nangabs
            theta1=angabs(ib)
            sec1=abs(1./cos(rot))
            tsec=theta1
!c
!c tsec is the angle the scattered beam makes with the normal to the sample
!c surface.  if abs(tsec) is close to 90 deg. calculation of absorption
!c coefficients is unreliable
!c
            if(abs(abs(tsec)-90.).gt.0.1) then
                tsec=tsec*piconv
                sec2=1./cos(tsec)
!c
!c first form the sum of path lengths (pl) through the other slabs into and out
!c of each slab. the sample and containers are divided into two slabs, dn and up.
!c up refers to the slab on the upstream side of the centre of the sample (i.e.
!c facing towards the neutron source). dn refers to slabs on the downstream side 
!c of the sample centre. in refers to the flight path into the slab. out refers
!c to the flight out of the slab towards the detector.
!c
                do il=1,nwavabs
                    do ic=1,ncan
                        pl_in_up(ic)=0.0
                        pl_out_up(ic)=0.0
                        pl_in_dn(ic)=0.0
                        pl_out_dn(ic)=0.0
                    end do
                    do ic=1,ncan
                        amut(ic)=den(ic)*sigtl(il,ic)
                        amut_out(ic)=den(ic)*sigtl_out(il,ic)
                        ts(ic)=(rad2(ic)+rad1(ic))
                        ts_in_up=amut(ic)*rad1(ic)*sec1
                        ts_out_up=amut_out(ic)*rad1(ic)*abs(sec2)
                        ts_in_dn=amut(ic)*rad2(ic)*sec1
                        ts_out_dn=amut_out(ic)*rad2(ic)*abs(sec2)

                        do ic1=1,ncan
                            if(sec2.gt.0.0) then
                                if(ic1.lt.ic) then
                                    pl_in_up(ic1)=pl_in_up(ic1)+ts_in_up
                                    pl_out_up(ic1)=pl_out_up(ic1)+ts_out_dn
                                    pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                                    pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_dn
                                else if(ic1.eq.ic) then
                                    pl_out_up(ic1)=pl_out_up(ic1)+ts_out_dn
                                    pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                                else if(ic1.gt.ic) then
                                    pl_out_up(ic1)=pl_out_up(ic1)+ts_out_up+ts_out_dn
                                    pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up+ts_in_dn
                                endif
                            else if(sec2.lt.0.0) then
                                if(ic1.lt.ic) then
                                    pl_in_up(ic1)=pl_in_up(ic1)+ts_in_up
                                    pl_out_up(ic1)=pl_out_up(ic1)+ts_out_up
                                    pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                                    pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_up
                                else if(ic1.eq.ic) then
                                    pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                                    pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_up
                                else if(ic1.gt.ic) then
                                    pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up+ts_in_dn
                                    pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_up+ts_out_dn
                                endif
                            endif
                        end do
                    end do
!      if(ib.eq.1.and.il.eq.1) write(6,*) (pl_in_up(ic),pl_out_up(ic),pl_in_dn(ic),pl_out_dn(ic),ic=1,ncan)
!c
!c calculate absorption integrals for each slab and form attenuation correction.
!c
                    do ic=1,ncan
                        fs_up=f(amut(ic),amut_out(ic),rad1(ic),sec1,sec2)
                        fs_dn=f(amut(ic),amut_out(ic),rad2(ic),sec1,sec2)
                        exp_up=exp(-(pl_in_up(ic)+pl_out_up(ic)))
                        exp_dn=exp(-(pl_in_dn(ic)+pl_out_dn(ic)))
                        abscort(il,ib,1,ic)=(rad1(ic)*exp_up*fs_up+rad2(ic)*exp_dn*fs_dn)/(ts(ic))
                    end do
                end do
            else
!c
!c case where tsec is close to 90 - no attenuation correction can be calculated
!c
                do ic=1,ncan
                    do il=1,nwavabs
                        abscort(il,ib,1,ic)=1.0
                    end do
                end do
            endif
        end do
        return
    end subroutine fltabstof

    function f(amu,amu_out,t,sec1,sec2)
        
        real            :: f,amu,amu_out,t,sec1,sec2
        real*8          :: s,s1,ft
        integer         :: nterms
    
        s=t*(amu*sec1-amu_out*sec2)
        s1=amu_out*t*sec2
        if(abs(s).lt.0.01)then
!c use a simple series expansion
            nterms=5
            ft=1.0-s/real(nterms)
            do while (nterms.gt.2)
                nterms=nterms-1
                ft=1.0-s*ft/real(nterms)
            end do               
        else
            ft=(1.0-exp(-s))/s
        endif
        if(s1.gt.0.0) then
            ft=ft*exp(-s1)
        endif
        f=ft
        return
    end function f

END MODULE attenuation_correction_routines
