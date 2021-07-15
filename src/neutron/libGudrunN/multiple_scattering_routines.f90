!     
! File:   multiple_scattering_routines.f90
! Author: aks45
!
! Created on 10 November 2013, 07:19
!

MODULE multiple_scattering_routines
    

    implicit none
    
    CONTAINS
    
!***********************************************************************************
!*
!*      get_mul_corr.FOR
!*
!*      A K Soper, December 1999
!*
!*      gets multiple scattering correction for any run
!*
!***********************************************************************************
    subroutine get_mul_corr(ntype,run,nperrq,forcecalc)

    use run_par
    use reallocation_routines
    use inputfilestrings
    use calibration_routines
    use beam_routines
    use mul_corr_azi
    use van_mul_azi
    use sam_mul_azi
    use corrections_routines
!c
!c internal variables
!c
      character(len=256)                          :: fname            !dummy file name
      integer i,ib,im,j,il,ind,k,ierr      !internal indices
      character(len=256) run                  !filename
      integer nperrq            !period number of data in raw file
      integer ntype                  !type of file being read
      integer forcecalc            != 1 forces calculation of corrs.
      real                                        :: aa                  !temporary real values
      real, dimension(:), allocatable             :: bb,cc,dd!(mcorrang)                  !temporary real values
        
        ierr=1 !Signal that corrections need to be calculated in case corrections file has problems
        if(forcecalc.eq.0) then
!c
!c set up filename of file to be read
!c
            fname=run
            call change_ext(fname,newext='mul')
            write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
!c
!c try to open multiple scattering correction file. If an error occurs programme
!c redirects to calculate the multiple scattering correction 
!c
            open(10,file=fname,status='old',iostat=ierr)
            if(ierr.eq.0) then
                read(10,*,iostat=ierr) nwavmul,nangmul,nazimul
                write(6,*) nwavmul,nangmul,nazimul
            end if
            if(ierr.eq.0) then
                !allocate the required space
                call reallocate1d_r(wavmul,nwavmul)
                call reallocate1d_r(angmul,nangmul)
                call reallocate1d_r(azimul,nazimul)
                call reallocate1d_r(bb,nazimul)
                call reallocate1d_r(cc,nazimul)
                call reallocate3d_r(onescat,nwavmul,nangmul,nazimul)
                call reallocate3d_r(mulscat,nwavmul,nangmul,nazimul)
                read(10,*,iostat=ierr) (azimul(k),k=1,nazimul)
            endif
            im=0
            do while (im.lt.nangmul.and.ierr.eq.0)
                im=im+1
                read(10,*,iostat=ierr) ib,angmul(im)
            j=0
                do while (j.lt.nwavmul.and.ierr.eq.0)
                    j=j+1
                    read(10,*,iostat=ierr) aa,(bb(k),cc(k),k=1,nazimul)
                    if(ierr.eq.0) then
                        wavmul(j)=aa
                  do k=1,nazimul
                     onescat(j,im,k)=bb(k)
                     mulscat(j,im,k)=cc(k)
                  end do
                    end if
            end do
            end do
            close(10)
        end if
        if(ierr.ne.0) then
!c .MUL file does not exist, or has a problem. Calculate multiple scattering
!c correction from input sample parameters
            call calc_mul(ntype,run,nperrq)
        end if
!c
!c now save the multiple scattering correction in the appropriate arrays
!c
      if(ntype.eq.1) then
!Set up the required arrays
            nvwavmul=nwavmul
            nvangmul=nangmul
            nvazimul=nazimul
            call reallocate1d_r(vwavmul,nwavmul)
            call reallocate1d_r(vangmul,nangmul)
            call reallocate1d_r(vazimul,nazimul)
            call reallocate3d_r(vonescat,nwavmul,nangmul,nazimul)
            call reallocate3d_r(vmulscat,nwavmul,nangmul,nazimul)
            do il=1,nwavmul
            vwavmul(il)=wavmul(il)
            end do
            do ib=1,nangmul
                vangmul(ib)=angmul(ib)
                do il=1,nwavmul
                    do k=1,nazimul
                        vonescat(il,ib,k)=onescat(il,ib,k)
                        vmulscat(il,ib,k)=mulscat(il,ib,k)
                    end do
                end do
            end do
            do k=1,nazimul
            vazimul(k)=azimul(k)
            end do
      else
            ind=ntype-3
            nswavmul=nwavmul
            nsangmul=nangmul
            nsazimul=nazimul
            call reallocate1d_r(swavmul,nwavmul)
            call reallocate1d_r(sangmul,nangmul)
            call reallocate1d_r(sazimul,nazimul)
            call reallocate4d_r(sonescat,nwavmul,nangmul,nazimul,ind)
            call reallocate4d_r(smulscat,nwavmul,nangmul,nazimul,ind)
            do il=1,nwavmul
            swavmul(il)=wavmul(il)
            end do
            do ib=1,nangmul
            sangmul(ib)=angmul(ib)
            do il=1,nwavmul
                    do k=1,nazimul
                        sonescat(il,ib,k,ind)=onescat(il,ib,k)
                  smulscat(il,ib,k,ind)=mulscat(il,ib,k)
                    end do
                end do
            end do
            do k=1,nazimul
            sazimul(k)=azimul(k)
            end do
      endif
        if(allocated(wavmul)) deallocate(wavmul,angmul,azimul,onescat,mulscat)
      return
    end subroutine get_mul_corr
    
!***********************************************************************************
!*
!*      calc_mul.FOR
!*
!*      A K Soper, December 1999
!*
!*      calculates multiple scattering correction for any run
!*
!***********************************************************************************
    subroutine calc_mul(ntype,run,nperrq)
        
        use run_par
        use reallocation_routines
        use inputfilestrings
        use calibration_routines
        use beam_routines
        use mul_corr_azi
        use van_mul_azi
        use sam_mul_azi
        use corrections_routines
        use van_par
        use sam_par

!c
!c internal variables
!c
      character(len=256)                          :: fname            !file name for m.s.corrections
      character(len=32)                           :: formatout
      integer                                     :: i,ib,ic,j,ind,il,k      !internal indices
      character(len=256)                          :: run                  !filename
      integer                                     :: nperrq            !period number of data in raw file
      integer                                     :: ntype                  !type of file being calculated for
      integer                                     :: nblock            !number of blocks in output file
      integer                                     :: lrec            !record length 
      real                                        :: aa
!        real, dimension(:), allocatable             :: bb,cc,dd!(mcorrang)      !temporary values
      real                                        :: pi,pi4

        write(6,*) a,b,a1,b1
      write(6,*) hbdown,hbup,hsbdown,hsbup

        pi=4.0*atan(1.0)
        pi4=4.0*pi
!c
!c calculate the appropriate corrections
!c
      if(ntype.eq.1) then
            nan=1
            call reallocate1d_r(rad1,nan)
            call reallocate1d_r(rad2,nan)
            call reallocate1d_r(den,nan)
            call reallocate1d_r(captcs,nan)
            rad1(1)=vdimen1
            rad2(1)=vdimen2
            theight=vheight
            den(1)=vrho
            captcs(1)=vcaptav
            nwavmul=nlambv
            call reallocate1d_r(wavmul,nwavmul)
            call reallocate2d_r(sigtl,nwavmul,nan)
            call reallocate2d_r(sigsl,nwavmul,nan)
            do i=1,nwavmul
                wavmul(i)=vlamb(i)
                sigtl(i,1)=vtscat(i)
                sigsl(i,1)=vsscat(i)
!                write(6,*) 'calc_mul> ',wavmul(i),sigtl(i,1)
            end do
            ncval=0
            write(6,*) 'calc_mul> ',nan,rad1(nan),rad2(nan),den(nan),captcs(nan),sigtl(1,nan),theight
!c
!c set up angles for corrections
!c
            call set_corr_ang(vgeom,vrot,angmul,azimul,nangmul,nazimul)
            write(6,*) 'calc_mul> Set up correction angles'
            if(vgeom.eq.1) then
                call cylmultof()
                write(6,*) 'calc_mul> Calculated cylindrical corrections'
            else if(vgeom.eq.2) then
                call fltmultof(vrot,vwidth,nan)
                write(6,*) 'calc_mul> Calculated flat plate corrections'
            endif
      else if(ntype.gt.3) then
            ind=ntype-3
            nwavmul=nlambs
            call reallocate1d_r(wavmul,nwavmul)
            do i=1,nwavmul
                wavmul(i)=slamb(i)
            end do
            theight=sheight
            !Allocate the required number of arrays
            nan=ncont-ind+1
            write(6,*) nan
            call reallocate1d_r(rad1,nan)
            call reallocate1d_r(rad2,nan)
            call reallocate1d_r(den,nan)
            call reallocate1d_r(captcs,nan)
            call reallocate2d_r(sigtl,nwavmul,nan)
            call reallocate2d_r(sigsl,nwavmul,nan)
            ncval=0
            !Use the largest value of ncval for the set of sample plus containers
            do i=1,ncont
                ncval=max(ncval,ncvals(i))
            end do
            write(6,*) ncval
            if(ncval.gt.0) then
                call reallocate1d_r(czval,ncval)
                call reallocate3d_r(pc,ncval,nwavmul,nan)
            end if
            nan=0
            do i=ind,ncont
            nan=nan+1
                rad1(nan)=sdimen1(i)
                rad2(nan)=sdimen2(i)
                den(nan)=srho(i)
!c
!c for sample include tweak factor in sample density
!c
                if(i.eq.1) den(nan)=den(nan)/abs(tweak(i))
                captcs(nan)=scaptav(i)
                do il=1,nwavmul
                    sigtl(il,nan)=stscat(il,i)
                    sigsl(il,nan)=ssscat(il,i)
                    if(ncval.gt.0) then
!c If a mint file has been named for this sample or container grab the corresponding DCS data as a function of cos theta and wavelenght
                        if(index(smintname(i),'.mint').gt.0) then
                            do ic=1,ncval
                                pc(ic,il,nan)=pcs(ic,il,i)
                            end do
!c If ncval is > 0, but mint is not specified we still need to set the values of the DCS for each data file
                        else
                            do ic=1,ncval
                                pc(ic,il,nan)=(ssscat(il,i))/pi4
                            end do
                        endif
                    end if
                end do
                write(6,*) 'calc_mul> ',nan,rad1(nan),rad2(nan),den(nan),captcs(nan),sigtl(1,nan),theight
            end do
!c
!c set up angles for corrections
!c
            call set_corr_ang(sgeom,srot(1),angmul,azimul,nangmul,nazimul)
            write(6,*) 'calc_mul> Set up correction angles'
            if(sgeom.eq.1) then
            call cylmultof()
                write(6,*) 'calc_mul> Calculated cylindrical corrections'
            else if(sgeom.eq.2) then
            call fltmultof(srot(1),swidth(1),ind)
                write(6,*) 'calc_mul> Calculated flat plate corrections'
            endif
      endif
!c
!c set up filename of file to write
!c
        fname=run
        call change_ext(fname,newext='mul')
        write(fname,'(a,i2.2)') fname(1:len_trim(fname)),nperrq
!c
!c set the output format
!c
      nblock=2*nazimul+1
      write(formatout,332) '(',nblock,'(1x,e12.5))'
332      format(a1,i3.3,a11)
      write(6,*) formatout
      lrec=nblock*13
      if(lrec.lt.80) lrec=80
!        call reallocate1d_r(bb,nazimul)
!        call reallocate1d_r(cc,nazimul)
!c
!c Open the multiple scattering correction file and write it. 
!c
      open(10,file=fname,status='unknown',form='formatted',recl=lrec)
      write(10,*) nwavmul,nangmul,nazimul
      write(10,formatout) (azimul(k),k=1,nazimul)
      do ib=1,nangmul
            write(10,*) ib,angmul(ib)
            do j=1,nwavmul
!            aa=wavmul(j)
!            do k=1,nazimul
!                    bb(k)=onescat(j,ib,k)
!                    cc(k)=mulscat(j,ib,k)
!            end do
            write(10,formatout) wavmul(j),(onescat(j,ib,k),mulscat(j,ib,k),k=1,nazimul)
            end do
      end do
      close(10)
      return
    end subroutine calc_mul
    
    subroutine CYLMULTOF()
        
        use reallocation_routines
        use run_par
        use calibration_routines
        use beam_routines
        use mul_corr_azi
        use van_mul_azi
        use sam_mul_azi
        use corrections_routines
        use cylmultof_routines
        
        integer nw,i,il,k1,k2,k3,kk
        real astep,h2,up,down,fac
!c
!c      Original   from AKS
!c
!C PROGRAM MULSCA - CONTROL PROGRAM FOR THE Cylindrical MULTIPLE SCATTERING
!C CALCULATION 
!C               
!c Modified March 2001 into a subroutine to be linked into GUDRUN
!c
!c      modified 5/02/2003 to include situation with azimuthal detectors
!c
!C                                                                     
!C ARRAYS - PRIMARY SCATTERING:-                                       
!C        SUMC = C-C                                                     
!C        SUMS = C-P                                                    
!C        SUMT = P-C                                                     
!C        SUMP = P-P                                                    
!C                                                                      
!C  ARRAYS - SECONDARY SCATTERING:-                                      
!C                                                                      
!C        SUMA = C-C                                                     
!C        SUMD = C-P                                                    
!C        SUME = P-C                                                    
!C        SUMB = P-P                                                    
!C                                                                      
!c output unit for diagnostic messages
!c
      nw=6
      write(nw,*) 'Got to CylMultof'
!C                                                                      
!C SET NZ IN THE RATIO OF HEIGHT/ASTEP                                  
!C                                                                       
      astep=stepm
      HIGHT=theight
!c
!c set beam height parameters with reference to bottom of sample rather
!c than centre: this is to correspond to definitions used in cylmultof
!c
      h2=0.5*HIGHT
      hdown=hbdown+h2
      hup=hbup+h2
      hsdown=hsbdown+h2
      hsup=hsbup+h2
        NZ=(HIGHT+0.5*ASTEP)/ASTEP
!        write(6,*) 'cylmultof> nz,HIGHT,astep: ',nz,HIGHT,astep
!C                                                                       MUL01900
!C DEFINE INTEGERS CORRESPONDING AS CLOSE AS POSSIBLE                    MUL01910
!C TO THE HEIGHT OF THE INCIDENT BEAM                                    MUL01920
!C                                                                       MUL01930
      IF(HUP.GT.HIGHT) HUP=HIGHT
        IF(HDOWN.LT.0.0) HDOWN=0.0
        HIGHB=HUP-HDOWN
        ZSTEP=HIGHT/NZ
        MINB=(HDOWN+0.5*ZSTEP)/ZSTEP+1
        MAXB=1+(HUP-0.5*ZSTEP)/ZSTEP
!C                                                                       MUL02000
!C                                                                       MUL02010
!C DEFINE INTEGERS CORRESPONDING AS CLOSE AS POSSIBLE TO THE             MUL02020
!C HEIGHT OF THE SCATTERED BEAM                                          MUL02030
!C                                                                       MUL02040
        IF(HSUP.GT.HIGHT) HSUP = HIGHT
        IF(HSDOWN.LT.0.) HSDOWN=0.
        HIGHS=HSUP-HSDOWN
        MINS=(HSDOWN+0.5*ZSTEP)/ZSTEP+1
        MAXS=(HSUP-0.5*ZSTEP)/ZSTEP+1
!C                                                                       MUL02100
!C DETERMINE HEIGHT OF OVERLAP BETWEEN INCIDENT AND SCATTERED BEAMS      MUL02110
!C                                                                       MUL02120
        UP=min(HUP,HSUP)
        DOWN=max(HDOWN,HSDOWN)
        HIGHBS=UP-DOWN
        IF(HIGHBS.LT.0.) HIGHBS=0.
!          write(6,*) 'cylmultof> a,b,a1,b1',a,b,a1,b1
      write(6,*) 'cylmultof> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
      write(6,*) 'cylmultof> HIGHB,HIGHS,HIGHBS,HIGHT ',HIGHB,HIGHS,HIGHBS,HIGHT
      write(6,*) 'cylmultof> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
!C
!C CALCULATE  SCATTERING C/S AT EACH WAVELENGTH ASSUMING 1/LAMBDA ABSORPTION
!C
        call reallocate1d_r(sigs,nwavmul)
        call reallocate1d_r(sigt,nwavmul)
        call reallocate1d_r(mus,nwavmul)
        call reallocate3d_r(onescat,nwavmul,nangmul,nazimul)
        call reallocate3d_r(mulscat,nwavmul,nangmul,nazimul)
        DO IL=1,nwavmul
            do i=1,nan
                  SIGS(I)=SIGSL(IL,I)
                  SIGT(I)=SIGTL(IL,I)
                  MUS(I)=SIGT(I)*DEN(I)
            end do
!            write(6,*) 'cylmultof> wavmul,(sigs,sigt,mus) ',wavmul(il),(sigs(i),sigt(i),mus(i),i=1,nan)
!C
!C SET UP THE GEOMETRY OF THE PROBLEM AND ARRAYS ALEN1,ALEN2             MUL00480
!C                                                                       MUL00490
!            write(6,*) 'Going to LEN3'
!      write(6,*) 'cylmultof> HIGHB,HIGHS,HIGHBS,HIGHT ',HIGHB,HIGHS,HIGHBS,HIGHT
!      write(6,*) 'cylmultof> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
            CALL LEN3(nw,ASTEP)      
            MS=NRAD(1)
!            write(6,*) 'Going to LEN4'
!      write(6,*) 'cylmultof> HIGHB,HIGHS,HIGHBS,HIGHT ',HIGHB,HIGHS,HIGHBS,HIGHT
!      write(6,*) 'cylmultof> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
            CALL LEN4(nw)!,OM,ALEN1,NTEST,NBLEN1)
!            write(6,*) 'Going to LEN5'
!      write(6,*) 'cylmultof> HIGHB,HIGHS,HIGHBS,HIGHT ',HIGHB,HIGHS,HIGHBS,HIGHT
!      write(6,*) 'cylmultof> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
            CALL LEN5(nw)!,OM,ALEN2,NT2,NBLEN1,NBLEN2)
!C                                                                       MUL00530
!C PERFORM THE PRIMARY AND MULTIPLE SCATTERING CALCULATION               MUL00540
!C                                                                       MUL00550
!            write(6,*) 'Going to SUMMU1'
!      write(6,*) 'cylmultof> HIGHB,HIGHS,HIGHBS,HIGHT ',HIGHB,HIGHS,HIGHBS,HIGHT
!      write(6,*) 'cylmultof> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
            CALL SUMMU1()
      write(6,*) (SUMP(1,kk,1),kk=1,nazimul)
      write(6,*) (sumB(1,kk,1,1),kk=1,nazimul)
!C                                                                       MUL00580
!C OUTPUT THE RESULTS                                                    MUL00590
!C                                                                       MUL00600
!            write(6,*) 'Going to DATAOP'
            CALL DATAOP(IL)
      end do
        return
    end subroutine cylmultof
    
    subroutine fltmultof(rot1,width,ind)
        
        use reallocation_routines
        use run_par
        use beam_routines
        use mul_corr_azi
        use van_mul_azi
        use sam_mul_azi
        use corrections_routines
        use fltmultof_routines
!c
!c      original   by aks
!c      modified   by wsh   1-3-89    for coral input
!c 14-5-90  : modified by aks and ach to increase size of array thetab for
!c            detector angles from 8 to 200 since lad requires 14 for a flat
!c            plate with standard grouping.
!c 16-7-90  : for secondary angles near 90 degrees the program puts
!c out unity for the primary scattering and zero for the multiple scattering
!c
!c completely revamped by aks 17-05-01 to become a subroutine of the gudrun
!c suite
!c
!c      modified 5/02/2003 to include situation with azimuthal detectors
!c
!c Modified again on 10/12/2010 to allow for explicit cos theta dependence of scattering
!c cross section. This applies to the second scattering - higher orders are treated within
!c the isotropic approximation as previously.

      integer                                         :: i,ind,it,i1,iref,ic,ic1,icref,ib,il,ip,ierr
        integer                                         :: iz,iz1,iz2,izref,j,k,nr,nz,nc,nslice2
      real                                            :: rot,rot1            !angle of flat plate normal to beam
      real                                            :: width,width2      !lateral width and half width of sample
      real, dimension(:), allocatable                 :: amus,amut!(mcorrwav)      !average scattering and total attenuation coeffsicients
      real                                            :: amuabs            !average capture cross section
      real                                            :: tsum                  !total thickness of slab
      real                                            :: theta,thetad      !temporary angle values
      real                                            :: pi,piconv      
      real                                            :: tsec,sec1,sec2      !secants of angles
      real, dimension(:), allocatable                 :: ts!(mcont)            !wall thickness of each container
      real                                            :: leakfac            !leakage factor to correct for finite width
        real                                            :: aa,bb
        real                                            :: accur,al0,al1,al2,amu,amusct,amuxs,analytic,ctheta
        real                                            :: crot1,delphi,delr,delz,fac,gamma,sum
        real                                            :: llimit,mutsec1,phi,pi2,pi4,r,requiv,rest,srot1
        real                                            :: sumdcsphi,sumphi,sumr,sumz,t,tes,unit,volel
        real                                            :: x,x1,xh,xx,z,z1,z2,zval,phid,psec,thetar,phir
        real, dimension(:,:),allocatable                :: xd,yd,zd!(mcorrang,mcorrang),
        real                                            :: xe,ye,ze,term,termd
        real                                            :: cthetar,sthetar,sthetae,cthetae,expfac
        real, dimension(:,:), allocatable               :: first,sumtot,sumrdet,sumdcsphidet,sumzdet!(mcorrang,mcorrang)
        real, dimension(:,:,:), allocatable             :: eout,sumnx2nd,aintendet!(mslice,mcorrang,mcorrang)
!c
!c internal arrays used by the integration routines
!c
      real, dimension(:), allocatable                 :: ainten,sumint,sumnx!(mslice)
        real, dimension(:), allocatable                 :: aintenc!(2*mslice+1)
        real, dimension(:,:), allocatable               :: pct!Temporary storage for differential cross section
        real                                            :: l,lsq,mutl

        character(len=256)                              :: fnameout
        logical                                         :: test

        nslice2=2*nslice+1
        call reallocate1d_r(amus,nwavmul)
        call reallocate1d_r(amut,nwavmul)
        call reallocate1d_r(ts,nan)
        call reallocate1d_r(ainten,nslice)
        call reallocate1d_r(sumint,nslice)
        call reallocate1d_r(sumnx,nslice)
        call reallocate2d_r(first,nangmul,nazimul)
        call reallocate2d_r(sumtot,nangmul,nazimul)
        call reallocate2d_r(sumrdet,nangmul,nazimul)
        call reallocate2d_r(sumdcsphidet,nangmul,nazimul)
        call reallocate2d_r(sumzdet,nangmul,nazimul)
        call reallocate3d_r(eout,nslice,nangmul,nazimul)
        call reallocate3d_r(sumnx2nd,nslice,nangmul,nazimul)
        call reallocate3d_r(aintendet,nslice2,nangmul,nazimul)
        call reallocate3d_r(onescat,nwavmul,nangmul,nazimul)
        call reallocate3d_r(mulscat,nwavmul,nangmul,nazimul)
        if(ncval.gt.0) then
            call reallocate2d_r(pct,ncval,nwavmul)
            call reallocate1d_r(aintenc,nslice2)
        end if
            
        pi=4.0*atan(1.0)
      piconv=pi/180.
      pi4=4.*pi
      pi2=2.*pi
      rot=rot1*piconv
        crot1=cos(rot)
        srot1=sin(rot)
        sec1=1./crot1
!c Correct for round off errors?
        if(sec1.lt.1..and.sec1.gt.0.) sec1=1.
        if(sec1.gt.-1..and.sec1.lt.0.) sec1=-1.
      gamma=0.5772156649
      width2=0.5*width

!c
!c specified accuracy of integrals
!c
      accur=0.0001
!C
!C form average total and capture cross sections for sample and containers
!C
      tsum=0.0
      do ic=1,nan
            ts(ic)=rad2(ic)+rad1(ic)
            tsum=tsum+ts(ic)
      end do
!c
!c total thickness of slab
!c
      t=tsum
      write(6,*)
      write(6,'(1x,a,f10.5,1x,a,i5)') 'fltmultof> Total thickness (cm): ',t,'No. of cos theta values: ',ncval

!c
!c define thickness of step for integration
!c
      x=t/real(nslice)
      do il=1,nwavmul
            amut(il)=0.0
            amus(il)=0.0
            amuabs=0.0
!c
!c Zero the DCS values for this sequence of samples and/or containers
!c
            if(ncval.gt.0) then
                do iz=1,ncval
                    pct(iz,il)=0.0
                end do
            endif
            do ic=1,nan
                amut(il)=amut(il)+ts(ic)*den(ic)*sigtl(il,ic)/tsum
                amus(il)=amus(il)+ts(ic)*den(ic)*sigsl(il,ic)/tsum
                amuabs=amuabs+ts(ic)*den(ic)*captcs(ic)/tsum
!c
!c Contribution to DCS if defined
!c
                if(ncval.gt.0) then
                    do iz=1,ncval
                        pct(iz,il)=pct(iz,il)+ts(ic)*den(ic)*pc(iz,il,ic)/tsum
                    end do
                endif
            end do
        end do
        
!        do il=1,nwavmul
!            fac=wavmul(il)*amuabs/1.7979
!            write(6,*) 'fltmultof> ',il,wavmul(il),amut(il),fac,amuabs
!            amus(il)=amut(il)-fac
!            write(6,*) 'fltmultof> ',il,wavmul(il),amus(il),amut(il),fac,amuabs
!        end do

!c Write out pc to check it was generated correctly

        if(ncval.gt.0) then
            do it=1,nangmul
                thetar=angmul(it)*piconv
                cthetar=cos(thetar)
                sthetar=sin(thetar)
                do ip=1,nazimul
                    phir=azimul(ip)*piconv
                    xd(it,ip)=sthetar*cos(phir)
                    yd(it,ip)=sthetar*sin(phir)
                    zd(it,ip)=cthetar
                end do
            end do
            write(fnameout,'(a,i1,a)') 'fltmultof',ind,'.pofc'
            open(10,file=fnameout,status='unknown',iostat=ierr)
            if(ierr.eq.0) then
                do il=1,nwavmul
                    write(10,'(/a,1x,i5,1x,i5,1x,f10.5)') '#',ncval,nwavmul,wavmul(il)
                    do iz=1,ncval
                        write(10,'(f10.5,1x,i5,1x,f10.5,1x,e12.5)') wavmul(il),iz,czval(iz),pct(iz,il)
                    end do
                end do
            endif
            close(10)
        endif
!c
!c step through wavelengths
!c
        do il=1,nwavmul
!c
!c calculate abs c/s
!c
            amu=amut(il)
            amusct=amus(il)
            unit=amusct/pi4
            fac=amusct*t/pi4
!            write(6,*) 'fltmultof> 1 ',il,wavmul(il),amu,amusct,unit,fac
            if(amu.lt.0.0.or.amusct.lt.0.0) then
                write(6,*) 'fltmultof> Cannot process negative cross sections!'
                stop
            end if
!c
!c calculate leakage factor for neutrons out the side of the slab
!c
            leakfac=width2*amu
            if(leakfac.lt.20.0) then
                leakfac=1.0-exp(-leakfac)
            else
                    leakfac=1.0
            endif
!c
!c calculate scattering from a point to each plane
!c
            xx=0.
            ainten(1)=pi2*unit*binte1(amu,xx,x)
            do i=2,nslice
            x1=(i-1)*x
            ainten(i)=unit*pi2*binte1(amu,x1,x)
            end do
!           write(6,*) 'fltmultof> 2 ',il,wavmul(il),amu,amusct,unit,fac
!c
!c calculate primary intensity at each plane
!c
            xh=x*0.5
!c 
!c If a DCS file has been defined for one or more of the containers, we need to calculate an alternative version of
!c sumnx which includes the dcs for scattering from the first point to the scattering element.
!c
            if(ncval.gt.0) then
                nc=ncval/2
                mutsec1=amu*sec1
                zval=t
                nz=nslice
                delz=zval/nz
                llimit=12.0
                nr=100

!c Determine radial step size based on llimit

                delr=llimit/(amu*nr)
!c Approximate volume of integration element
                volel=delz*delr*delr
!c Calculate the radius of a sphere with the same volume as an element
                requiv=(3*volel/(4*pi))**(1.0/3.0)

                iz=0
!c NB for the isotropic calculation -z gives identical results to +z, but when the DCS is included, this is no longer true
!c so we need to keep +z separate from -z
                z=-zval-delz
                do while (iz.lt.2*nz+1)

                    iz=iz+1
                    z=z+delz
                    sumr=0.0
                    do it=1,nangmul
                        do ip=1,nazimul
                            sumrdet(it,ip)=0.0
                        end do
                    end do
                    r=-0.5*delr
                    l=0
                    mutl=amu*l
                    do while (mutl.lt.llimit)

                        r=r+delr
                        lsq=z*z+r*r
                        l=sqrt(lsq)
                        mutl=amu*l
                        sumphi=0.0
                        sumdcsphi=0.0
                        do it=1,nangmul
                            do ip=1,nazimul
                                sumdcsphidet(it,ip)=0.0
                            end do
                        end do
!c
!c Sine and Cosine (theta) of current element
!c
                        cthetae=z/l
                        cthetae=min(cthetae,1.0)
                        cthetae=max(cthetae,-1.0)
                        sthetae=sqrt(1-sthetae*sthetae)
                        delphi=delr/r
                        phi=-0.5*delphi
                        do while (phi.lt.pi2)
                            phi=phi+delphi
!c Direction cosines of current element
                            xe=sthetae*cos(phi)
                            ye=sthetae*sin(phi)
                            ze=cthetae
!c
!c Get the rotated z coordinate and use to calculate the new cos(scattering angle).
!c Note that in this case we are rotate the slab relative to the coordinate system.
!c This is the opposite to the case in get_corr_ang where we rotate the coordinates of
!c the detector relative to the slab axes. This because the point defined by ctheta,phi is in the slab itself -
!c it is not a detector fixed in the lab coordinate system. As a result the sign of the rotation is
!c opposite in the two cases. As elsewhere the rotation is about the laboratory (or slab) y-axis.
!c
                            ctheta=(ze*crot1-xe*srot1)
                            icref=int((1.0+ctheta)*nc)+1
                            icref=min(icref,ncval)
                            term=delphi*pct(icref,il)
                            sumdcsphi=sumdcsphi+term
!c 
!c Calculate the contribution for scattering to each detector
!c
                            do it=1,nangmul
                                do ip=1,nazimul
!c Cosine of scattering angle to detector
                                    ctheta=xe*xd(it,ip)+ye*yd(it,ip)+ze*zd(it,ip)
                                    ctheta=min(ctheta,1.0)
                                    ctheta=max(ctheta,-1.0)
                                    icref=int((1.0+ctheta)*nc)+1
                                    icref=min(icref,ncval)
                                    termd=term*pct(icref,il)
                                    sumdcsphidet(it,ip)=sumdcsphidet(it,ip)+termd
                                end do
                            end do
                            sumphi=sumphi+delphi
                        end do
!c Since the total phi covered might be more or less than 2*pi exactly, we renormalise the sums accordingly.
                        term=pi2/sumphi
                        sumdcsphi=term*sumdcsphi
                        do it=1,nangmul
                            do ip=1,nazimul
                                sumdcsphidet(it,ip)=term*sumdcsphidet(it,ip)
                            end do
                        end do
                        term=r*avl2(requiv,requiv,l)*exp(-mutl)/lsq
                        sumr=sumr+sumdcsphi*term
                        do it=1,nangmul
                            do ip=1,nazimul
                                sumrdet(it,ip)=sumrdet(it,ip)+sumdcsphidet(it,ip)*term
                            end do
                        end do
                    end do
                    aintenc(iz)=sumr*delr
                    do it=1,nangmul
                        do ip=1,nazimul
                            aintendet(iz,it,ip)=sumrdet(it,ip)*delr
                        end do
                    end do
!c Calculate analytic version of this integral
                    analytic=unit*pi2*binte1(amu,z,delz)     
                end do
!c            write(6,*) 'fltmultof3> ',il,wavmul(il),z,amusct,amu,aintenc(iz)
!c     +,aintendet(iz,1,1),analytic
!c
!c At this stage aintenc represents the intensity at each slice due to scattering from a
!c point a perpendicular distance z away. Now integrate this for the incident beam to give the second intensity
!c at each z value. Compare with analytical calculation
!c
                z1=-0.5*delz
                iz1=0
                do while (iz1.lt.nz)
                    iz1=iz1+1
                    z1=z1+delz
                    sumz=0.0
                    do it=1,nangmul
                        do ip=1,nazimul
                            sumzdet(it,ip)=0.0
                        end do
                    end do
                    z2=-0.5*delz
                    iz2=0
                    do while(iz2.lt.nz)

                        iz2=iz2+1
                        z2=z2+delz
                        izref=iz2-iz1+nz+1
                        expfac=exp(-z2*mutsec1)
                        sumz=sumz+aintenc(izref)*expfac
                        do it=1,nangmul
                            do ip=1,nazimul
                                sumzdet(it,ip)=sumzdet(it,ip)+aintendet(izref,it,ip)*expfac
                            end do
                        end do
                    end do
                    term=sec1*delz
                    sumnx(iz1)=term*sumz*amusct
                    term=term*pi4
                    do it=1,nangmul
                        do ip=1,nazimul
                            sumnx2nd(iz1,it,ip)=term*sumzdet(it,ip)
                        end do
                    end do
!c Compare this with the analytical result

                    al1=amu*z1
                    al2=amu*zval-al1
                    amuxs=al1*sec1
                    aa=e1int1(sec1,al1)
                    bb=e1int2(sec1,al2)
                    analytic=(aa+bb*exp(-amuxs))*0.5*amusct**2/amu

                end do
!c            write(6,*) wavmul(il),z1,sumnx(iz1),sumnx2nd(iz1,1,1),analytic
            else
!c Otherwise we simply assume isotropic first scattering
                do i=1,nslice
                    x1=(i-1)*x+xh
                    al1=amu*x1
                    al2=amu*t-al1
                    amuxs=al1*sec1
                    aa=e1int1(sec1,al1)
                    bb=e1int2(sec1,al2)
                    sumnx(i)=(aa+bb*exp(-amuxs))*0.5*amusct**2/amu
                end do
            end if
!            write(6,*) 'fltmultof> 3 ',il,wavmul(il),amu,amusct,unit,fac
!c      write(6,*) 'fltmultof4> ',il,wavmul(il),amu,amusct,unit,fac
!c      write(6,*) (sumnx(i),i=1,nslice)
!c
!c step through banks and calculate corrections at each wavelength and
!c scattering angle
!c
            do it=1,nangmul
                do ip=1,nazimul
                    thetad=angmul(it)
                    phid=azimul(ip)
!c
!c tsec is the angle the scattered beam makes with the normal to the sample
!c surface.  if abs(tsec) is close to 90 deg. calculation of multiple scattering
!c coefficients is unreliable
!c
                    if(abs(abs(thetad)-90.).gt.1.) then
                        tsec=thetad*piconv
                        psec=phid*piconv
                        sec2=1./cos(tsec)
                        tes=t
                        if(sec2.lt.0.) tes=0.
!c
!c Calculate the attenuation factor to the detector
!c
                        xh=x*0.5
                        do i=1,nslice
                            x1=(i-1)*x+xh
                            al1=(tes-x1)
                            eout(i,it,ip)=aintex(al1,sec2,amu,x)/pi4
                        end do
!c Include the scattering angle in the primary scattering if the scattered dcs is defined
!c 23/12/2010 In practice this does not work very well since the single scattering may be used to
!c estimate the sample dependent background, which can come from different scattering angles,
!c so better to use simply the total scattering cross section at each energy.
                        first(it,ip)=prime(sec1,sec2,amu,t)*fac*sec1
                    else
                        do i=1,nslice
                            eout(i,it,ip)=0.0
                        end do
                        first(it,ip)=1.0
                    endif
                    sumtot(it,ip)=0.0
!                   write(6,*) 'fltmultof> 4 ',il,it,ip,wavmul(il),first(it,ip)
                end do
            end do
!c
!c calculate nord orders of scattering
!c
            i=1
            test=.true.
            do while(test)
!c
!c Iterate over the scattering angles
!c
                test=.false.
                do it=1,nangmul
                    do ip=1,nazimul
                        sum=0.
                        if(i.eq.1) then
                            sum=first(it,ip)
                        else if(i.eq.2.and.ncval.gt.0) then
                            do j=1,nslice
                                sum=sum+sumnx2nd(j,it,ip)*x*leakfac*eout(j,it,ip)
                            end do
                        else
                            do j=1,nslice
                                sum=sum+sumint(j)*eout(j,it,ip)
                            end do
                        endif
                        sumtot(it,ip)=sumtot(it,ip)+sum
                        test=test.or.sum.ge.accur*sumtot(it,ip)
!           write(6,*) 'fltmultof> 5 ',il,it,ip,wavmul(il),first(it,ip),sumtot(it,ip)
                    end do
                end do
!c      write(6,*) 'fltmultof6> ',il,wavmul(il),it,ip,i,sum,sumtot(nangmul,1)
                if(test) then
                    if(i.gt.1) then
!c calculate next order
                        do j=1,nslice
                            sumnx(j)=0.0
                            do k=1,nslice
                                iref=iabs(j-k)+1
                                sumnx(j)=sumnx(j)+ainten(iref)*sumint(k)
                            end do
                        end do
                    endif
!c Correct for leakage
                    do j=1,nslice
                        sumint(j)=sumnx(j)*x*leakfac
                    end do
                endif
                i=i+1
            end do
!c
!c save the results
!c
            do it=1,nangmul
                do ip=1,nazimul
                    onescat(il,it,ip)=first(it,ip)
                    mulscat(il,it,ip)=sumtot(it,ip)-first(it,ip)
                end do
            end do
        end do
      return
    
    end subroutine fltmultof
    
END MODULE multiple_scattering_routines
