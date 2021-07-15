      program calc_corrsx

      include 'dimension.inc'
      include 'beam.inc'
       include 'inputfilestrings.inc'
      include 'abs_corr_azi.inc'
      include 'mul_corr_azi.inc'
      include 'sam_par.inc'
      character*256 text

c internal variables

      character*32 formatout
      integer i,ib,j,ind,il,iang,k      !internal indices
      character*256 run,filenamebase,fname                  !run number and filenames
      integer nperrq            !period number of data in raw file
      integer nan                  !temporary array of number of annuli
      integer ncorr                  !no. of absorption corrections per lamb.
      integer nblock            !number of blocks in output file
      integer lrec            !record length
      integer ilenfilenamebase    !character length of filenamebase
      real rad1(mcont),rad2(mcont)      !temporary dimension values
      real den(mcont)            !density of each sample and cont.
      real captcs(mcont)                  !capture cross section of sample or cont
      real captcs_out(mcont)                  !capture cross section of sample or cont
      real sigtl(mcorrwav,mcont)      !total cross section of sample or cont.
      real sigtl_out(mcorrwav,mcont)      !total cross section of sample or cont.
      real abstemp(mcorrwav,mcorrang,mcorrang,mcont) !temporary attenuation correction
      real onecor(mcorrwav,mcorrang,mcorrang,mcont) !array to store single scattering for each set of annuli
      real mulcor(mcorrwav,mcorrang,mcorrang,mcont) !array to store multiple scattering for each set of annuli
      real azidet,afirst,alast
      logical equalinout   !Signifies whether the incident and output cross sections are the same
      logical ignoresample  !For flat plate samples, tells cylabstof not to include the sample when calculating As,sc.

c Get the beam parameters

      call getBeamParameters()

c set up some initial arrays

      nwavabs=1
      nwavmul=nwavabs
      do i=1,nwavabs
            wavabs(i)=1.0
            wavmul(i)=wavabs(i)
c      write(6,*) i,slamb(i)
      end do

c Get input from any supplied arguments

      narg=iargc()

c First argument is the number of annuli.
c Second is the filename base to which to write the corrections
c Third is the sample angle (only used by flatabs)
c Fourth is the sample height (cylindrical) or width (flatplate).
c Fifth is the azimuthal angle
c After this, the arguments to be read in are (for each annulus):-
c      inner radius
c      outer radius
c      sample density (in atoms/A**3)
c      sample total cross section (in barns/atom)
c      sample capture cross section (in barns/atom)
c
c Hence there should be a total of 5*N+5 arguments, where N is the number of annuli,
c otherwise the procedure will exit with an error message

      iarg=1
      call getarg(iarg,text)
      read(text,*,iostat=ierr) ncont
      if(ierr.eq.0) then

c Check the total number of arguments is correct

         if(narg.ge.(7*ncont+5)) then

            iarg=iarg+1
            call getarg(iarg,filenamebase)
            ilenfilenamebase=index(filenamebase," ")-1

            iarg=iarg+1
            call getarg(iarg,text)
            read(text,*,iostat=ierr) srot(1)

            iarg=iarg+1
            call getarg(iarg,text)
            read(text,*,iostat=ierr) sheight

c Get the azimuthal angle

             iarg=iarg+1
            call getarg(iarg,text)
            read(text,*,iostat=ierr) azidet

                 icont=0
            do while (icont.lt.ncont)

                icont=icont+1
                  iarg=iarg+1
                  call getarg(iarg,text)
                read(text,*) sdimen1(icont)
                  iarg=iarg+1
                  call getarg(iarg,text)
                read(text,*) sdimen2(icont)
                  iarg=iarg+1
                  call getarg(iarg,text)
                read(text,*) srho(icont)
                  iarg=iarg+1
                  call getarg(iarg,text)
                read(text,*) stscat(1,icont)
                  iarg=iarg+1
                  call getarg(iarg,text)
                read(text,*) scaptav(icont)
                  iarg=iarg+1
                  call getarg(iarg,text)
                read(text,*) stscat_out(1,icont)
                  iarg=iarg+1
                  call getarg(iarg,text)
                read(text,*) scaptav_out(icont)

              end do
         else

            write(6,*) 'Not enough parameters in call '
     $,'to calc_corrsx_in_out'

         endif
      endif

c
c calculate the corrections
c
      ind=0
      do j=1,ncont
         nan=0
         equalinout = .false.
         do i=j,ncont
            nan=nan+1
            rad1(nan)=sdimen1(i)
            rad2(nan)=sdimen2(i)
            den(nan)=srho(i)
            captcs(nan)=scaptav(i)
            captcs_out(nan)=scaptav_out(i)
            do il=1,nwavabs
                  sigtl(il,nan)=stscat(il,i)
                  sigtl_out(il,nan)=stscat_out(il,i)

            end do
            equalinout = (abs(sigtl(1,nan)-sigtl_out(1,nan)).lt.0.01)
            write(6,*) den(nan),sigtl(1,nan), sigtl_out(1,nan),equalinout
            
         end do
c
c set up angles for corrections
c
         afirst=0.0
         alast=180.0
         if(theta2theta) then
            afirst=real(ndeg)
            alast=180.0-real(ndeg)
         endif
         call set_corr_ang(sgeom,afirst,alast,azidet,srot(1)
     1,angabs
     1,aziabs,rotabsang,nangabs,naziabs,tthetaabs,nanglesabs)
         call set_corr_ang(sgeom,afirst,alast,azidet,srot(1)
     1,angmul
     1,azimul,rotmulang,nangmul,nazimul,tthetamul,nanglesmul)
c         do ib=1,nanglesabs
c            write(6,*) ib,tthetaabs(ib),rotabsang(ib),angabs(ib)
c     *,aziabs(ib)
c            write(6,*) ib,tthetamul(ib),rotmulang(ib),angmul(ib)
c     *,azimul(ib)
c         end do
             if(sgeom.eq.1) then
            call cylabstof(nan,abstemp
     1,rad1,rad2,sheight,den,sigtl,captcs
     1,sigtl_out,captcs_out)
            if(equalinout) 
     1call cylmultof(nan,rad1,rad2,sheight,den,sigtl,captcs)
         else if(sgeom.eq.2) then
c set the flag for when theta2theta scanning is being used: in this case the final sample absorption correction is
c performed as part of the data normalisation. This is a fix for the fact that the current attenuation correction
c (flat plates) ignores the beam dimensions and its potential interaction with the edges of the sample.
            ignoresample=j.eq.1.and.theta2theta
            call fltabstof(nan,abstemp
     1,rad1,rad2,sheight,den,sigtl,captcs
     1,sigtl_out,captcs_out,ignoresample)
            if(equalinout)
     1call fltmultof(nan,rad1,rad2,sheight,den
     1,sigtl,captcs)
         endif

c save the attenuation correction for this set of sample and containers
c
c the corrections are always stored in order of increasing radius.
c
         do i=1,nan
            ind=ind+1
            do il=1,nwavabs
               do iang=1,nangabs
                  do k=1,naziabs
                     abscor(il,iang,k,ind)=abstemp(il,iang,k,i)
                  end do
               end do
            end do
        write(6,*) ind,abscor(1,1,1,ind)
         end do

c Save the multiple scattering corrections for this annulus

         if(equalinout) then
            do il=1,nwavmul
               do iang=1,nangmul
                  do k=1,nazimul
                     onecor(il,iang,k,j)=onescat(il,iang,k)
                     mulcor(il,iang,k,j)=mulscat(il,iang,k)
                  end do
               end do
            end do
         endif
      end do
      ncorr=ind
      write(6,*) ncorr
     $,' attenuation corrections created in calc_corrsx'
      if(equalinout) write(6,*) ncont
     1,' multiple scattering corrections created in calc_corrsx'
c
c set the output format
c
      nblock=naziabs*ncorr+1
      write(formatout,332) '(',nblock,'(1x,e12.5))'
332      format(a1,i3.3,a11)
      write(6,*) formatout
      lrec=nblock*13
      if(lrec.lt.80) lrec=80
c
c set up filename of file to write
c
      write(fname,103) filenamebase(1:ilenfilenamebase)
103      format(a,'.abs')
c
c Open the attenuation correction file and write it. 
c
      open(10,file=fname,status='unknown',form='formatted'
     *,recl=lrec)
      write(10,*) nwavabs,nangabs,naziabs
      write(10,formatout) (tthetaabs(k),k=1,naziabs)
      do ib=1,nangabs
         write(10,*) ib,tthetaabs(ib)
         do j=1,nwavabs
            write(10,formatout) wavabs(j)
     1,((abscor(j,ib,k,i),i=1,ncorr),k=1,naziabs)
         end do
      end do
      close(10)
c
c Repeat this for the multiple scattering corrections
c
      if(equalinout) then

c set the output format
c
         nblock=nazimul*ncont*2+1
         write(formatout,332) '(',nblock,'(1x,e12.5))'
         write(6,*) formatout
         lrec=nblock*13
         if(lrec.lt.80) lrec=80
c
c set up filename of file to write
c
         write(fname,104) filenamebase(1:ilenfilenamebase)
104      format(a,'.mul')
c
c Open the attenuation correction file and write it. 
c
         open(10,file=fname,status='unknown',form='formatted'
     *,recl=lrec)
         write(10,*) nwavmul,nangmul,nazimul
         write(10,formatout) (tthetamul(k),k=1,nazimul)
         do ib=1,nangmul
            write(10,*) ib,tthetamul(ib)
            do j=1,nwavabs
               write(10,formatout) wavmul(j)
     1,((onecor(j,ib,k,i),mulcor(j,ib,k,i),i=1,ncont),k=1,nazimul)
            end do
         end do
         close(10)
      endif

      return
      end

      subroutine CYLABSTOF(nant,abscort,radt1,radt2,theight,dent
     1,sigtl,abscs,sigtl_out,abscs_out)
c
c      ORIGINAL from AKS
c      modified 27/03/2001 to become a subroutine within GUDRUN
c
c      modified 5/02/2003 to include corrections for azimuthal detectors
c
c      Note that the arrays angabs and aziabs are in sample cylindrical coordinates
c       (as opposed to detector coordinates) so that angabs is relative a z-axis in the 
c      vertically upwards direction
c
      include 'dimension.inc'
      include 'beam.inc'
      include 'abs_corr_azi.inc'
      
      DIMENSION RHO(mcont),SIGT(mcont),abscs(mcont)
      real*4 MUS(mcorrwav,mcont),MUT(mcorrwav,mcont)
     1,SIGTL(mcorrwav,mcont),SIGSL(mcorrwav,mcont)
     1,sigtl_out(mcorrwav,mcont),abscs_out(mcont)
     1,MUT_OUT(mcorrwav,mcont)
      real*4 radt1(mcont),radt2(mcont),dent(mcont)
      real*4 theta,azi
      real*4 abscort(mcorrwav,mcorrang,mcorrang,mcont)
      COMMON /cylabs/ NAN,RAD1(mcont),rad2(mcont),THETA,azi,PI
     1,AMUS(mcont),AMUT(mcont),DEN(mcont),abstemp(mcont)
     1,AMUT_OUT(mcont)
      PI=4.0*atan(1.0)
      astep=stepa
      height=theight
      nan=nant
      do i=1,nan
            den(i)=dent(i)
            rad1(i)=radt1(i)
            rad2(i)=radt2(i)
C
C CALCULATE  SCATTERING C/S AT EACH WAVELENGTH ASSUMING 1/LAMBDA ABSORPTION
C
             DO IL=1,nwavabs
                  SIGSL(IL,I)=SIGTL(IL,I)-ABSCS(I)
                  MUS(IL,I)=DEN(I)*SIGSL(IL,I)
                  MUT(IL,I)=DEN(I)*SIGTL(IL,I)
                  MUT_OUT(IL,I)=DEN(I)*SIGTL_OUT(IL,I)
            end do
!            write(6,*) i,den(i),rad1(i),rad2(i),mut(1,i),mut_out(1,i)
      end do
      MS=(RAD2(1)-RAD1(1)+0.0001)/ASTEP
      IF (MS.LT.1)MS=1
      do iazi=1,naziabs
         azi=aziabs(iazi)*pi/180.0
         DO I=1,nangabs
            THETA=angabs(I)*pi/180.0
            DO J=1,nwavabs
               DO IR=1,NAN
                  AMUS(IR)=MUS(J,IR)
                  AMUT(IR)=MUT(J,IR)
                  AMUT_OUT(IR)=MUT_OUT(J,IR)
               end do
               nanw=nan-1
               CALL ACYL(MS,AREAS,AREAC)
               do k=1,nan
                  abscort(j,i,iazi,k)=abstemp(k)
!                  write(6,*) ms,aziabs(iazi),k,abstemp(k)
               end do
            end do
         end do
      end do
      return
      END

      subroutine fltabstof(ncan,abscort,rad1,rad2,width,rho,sigtl,abscs
     1,sigtl_out,abscs_out,ignoresample)
c
c to calculate flat plate absorption factors
c
c      original   by aks
c      modified   by wsh    10-3-88   for use with coral
c      modified by aks 16-7-90 to put out a unity absorption factor
c if secondary angle is close to 90 degrees.
c
c      completely revamped by aks 17-05-01 to become a subroutine of the gudrun
c suite
c
c      modified 5/02/2003 to include situation with azimuthal detectors
c      modified January 2012 to replace inner and outer dimensions with upstream and downstream slabs
c      also corrected logic which was flawed in some cases. 
c
      implicit none
      include 'dimension.inc'
      include 'beam.inc'
      include 'abs_corr_azi.inc'
      integer ncan            !no. of containers + sample
      integer ic,ic1,ib,il
      real rot              !angle of flat plate normal to beam
      real rad1(mcont)      !dimension of upstream slab
      real rad2(mcont)      !dimension of downstream slab
      real width,width2            !lateral half width of sample and containers
      real rho(mcont)            !atomic number density of sample and containers
      real sigtl(mcorrwav,mcont)      !total cross-sections
      real abscs(mcont)            !capture cross section (not used)
      real mut(mcont),amut      !temporary atten. coeff. 
      real theta,theta1      !temporary angle values
      real fs_up,fs_dn,exp_dn,exp_up      !temporary exponential values
      real pi,piconv      
      real tsec,sec1,sec2      !secants of angles
      real ts(mcont)            !wall thickness of each container
      real ts_in_up,ts_out_up,ts_in_dn,ts_out_dn      !in and out flight path through each container
      real pl_in_dn(mcont)      !total flight path into downstream slab
      real pl_out_dn(mcont)      !total flight path out of downstream slab
      real pl_in_up(mcont)      !total flight path into uptream slab
      real pl_out_up(mcont)      !total flight path out of upstream slab
      real f            !function to calculate integral of exponential
      real abscort(mcorrwav,mcorrang,mcorrang,mcont)
      real sigtl_out(mcorrwav,mcont)
      real abscs_out(mcont)
      real mut_out(mcont),amut_out      !temporary atten. coeff. 
      real footprintfac,diskfootprint
      logical ignoresample  !For flat plate samples, tells cylabstof not to include the sample when calculating As,sc.

      pi=3.141592653
      piconv=pi/180.
      width2=0.5*width
c      rot=rotabsang(1)*piconv
c      sec1=abs(1./cos(rot))
c
c step through banks and calculate corrections at each wavelength and
c scattering angle
c
      do ib=1,nangabs
         theta1=angabs(ib)
         rot=rotabsang(ib)*piconv
         sec1=abs(1./cos(rot))
c         footprintfac=0.5*(diskfootprint(width2,hbup,a,sec1)
c     1+diskfootprint(width2,hbup,b,sec1))
c         write(6,*) ib,width2,hbup,a,b,sec1,footprintfac
         tsec=theta1
c
c tsec is the angle the scattered beam makes with the normal to the sample
c surface.  if abs(tsec) is close to 90 deg. calculation of absorption
c coefficients is unreliable
c
         if(abs(abs(tsec)-90.).gt.0.1) then
            tsec=tsec*piconv
            sec2=1./cos(tsec)
c
c first form the sum of path lengths (pl) through the other slabs into and out
c of each slab. the sample and containers are divided into two slabs, dn and up.
c up refers to the slab on the upstream side of the centre of the sample (i.e.
c facing towards the neutron source). dn refers to slabs on the downstream side 
c of the sample centre. in refers to the flight path into the slab. out refers
c to the flight out of the slab towards the detector.
c
            do il=1,nwavabs
               do ic=1,ncan
                  pl_in_up(ic)=0.0
                  pl_out_up(ic)=0.0
                  pl_in_dn(ic)=0.0
                  pl_out_dn(ic)=0.0
               end do
               do ic=1,ncan
                  mut(ic)=rho(ic)*sigtl(il,ic)
                  mut_out(ic)=rho(ic)*sigtl_out(il,ic)
                  ts(ic)=(rad2(ic)+rad1(ic))
                  ts_in_up=mut(ic)*rad1(ic)*sec1
                  ts_out_up=mut_out(ic)*rad1(ic)*abs(sec2)
                  ts_in_dn=mut(ic)*rad2(ic)*sec1
                  ts_out_dn=mut_out(ic)*rad2(ic)*abs(sec2)

                  do ic1=1,ncan
                     if(sec2.gt.0.0) then
                        if(ic1.lt.ic) then
                           pl_in_up(ic1)=pl_in_up(ic1)+ts_in_up
                           pl_out_up(ic1)=pl_out_up(ic1)+ts_out_dn
                           pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                           pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_dn
                        else if(ic1.eq.ic.and.(ic.eq.1.and.
     *.not.ignoresample)) then
                           pl_out_up(ic1)=pl_out_up(ic1)+ts_out_dn
                           pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                        else if(ic1.gt.ic) then
                           pl_out_up(ic1)=pl_out_up(ic1)+ts_out_up
     *+ts_out_dn
                           pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
     *+ts_in_dn
                        endif
                     else if(sec2.lt.0.0) then
                        if(ic1.lt.ic) then
                           pl_in_up(ic1)=pl_in_up(ic1)+ts_in_up
                           pl_out_up(ic1)=pl_out_up(ic1)+ts_out_up
                           pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                           pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_up
                        else if(ic1.eq.ic.and.(ic.ne.1.or.
     *.not.ignoresample)) then
                           pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
                           pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_up
                        else if(ic1.gt.ic) then
                           pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in_up
     *+ts_in_dn
                           pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out_up
     *+ts_out_dn
                        endif
                     endif
                  end do
               end do
      if(ib.eq.1.and.il.eq.1) write(6,*) 
     *(pl_in_up(ic),pl_out_up(ic),pl_in_dn(ic),pl_out_dn(ic),ic=1,ncan)
c
c calculate absorption integrals for each slab and form attenuation correction.
c
               do ic=1,ncan
                  amut=mut(ic)
                  amut_out=mut_out(ic)
                  if(ic.ne.1.or..not.ignoresample) then
                     fs_up=f(amut,amut_out,rad1(ic),sec1,sec2)
                     fs_dn=f(amut,amut_out,rad2(ic),sec1,sec2)
                  else
                     fs_up=1.0
                     fs_dn=1.0
                  endif
                  exp_up=exp(-(pl_in_up(ic)+pl_out_up(ic)))
                  exp_dn=exp(-(pl_in_dn(ic)+pl_out_dn(ic)))
                  abscort(il,ib,1,ic)=(rad1(ic)*exp_up*fs_up
     *+rad2(ic)*exp_dn*fs_dn)/(ts(ic))
c If no beam compensation is occurring, then multiply these values by sec1 to represent the beam footprint at this angle
c                  if(.not.beamcompensation) then
c                     abscort(il,ib,1,ic)=
c     *footprintfac*abscort(il,ib,1,ic) 
c                  endif
c      write(6,*) il,ib,ic,amut,amut_out,fs,exp_dn,exp_up,ts(ic)
               end do
            end do
c
c case where tsec is close to 90 - no attenuation correction can be calculated
c
         else

            do ic=1,ncan
               do il=1,nwavabs
                  abscort(il,ib,1,ic)=1.0
               end do
            end do
         endif
c      write(6,*) ib,angabs(ib),rotabsang(ib),sec1,sec2
      end do
      return
      end

      
