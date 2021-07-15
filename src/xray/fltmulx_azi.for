      subroutine fltmultof(ncan,rad1,rad2,width,rho,sigtl,captcs)
c
c      original   by aks
c      modified   by wsh   1-3-89    for coral input
c 14-5-90  : modified by aks and ach to increase size of array thetab for
c            detector angles from 8 to 200 since lad requires 14 for a flat
c            plate with standard grouping.
c 16-7-90  : for secondary angles near 90 degrees the program puts
c out unity for the primary scattering and zero for the multiple scattering
c
c completely revamped by aks 17-05-01 to become a subroutine of the gudrun
c suite
c
c      modified 5/02/2003 to include situation with azimuthal detectors
c
      include 'dimension.inc'
      include 'beam.inc'
      include 'mul_corr_azi.inc'
      integer ncan            !no. of containers + sample
      integer ic,ic1,ib,il
      real rot
      real rad1(mcont)      !dimension of upstream slab
      real rad2(mcont)      !dimension of downstream slab
      real width,width2      !lateral width and half width of sample
      real rho(mcont)            !atomic number density of sample and containers
      real sigtl(mcorrwav,mcont)      !total cross-sections
      real captcs(mcont)      !capture cross section
      real amut(mcorrwav)      !average total attenuation coeffsicients
      real amus(mcorrwav)      !average scattering attenuation coefficients
      real amuabs            !average capture cross section
      real tsum                  !total thickness of slab
      real theta,theta1      !temporary angle values
      real pi,piconv      
      real tsec,sec1,sec2      !secants of angles
      real ts(mcont)            !wall thickness of each container
      real leakfac            !leakage factor to correct for finite width
c
c internal arrays used by the integration routines
c
      real*4 ainten(mslice),sumint(mslice),sumnx(mslice),eout(mslice)

      open(20,file='fltmuldiag.txt',status='unknown')

      pi=3.141592653
      piconv=pi/180.
      pi4=4.*pi
      pi2=2.*pi
      gamma=0.5772156649
      width2=0.5*width

c
c specified accuracy of integrals
c
      accur=0.0001
C
C form average total and capture cross sections for sample and containers
C
      tsum=0.0
      do ic=1,ncan
            ts(ic)=rad2(ic)+rad1(ic)
            tsum=tsum+ts(ic)
      end do
c
c total thickness of slab
c
c      t=2.0*tsum
      t=tsum

      write(6,*)
      write(6,*) tsum

c
c define thickness of step for integration
c
      x=t/real(nslice)
      do il=1,nwavmul
         amut(il)=0.0
         amuabs=0.0
         do ic=1,ncan
            amut(il)=amut(il)+ts(ic)*rho(ic)*sigtl(il,ic)/tsum
            amuabs=amuabs+ts(ic)*rho(ic)*captcs(ic)/tsum
         end do
         amus(il)=amut(il)-amuabs
      end do
c
c step through banks and calculate corrections at each wavelength and
c scattering angle
c
c      write(20,'(a)') 'if(ib.eq.1) write(20,*) i,ainten(i),eout(i)'
c     *//',sec1,al1,al2,a,b,amuxs,sumnx(i)'
      do ib=1,nangmul
         theta1=angmul(ib)
         rot=rotmulang(ib)*piconv
         sec1=1./cos(rot)
         tsec=theta1
c
c tsec is the angle the scattered beam makes with the normal to the sample
c surface.  if abs(tsec) is close to 90 deg. calculation of multiple scattering
c coefficients is unreliable
c
         if(abs(abs(tsec)-90.).gt.0.1) then
            tsec=tsec*piconv
            sec2=1./cos(tsec)
            if(sec1.lt.1..and.sec1.gt.0.) sec1=1.
            if(sec1.gt.-1..and.sec1.lt.0.) sec1=-1.
            tes=t
c            if(sec2.lt.0.) tes=0.
c
c step through wavelengths
c
            do il=1,nwavmul
c
c calculate abs c/s
c
               amu=amut(il)
               amusct=amus(il)
               unit=amusct/pi4
               fac=amusct*t/pi4
c
c calculate leakage factor for neutrons out the side of the slab
c
               leakfac=width2*amu
               if(leakfac.lt.20.0) then
                  leakfac=1.0-exp(-leakfac)
               else
                  leakfac=1.0
               endif
c
c calculate scattering from a point to each plane
c
               xx=0.
               ainten(1)=pi2*unit*binte1(amu,xx,x)
               do i=2,nslice
                  x1=(i-1)*x
                  ainten(i)=unit*pi2*binte1(amu,x1,x)
c      write(20,*) i,amu,x1,x,ainten(i)
              end do
c
c calculate primary intensity at each plane
c
               xh=x*0.5
               do i=1,nslice
                  x1=(i-1)*x+xh
                  al0=x1
                  if(sec2.gt.0.0) then
                     al1=tes-x1
                  else
                     al1=x1
                  endif
                  eout(i)=aintex(al1,abs(sec2),amu,x)/pi4
                  al1=amu*x1
c                  al2=abs(amu*tes-al1)
                  al2=amu*tes-al1
                  amuxs=al1*sec1
                  aa=e1int1(sec1,al1)
                  bb=e1int2(sec1,al2)
                  sumnx(i)=(aa+bb*exp(-amuxs))*0.5*amusct**2/amu
c      if(ib.eq.1) write(20,*) i,ainten(i),eout(i),sec1,al1,al2
c     *,aa,bb,amuxs,sumnx(i)
               end do
c
c calculate nord orders of scattering
c
c      write(20,'(a)') 'if(ib.eq.1) write(20,*) ib,sec1,i,sum,sumtot'
               sumtot=0.
               i=1
               sum=0.
               do while (sum.ge.accur*sumtot)
                  sum=0.0
                  if(i.eq.1) then
                     sum=prime(sec1,sec2,amu,t)*fac
                  else
                     do j=1,nslice
                        sum=sum+sumint(j)*eout(j)
                     end do
                  endif
                  i1=i
                  sumtot=sumtot+sum
                  if(i.eq.1) first=sumtot
                  if(sum.ge.accur*sumtot) then
                     if(i.gt.1) then
c
c calculate next order
c
                        do j=1,nslice
                           sumnx(j)=0.
                           do k=1,nslice
                              iref=iabs(j-k)+1
                              sumnx(j)=sumnx(j)+ainten(iref)*sumint(k)
                           end do
                        end do
                     endif
                     do j=1,nslice
                        sumint(j)=sumnx(j)*x*leakfac
                     end do
                  endif
                  if(ib.eq.1) write(20,*) ib,sec1,i,sum,sumtot
                  i=i+1
               end do
               rest=sumtot-first
               onescat(il,ib,1)=first
               mulscat(il,ib,1)=rest
            end do
         else
            do il=1,nwavmul
               first=1.0
               rest=0.0
               onescat(il,ib,1)=first
               mulscat(il,ib,1)=rest
            end do
         endif
c If no beam compensation is occurring, then multiply these values by sec1 to represent the beam footprint at this angle
c         if(.not.beamcompensation) then
c            fac=0.5*(diskfootprint(width2,hbup,a,sec1)
c     1+diskfootprint(width2,hbup,b,sec1))
c            write(20,*) ib,width2,hbup,a,b,sec1,fac
c            do il=1,nwavmul
c               onescat(il,ib,1)=fac*onescat(il,ib,1)
c               mulscat(il,ib,1)=fac*mulscat(il,ib,1)
c            end do
c         endif
      end do
      close(20)
      return
      end
