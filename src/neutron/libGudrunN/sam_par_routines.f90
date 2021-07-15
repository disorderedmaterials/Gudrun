!     
! File:   get_sample_par.f90
! Author: aks45
!
! Created on 05 November 2013, 08:48
!
MODULE sam_par_routines
    
    
    implicit none
    
    CONTAINS
!***********************************************************************************
!*
!*   get_sam_par.FOR
!*
!*   A K Soper, March 2001
!*   
!*   Reads the vanadium and sample parameters
!*
!***********************************************************************************
    subroutine get_sam_par(nin,nread,run,perno,wavemin,wavemax,wavestep)

      use inputfilestrings
!      use reallocation_routines
      use run_par
      use groups_routines
      use beam_routines
      use van_par
      use sam_par
      use write_routines
      use math_routines
      use interpolation_routines
      use cylmultof_routines
        
        integer nin      !unit number to read data from 
   integer nread      !0 for vanadium, 1 for sample, 2 for containers
   integer ilenfname      !no. of lines to skip at start of c/s tables
   character(len=256) run      !run number of file for mut, mul, abs
   integer perno      !period number to create
   integer il,iv,is,ierr      !do loop counter
   integer nlambt   !temporary counter of no. of wavelengths
   integer massno   !mass number of element (0 for natural element)
   integer atomicno,atomicnosave   !atomic number of element
   integer nsum      !check integer on number of elements found
        integer sgeomtrial      !Temporary store of sample or container geometry - not used in evaluating corrections - just in case sample or container geometry is incorrectly specified.
   real, dimension(:), allocatable     :: templa,tempcs,tempscs!(mcorrwav)   !temporary array of wavelengths
   real xx      !temporary value
   real area,vol,sthick   !dummy variables
   real pi,fsum,asum,bsum,b2sum,ssum,csum,wave
   real atwt      !atomic weight of element or isotope
   real bscat      !bound scattering length of element or isotope
   real scatcs      !scattering cross section and average
   real captcs      !capture cross section and average
   real dummy      !dummy read variable
   real wavemin,wavemax   !minimum and maximum wavelengths to be used
   real wavestep,wavestep2      !step in wavelength for corrections
   real wave1,wave2
   real amass,afrac
   real AR,ARB,ARS
   character*256 fname   !filename for c/s data
   character*2 symbol,symbolsave   !elemental symbol
   character*6 spin   !nuclear spin and parity
   character*1 hash

        pi=4.0*atan(1.0)
   wavestep2=0.5*wavestep
        if(nread.eq.0) then
!c
!c read composition of vanadium
!c
            iv=0
            mvelement=0
            fsum=0
            symbol='ZZ'
            afrac=1.0
            do while (symbol.ne.'*')
                read(nin,*) symbol,amass,afrac
                if(symbol.ne.'*') then
                    iv=iv+1
                    do while (iv.gt.mvelement)
                        mvelement=mvelement+5
              call reallocate1d_c(vsymbol,len(vsymbol),mvelement)
              call reallocate1d_i(vmass,mvelement)
              call reallocate1d_r(vfrac,mvelement)
                    end do
                    vsymbol(iv)=symbol
                    vmass(iv)=amass
                    vfrac(iv)=afrac
                    fsum=fsum+vfrac(iv)
                endif
            end do
            nvelement=iv
            if(nvelement.gt.0.and.fsum.gt.0.0) then
                write(6,*)
      write(6,*) 'No of elements found for vanadium: ',nvelement
      do iv=1,nvelement
                    write(6,*) vsymbol(iv),vmass(iv)
                end do
      write(6,*)
            else
      write(6,*)
      write(6,*) 'Must be at least 1 element for vanadium'
      write(6,*)
      stop
            endif
!
!c calculate atomic fractions of components
!c
            do iv=1,nvelement
                vfrac(iv)=vfrac(iv)/fsum
            end do
            call reallocate1d_r(vatwt,nvelement)
            call reallocate1d_r(vscatlen,nvelement)
            call reallocate1d_r(vscatcs,nvelement)
            call reallocate1d_r(vcaptcs,nvelement)

!c
!c open the cross section and mass tables and find the scattering cross section
!c and capture cross section for each element in vanadium. 
!c
            open(99,file=neutronparfilename,status='old',iostat=ierr)
            if(ierr.ne.0) then
                neutronparfilename=convertfilename(neutronparfilename)
                open(99,file=neutronparfilename,status='old',iostat=ierr)
            endif
            if(ierr.eq.0) then
!c
!c counter for the number of element matches found in the cross-section table
!c
      nsum=0
      asum=0.0
      bsum=0.0
      ssum=0.0
      csum=0.0
      symbolsave='  '
      atomicnosave=0
!c
!c format to read the c/s tables
!c
99   format(a2,1x,i2,1x,i3,1x,a6,1x,f9.5,4(1x,f8.4),1x,f9.4,1x,e9.3)
!c
!c loop through the cross section file looking for matches for the elements
!c in the list
!c
                do while (ierr.ge.0.and.nsum.lt.nvelement)
                    read(99,99,iostat=ierr) symbol,atomicno,massno,spin,atwt,bscat,dummy,dummy,dummy,scatcs,captcs
                    if(ierr.eq.0) then
                        if(symbol.ne.symbolsave.and.symbol.ne.'  ') symbolsave=symbol
                        if(atomicno.ne.atomicnosave.and.atomicno.ne.0) atomicnosave=atomicno
!c
!c step through the elements and see if any match this one.
!c
                        bscat=0.1*bscat
                        do iv=1,nvelement
                            if(symbolsave.eq.vsymbol(iv).and.massno.eq.vmass(iv)) then
                                nsum=nsum+1
                                vatwt(iv)=atwt
                                vscatlen(iv)=bscat
                                vscatcs(iv)=scatcs
                                vcaptcs(iv)=captcs
                                asum=asum+atwt*vfrac(iv)
                                bsum=bsum+bscat*vfrac(iv)
                                ssum=ssum+scatcs*vfrac(iv)
                                csum=csum+captcs*vfrac(iv)
                            endif
                        end do
                    end if
                end do
      close(99)
      if(nsum.eq.nvelement) then
                    vatwtav=asum
                    vscatlenav=bsum
                    vscatav=ssum
                    vcaptav=csum
      else
                    write(6,*) 'No. of elements input ',nvelement &
        ,' NOT EQUAL to no. of elements found in c/s file', nsum
                    write(6,*) 'Check input file'
                    stop
      endif
            else
                write(6,*) 'Neutron cross section file not found'
                stop
            end if
!c Next line will specify the geometry of this sample or container

            call getaline(nin)
            write(6,*) 'get_van_par> ',nwords,ncf(1),ncl(1)
            if(line(ncf(1):ncl(nwords)).eq.'CYLINDRICAL') then
               vgeom=1

            else if(line(ncf(1):ncl(nwords)).eq.'FLATPLATE') then

               vgeom=2

            endif
!c
!c vanadium inner and outer dimension
!c
            read(nin,*) vdimen1,vdimen2
!c
!c for cylinders, read vanadium overall height
!c
            call getaline(nin)
            if(vgeom.eq.1) then
                if(nwords.eq.2) then
                    read(line(ncf(2):ncl(2)),*) vheight
                else
                    read(line(ncf(1):ncl(1)),*) vheight
                end if
                vrot=0.0
                vwidth=vheight
            else
!c
!c otherwise for flat plates only, read rotation angle and full width of vanadium
!c
                if(nwords.eq.2) then
                    read(line(ncf(1):ncl(2)),*) vrot,vwidth
                else
                    read(line(ncf(1):ncl(1)),*) vwidth
                    vrot=0.0
                end if
                vheight=vwidth
            endif
!c
!c vanadium density in gm/cm**3. If -ve, then it is assumed to be in atoms/A**3

!c
            read(nin,*) vrho
            if(vrho.lt.0.0) then
                vrho=abs(vrho)
            else
                vrho=vrho*0.60221/vatwtav
            endif
!c
!c vanadium temperature (K) for Placzek correction
!c --- if zero, no correction will be performed.
!c
            read(nin,*) vtemp
!c 
!c vanadium c/s file
!c
            call getaline(nin)
            if(nwords.gt.0) then
                ilenfname=ncl(1)-ncf(1)+1
                fname=line(ncf(1):ncl(1))
            endif
!c
!c if the word 'TABLES' is read, then programme generates the total cross-section
!c from the table of cross sections in steps of wavestep.
!c
            if(fname.eq.'TABLES') then
                nlambv=nint((wavemax-wavemin)/wavestep)
                call reallocate1d_r(vlamb,nlambv)
                call reallocate1d_r(vtscat,nlambv)
                call reallocate1d_r(vsscat,nlambv)
                do il=1,nlambv
                    vlamb(il)=wavemin+real(il)*wavestep-wavestep2 !Need point format if using the tables
                    vtscat(il)=vscatav+vlamb(il)*vcaptav/1.798
                    vsscat(il)=vscatav
                end do
!c
!c write the data to an appropriate .mut file
!c
      fname=run
                call change_ext(fname,newext='mut')
                write(fname,'(a,i2.2)') fname(1:len_trim(fname)),perno
                open(10,file=fname,status='unknown')
                write(10,*) nlambv
                do il=1,nlambv
                    write(10,101) vlamb(il),vtscat(il),vsscat(il)
                end do
101             format(3(1x,e13.6))
      close(10)
      vtmon=0
!c
!c if the word 'TRANSMISSION' is read then that is used to set a flag
!c to calculate the total cross-section from the transmission monitor when that
!c becomes available
!c
            else if(fname.eq.'TRANSMISSION') then
                vtmon=1
                nlambv=nint((wavemax-wavemin)/wavestep)+1
                call reallocate1d_r(vlamb,nlambv)
                call reallocate1d_r(vtscat,nlambv)
                call reallocate1d_r(vsscat,nlambv)
                do il=1,nlambv
                    vlamb(il)=wavemin+real(il-1)*wavestep !Need histogram for tranmission monitor
                end do
            else
                open(16,file=fname,status='old')
                read(16,*) hash,nlambv
                call reallocate1d_r(vlamb,nlambv)
                call reallocate1d_r(vtscat,nlambv)
                call reallocate1d_r(vsscat,nlambv)
                do il=1,nlambv
                    read(16,'(a)') line
                    call parse()
                    if(nwords.gt.2) then
                        read(line,*) vlamb(il),vtscat(il),vsscat(il)
                    else if(nwords.eq.2) then
                        read(line,*) vlamb(il),vtscat(il)
                        vsscat(il)=vtscat(il) !Set the scattering cross section to the total cross-section if not specified.
                    end if
                end do
                close(16)
                vtmon=0
            endif
   else 
            if(nread.eq.1) then
      ncont=1
                mselement=0
            else
      ncont=ncont+1
            endif
            call reallocate1d_i(ncvals,ncont)
            call reallocate2d_c(ssymbol,len(ssymbol),mselement,ncont)
            call reallocate2d_i(smass,mselement,ncont)
            call reallocate2d_r(sfrac,mselement,ncont)
            ncvals(ncont)=0

!c
!c read composition of sample
!c
            is=0
            fsum=0
            symbol='ZZ'
            afrac=1.0
            do while (symbol.ne.'*')
               read(nin,*) symbol,amass,afrac
               if(symbol.ne.'*') then
                    is=is+1
                    do while(is.gt.mselement)
                        mselement=mselement+5
              call reallocate2d_c(ssymbol,len(ssymbol),mselement,ncont)
              call reallocate2d_i(smass,mselement,ncont)
              call reallocate2d_r(sfrac,mselement,ncont)
                    end do
                    ssymbol(is,ncont)=symbol
                    smass(is,ncont)=amass
                    sfrac(is,ncont)=afrac
                    fsum=fsum+sfrac(is,ncont)
                endif
            end do
            call reallocate1d_i(nselement,ncont)
            nselement(ncont)=is
            if(nselement(ncont).gt.0.and.fsum.gt.0.0) then
      write(6,*)
      write(6,*) 'No of elements found for sample/container ',ncont,' : ',nselement(ncont)
      write(6,*)
            else
                write(6,*)
      write(6,*) 'Must be at least 1 element defined'
      write(6,*)
      stop
            endif
!c
!c calculate atomic fractions of components
!c
            do is=1,nselement(ncont)
                sfrac(is,ncont)=sfrac(is,ncont)/fsum
            end do
            call reallocate2d_r(satwt,mselement,ncont)
            call reallocate2d_r(sscatlen,mselement,ncont)
            call reallocate2d_r(sscatcs,mselement,ncont)
            call reallocate2d_r(scaptcs,mselement,ncont)
!c
!c open the cross section and mass tables and find the scattering cross section
!c and capture cross section for each element in vanadium. 
!c
            open(99,file=neutronparfilename,status='old',iostat=ierr)

            if(ierr.ne.0) then

                neutronparfilename=convertfilename(neutronparfilename)
                open(99,file=neutronparfilename,status='old',iostat=ierr)

            endif

!c
!c counter for the number of element matches found in the cross-section table
!c
            nsum=0
            asum=0.0
            bsum=0.0
            b2sum=0.0
            ssum=0.0
            csum=0.0
            symbolsave='  '
            atomicnosave=0.0
!c
!c loop through the cross section file looking for matches for the elements
!c in the list
!c
            do while (ierr.ge.0.and.nsum.lt.nselement(ncont))
                read(99,99,iostat=ierr) symbol,atomicno,massno,spin,atwt,bscat,dummy,dummy,dummy,scatcs,captcs
                if(ierr.eq.0) then
                    if(symbol.ne.symbolsave.and.symbol.ne.'  ') symbolsave=symbol
                    if(atomicno.ne.atomicnosave.and.atomicno.ne.0) atomicnosave=atomicno
!c
!c step through the elements and see if any match this one.
!c
!c
!c convert bscat to 10**-12cm
!c
                    bscat=0.1*bscat
                    do is=1,nselement(ncont)
                        if(symbolsave.eq.ssymbol(is,ncont).and.massno.eq.smass(is,ncont)) then
                            nsum=nsum+1
                            satwt(is,ncont)=atwt
                            sscatlen(is,ncont)=bscat
                            sscatcs(is,ncont)=scatcs
                            scaptcs(is,ncont)=captcs
                            asum=asum+atwt*sfrac(is,ncont)
                            bsum=bsum+bscat*sfrac(is,ncont)
                            b2sum=b2sum+(bscat*bscat)*sfrac(is,ncont)
                            ssum=ssum+scatcs*sfrac(is,ncont)
                            csum=csum+captcs*sfrac(is,ncont)
                        endif
                    end do
                end if
            end do
            close(99)
            call reallocate1d_r(satwtav,ncont)
            call reallocate1d_r(sscatlenav,ncont)
            call reallocate1d_r(sscatlensqav,ncont)
            call reallocate1d_r(sscatav,ncont)
            call reallocate1d_r(scaptav,ncont)
            if(nsum.eq.nselement(ncont)) then
                satwtav(ncont)=asum
                sscatlenav(ncont)=bsum
                sscatlensqav(ncont)=b2sum
                sscatav(ncont)=ssum
                scaptav(ncont)=csum
            else
      write(6,*) 'No. of elements input ',nselement(ncont),' NOT EQUAL to no. of elements found in c/s file', nsum
                write(6,*) 'Check input file'
                stop
            end if
            !c Next line will specify the geometry of this sample or container
            call getaline(nin)
            sgeomtrial=sgeom
            if(line(ncf(1):ncl(nwords)).eq.'CYLINDRICAL') then
               sgeomtrial=1
            else if(line(ncf(1):ncl(nwords)).eq.'FLATPLATE') then
               sgeomtrial=2
            endif
!c
!c sample inner and outer dimension
!c
            call reallocate1d_r(sdimen1,ncont)
            call reallocate1d_r(sdimen2,ncont)
            call reallocate1d_r(srot,ncont)
            call reallocate1d_r(swidth,ncont)
            call reallocate1d_r(srho,ncont)
            call reallocate1d_r(stemp,ncont)
            call reallocate1d_i(stmon,ncont)
            call reallocate1d_r(tweak,ncont)
            call reallocate1d_r(snorm,ncont)
            read(nin,*) sdimen1(ncont),sdimen2(ncont)
!c
!c for flat plates only, read rotation angle of sample and full width of sample, but set the sample height in case the specified geometry is incorrect.
!c
            call getaline(nin)
            if(sgeomtrial.eq.1) then
                if(nwords.eq.2) then
                    read(line(ncf(2):ncl(2)),*) sheight
                else
                    read(line(ncf(1):ncl(1)),*) sheight
                end if
                srot(ncont)=0.0
                swidth(ncont)=sheight
            else if(sgeomtrial.eq.2) then
!c
!c otherwise get the overall sample height - this is assumed to be the same
!c for all the subsequent containers
!c
                if(nwords.eq.2) then
                    read(line(ncf(1):ncl(2)),*) srot(ncont),swidth(ncont)
                else
                    read(line(ncf(1):ncl(1)),*) swidth(ncont)
                    srot(ncont)=0.0
                end if
                sheight=swidth(ncont)
            endif
!c
!c sample number density. +ve means gm/cm**3, -ve means atoms/A**3
!c
            read(nin,*) srho(ncont)
            if(srho(ncont).lt.0.0) then
                srho(ncont)=abs(srho(ncont))
            else
      srho(ncont)=srho(ncont)*0.60221/satwtav(ncont)
            endif
            write(6,*) ncont,sscatav(ncont),scaptav(ncont),srho(ncont)
!c
!c sample temperature (K) for Placzek correction (only for ncont = 1)
!c --- if zero, no correction will be performed.
!c
            if(ncont.eq.1) read(nin,*) stemp(ncont)
!c 
!c sample c/s file
!c
            call getaline(nin)
            if(nwords.gt.0) then
                ilenfname=ncl(1)-ncf(1)+1
                fname=line(ncf(1):ncl(1))
            endif
!c
!c if the word 'TABLES' is read, then programme generates the total cross-section
!c from the table of cross sections in steps of wavestep
!c
            if(fname.eq.'TABLES'.or.index(fname,'.mint').ne.0) then
                if(ncont.eq.1) then
                    nlambs=nint((wavemax-wavemin)/wavestep) !Only the sample can change the wavelength scale
                    call reallocate1d_r(slamb,nlambs)
                endif
                call reallocate2d_r(stscat,nlambs,ncont)
                call reallocate2d_r(ssscat,nlambs,ncont)
                do il=1,nlambs
                    if(ncont.eq.1) slamb(il)=wavemin+real(il)*wavestep-wavestep2 !Point format if using tables
                    stscat(il,ncont)=sscatav(ncont)+slamb(il)*scaptav(ncont)/1.798
                    ssscat(il,ncont)=sscatav(ncont)
                end do

!c If this is a mint file, read the data and process it

                call reallocate1d_c(smintname,len(smintname),ncont)
                if(index(fname,'.mint').ne.0) then
                    smintname(ncont)=fname
                    call readandprocessmint(fname,nlambs,slamb,tempcs,sscatav(ncont),scaptav(ncont),srho(ncont),ncont)
                    do il=1,nlambs
                        stscat(il,ncont)=tempcs(il)
                    end do
                else
!c Store the default in smintname
                    smintname(ncont)='*'
                    ncvals(ncont)=0
                endif
!c
!c write the data to an appropriate .mut file
!c
                fname=run
                call change_ext(fname,newext='mut')
                write(fname,'(a,i2.2)') fname(1:len_trim(fname)),perno
                open(10,file=fname,status='unknown')
                write(10,*) nlambs
                do il=1,nlambs
                    write(10,101) slamb(il),stscat(il,ncont),ssscat(il,ncont)
                end do
                close(10)
                stmon(ncont)=0
!c
!c if the word 'TRANSMISSION' is read then that is used to set a flag
!c to calculate the total cross-section from the transmission monitor when that
!c becomes available
!c
            else if(fname.eq.'TRANSMISSION') then
!c
!c set the flag for using the transmission monitor data
!c
                stmon(ncont)=1
                if(ncont.eq.1) then
                    nlambs=nint((wavemax-wavemin)/wavestep)+1 !Only the sample can set the wavelength scale
                    call reallocate1d_r(slamb,nlambs)
                endif
                call reallocate2d_r(stscat,nlambs,ncont)
                call reallocate2d_r(ssscat,nlambs,ncont)
                if(ncont.eq.1) then
                    do il=1,nlambs
                        slamb(il)=wavemin+real(il-1)*wavestep !Need histogram for tranmission monitor
                    end do
                end if
            else
                write(6,*) fname(1:ilenfname)
                open(16,file=fname,status='old')
                if(ncont.eq.1) then
                    read(16,*) hash,nlambs
                    write(6,*) hash,nlambs
                    call reallocate1d_r(slamb,nlambs)
                    call reallocate2d_r(stscat,nlambs,ncont)
                    call reallocate2d_r(ssscat,nlambs,ncont)
                    do il=1,nlambs
                        read(16,'(a)') line
                        call parse()
                        if(nwords.gt.2) then
                            read(line,*) slamb(il),stscat(il,ncont),ssscat(il,ncont)
                        else if(nwords.gt.1) then
                            read(line,*) slamb(il),stscat(il,ncont)
                            ssscat(il,ncont)=stscat(il,ncont)
                        end if
                    end do
                    close(16)
                else
                    read(16,*) hash,nlambt
                    call reallocate2d_r(stscat,nlambs,ncont)
                    call reallocate2d_r(ssscat,nlambs,ncont)
                    call reallocate1d_r(templa,nlambt)
                    call reallocate1d_r(tempcs,nlambt)
                    call reallocate1d_r(tempscs,nlambt)
                    do il=1,nlambt
                        read(16,'(a)') line
                        call parse()
                        if(nwords.gt.2) then
                            read(line,*) templa(il),tempcs(il),tempscs(il)
                        else if(nwords.gt.1) then
                            read(line,*) templa(il),tempcs(il)
                            tempscs(il)=tempcs(il)
                        end if
                    end do
                    close(16)
!c
!c interpolate tempcs onto same wavelength scale as sample
!c
                    do il=1,nlambs
                        xx=slamb(il)
                        stscat(il,ncont)=ainter(templa,tempcs,xx,nlambt)
                        ssscat(il,ncont)=ainter(templa,tempscs,xx,nlambt)
                    end do
                endif
                stmon(ncont)=0
            endif
!c
!c read sample tweak factor
!c
            read(nin,*) tweak(ncont)
!c
!c calculate the sample calibration constant for first item
!c
            if(ncont.eq.1) then
                if(sgeom.eq.1) then
!c
!c area and volume of cylinder in neutron beam
!c
                    call EQUIV(sdimen1(ncont),sdimen2(ncont),300,AR,ARB,ARS)
                    area=arb
                    if((hbup-hbdown).le.sheight) then
                        vol=area*(hbup-hbdown)
                    else
                        vol=area*sheight
                    endif
                    snorm(ncont)=(srho(ncont)*vol)
                else
!c
!c thickness of flat plate sample
!
                    sthick=sdimen2(ncont)+sdimen1(ncont)
                    snorm(ncont)=(srho(ncont)*sthick)/cosd(srot(1))
                endif
            else
                snorm(ncont)=1.0
            endif
        endif
    return
    end subroutine get_sam_par

    subroutine readandprocessmint(fname,nlamb,lambscale,totalcs,sigs,siga,rho,is)

      use inputfilestrings
      use reallocation_routines
      use run_par
      use groups_routines
      use beam_routines
      use sam_par
      use write_routines

        real, dimension(*), intent(in)                  :: lambscale
        real, dimension(:), allocatable, intent(out)    :: totalcs
        real, dimension(:), allocatable                 :: qdata,dcsdata,dcserr
        real                                            :: sigs,siga,rho,dummy
        real                                            :: q1,q2
        integer                                         :: ic,is,il,nlamb,ndata,nc,iq,ierr,mq !Internal variables

        character*256 fname,fnameout,sometext

        call reallocate1d_r(totalcs,nlamb)
        !c Read the data to be processed
        write(6,101) fname(1:len_trim(fname))
101   format(/'readandprocessmint> Name of the file to be processed: ',a)
!c Open the file and read it
   open(10,file=fname,status='old',iostat=ierr)
        !If an error occurred, we simply put out the standard total cross section file
        if(ierr.ne.0) then
            fname=convertfilename(fname)
            open(10,file=fname,status='old',iostat=ierr)
            if(ierr.ne.0) then
                do il=1,nlamb
                    totalcs(il)=sigs+lambscale(il)*siga/1.798
                end do
                ncvals(is)=0
                return
            end if
        end if
!c Ignore lines with # in them
   sometext='#'
   do while(index(sometext,'#').gt.0)
            read(10,*,iostat=ierr) sometext
   end do
!c Backspace 1 line
   backspace(10)
!c Now read the data. It is assumed that the last line containing data is terminated
!c with a new line character. Otherwise the last line will be lost.
   ierr=0
   ndata=0
        mq=0
        do while(ierr.eq.0)
            ndata=ndata+1
            if(ndata.gt.mq) then
                mq=mq+10
                call reallocate1d_r(qdata,mq)
                call reallocate1d_r(dcsdata,mq)
                call reallocate1d_r(dcserr,mq)
            end if
            read(10,*,iostat=ierr) qdata(ndata),dcsdata(ndata),dcserr(ndata)
   end do
   close(10)
!c Throw away the last line in case it contains garbage
   ndata=ndata-1
   write(6,102) qdata(ndata)
102   format(/'readandprocessmint> Maximum Q found is: ',f10.5)
!c Each bin boundary will be placed midway between the existing bins
        q1=qdata(1)
   do iq=2,ndata
      q2=qdata(iq)
      qdata(iq-1)=0.5*(q1+q2)
      q1=q2
   end do
        ndata=ndata-1
   fnameout=fname(1:index(fname,'.')-1)//'.temp'
   call w_diag_file(fnameout,ndata,qdata,dcsdata,dcserr)
!c Perform integration over 2*theta
        nc=50 
!        if(2*nc.gt.mcval) nc=mcval/2
        call integsofq(nlamb,lambscale,totalcs,ndata,qdata,dcsdata,sigs,siga,nc,ncvals(is),czvals,pcs,is)
!c As a check write out the results
   fnameout=fname(1:index(fname,'.')-1)//'.pofc'
        open(10,file=fnameout,status='unknown',iostat=ierr)
        if(ierr.eq.0) then
            do il=1,nlamb
                write(10,'(/a,1x,i5,1x,i5,1x,f10.5)') '#',ncvals(is),nlamb,lambscale(il)
                do ic=1,ncvals(is)
                    write(10,'(f10.5,1x,f10.5,1x,e12.5)') lambscale(il),czvals(ic),pcs(ic,il,is)
                end do
            end do
        endif
        close(10)
        return
    end

    subroutine integsofq(nlamb,lamb,integdata,ndata,qdata,sdata,sigs,siga,nz,nzval,zval,pz,is)


      use inputfilestrings
      use reallocation_routines
      use run_par
      use groups_routines
      use beam_routines
      use sam_par
      use write_routines

        real, dimension(*), intent(in)                      :: lamb,qdata,sdata
        real, dimension(*), intent(out)                     :: integdata
        real, dimension(:), allocatable, intent(out)        :: zval
        real, dimension(:,:,:), allocatable, intent(out)     :: pz
        integer, intent(out)                                :: nzval
        real q,q1,q2,data1,data2,z,zstep,zstep2,pi,sum,grad,term
        real sigs,siga,sigs4pi
        integer is,iq,nlamb,ndata,nz,ilamb,iz
        logical found

        pi=4.0*atan(1.0)
        sigs4pi=sigs/(4*pi)
!c Step in z
        nzval=2*nz
        zstep = 1.0/float(nz)
        zstep2=0.5*zstep
        ilamb=0
        call reallocate1d_r(zval,nzval)
        call reallocate3d_r(pz,nzval,nlamb,ncont)
        do while(ilamb.lt.nlamb)
            ilamb=ilamb+1
            sum=0.0
!c Start from z = -1+0.5/nz
            z=-1.0-zstep2
            iz=0
            do while (iz.lt.nzval)
                iz=iz+1
                z=z+zstep
                zval(iz)=z
!c Convert this to Q
                q=4*pi*sqrt(0.5*(1.0-z))/lamb(ilamb)
!c Get the nearest pair of data points to this q value
                q2=0
                found=.false.
                iq=1
!c Data are assumed to be in order of increasing Q
                do while(iq.lt.ndata.and..not.found)
                    iq=iq+1
                    q2=qdata(iq)
                    found=q.le.q2
                end do
                q1=qdata(iq-1)
!c interpolate to this q value
!c For points beyond the range of the data we simply use the first or last values
!c The dcs beyond the last value is assumed zero in this case.
                if(q.gt.qdata(ndata)) then
                    term=0.0
                else if(q.lt.qdata(1)) then
                    term=sdata(1)
                else
                    grad=(sdata(iq)-sdata(iq-1))/(q2-q1)
                    term=sdata(iq-1)+(q-q1)*grad
                endif
                term=sigs4pi+term
                sum=sum+term
                pz(iz,ilamb,is)=term
!c         write(6,*) z,q,ilamb,lamb(ilamb),iq,q1,q2,grad,term,sum
            end do
            integdata(ilamb)=lamb(ilamb)*siga/1.798+sum*2*pi*zstep
        end do
        return
    end subroutine integsofq

END MODULE sam_par_routines                  
