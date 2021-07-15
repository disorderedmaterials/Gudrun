!     
! File:   gudrun_dcs.f90
! Author: aks45
!
! Created on 23 October 2013, 15:57
!

!********************************************************************************
!*
!*      Gudrun_dcs.for
!*
!*      A K Soper, January 2003
!*      Incorporates the ISIS GET routines from Freddie Akeroyd
!*      and better smooth routine for vanadium, plus several other
!*       user requested features.
!*
!*      This produces a DCS file, reading data directly from 
!       ISIS RAW files
!*
!*      It was initially adapted from Alan Soper's program DCS.for 
!*      by Piers Buchanan
!*
!*      Modified further by AKS in March 2001, to accommodate the 
!*                     'GUDRUN Input File'
!*
!*      GUI_1 version, started 15th Feb 2005 to accommodate the GudrunGUI_1 (Java)
!*       platform.
!*
!*      GUI_1 version, modification 1, started 21st March, to accommodate D4 data files 
!
! GUI_4 version to incorporate NeXus format and reallocatable arrays. Started on 23/10/2013
!*
!********************************************************************************
      PROGRAM gudrun_dcs

        use ran1_routines
        use inputfilestrings
        use reallocation_routines
        use runfactor_list
        use run_par
        use calibration_routines
        use bad_detectors
        use spec_van
        use spec_bakv
        use spec_sam
        use spec_baks
        use monitor_routines
        use groups_routines
        use beam_routines
        use local_data
        use merge_routines
        use sam_par
        use van_par
        use get_data_routines
        use deadtime_routines
        use amutshield_routines
        use sam_par_routines
        use transmission_routines
        use multiple_scattering_routines
        use mul_corr_azi
        use placzek_routines
        use normalisation_routines
        use module_sums
        use attenuation_correction_routines
        use process_routines
        use write_routines

        implicit none

!      include 'dimension.inc'
!      include 'spec_van.inc'
!      include 'spec_bakv.inc'
!      include 'spec_sam.inc'
!      include 'spec_baks.inc'
!      include 'van_par.inc'
!      include 'sam_par.inc'
!      include 'mul_corr_azi.inc'
!      include 'mul_sam_azi.inc'
!      include 'abs_corr_azi.inc'
!      include 'abs_sam_azi.inc'
!      include 'merge_arrays.inc'
!c
!c internal variables
!c
        integer i,is,is3,iff,j,k,l,m,ivar,igood      !internal indices
      integer iinst,ibeam,ivan,ibs,iloop      !indicators that data have been read
      integer debug
      integer inform,nchanout
      integer nsmoom            !number of 3-point smooths for mon and v
      integer nsample            !number of containers plus sample
      integer nfilevold            !saves no. of vanadium runs
      integer runvold            !saves current vanadium run number
      integer nfilebsold            !saves no. of sample background runs
      integer runbsold            !saves current sample background run no.
      integer, dimension(:), allocatable :: nfilesold      !saves no. of sample runs
      integer, dimension(:), allocatable :: runsold      !save run numbers of sample & containers
      integer ngrpwrt            !0 = 1 group, 1 = ngroup groups
      integer nspecwrt      !diagnostic spectrum no. 0 means no diagnostics
      integer ispec                  !spectrum counter
      integer forcecalc            !=1 forces calc. of sample corrs.
      integer mergepwr            !power used to set Q weighting for merge
      integer nbacksub            !1 = background subtraction prior to merge
      integer lenteststring            !length of teststring used to control algorithm
      integer itest               !Used to test for exponential amplitude and decay values 
      integer ilenfnameselfe      !length of filename used to read in self scattering.
      integer iref,ierr,isoftedges
      integer nbacksubwavelength  !1 = wavelength data subtraction prior to merge
      integer :: mproc   !Number of spectra to process at one time.

      real wavemin,wavemax            !minimum and maximum wavelengths to use
      real wavestep                  !step in wavelength for corrections
      real qmin,qstep,qmax,qpower            !qmin, qstep, and qmax for final merge
      real acceptance            !acceptance fraction for final merge
      real aval,bval,cval            !min and max Q for individual groups, and back. factor
      real qshell,rshell,rmax            !Width of top hat function, minimum radius for background subtraction, and maximum radius for Fourier transform
      real qwindow                  !window function for final FT. 0 for no window function
      real normalisationfactor    !Factor to multiply normalisation by prior to dividing into sample.
      real alamb,dcssave          !wavelength and temporary dcs store
      real vanlimit                 !Lowest acceptable limit on the smoothed vanadium value.
      real smoovlimit               !Controls amount of smoothing on vanadium
      real amushield                !Absorption coefficient of shielding

      character*256 fname            !name of file to read input data from
      character*256 teststring      !string used to test values in input file
      character*256 fnameselfe      !name of file containing DCS as a function of energy:
!c this will be used to subtract the self scattering for each detector
      character*20 geometry

      logical ifoundselfe,ifoundline

      write(6,*) 'Let us begin: Gudrun_dcs NeXus dynamic allocation version GUI_4: November 2013'

!c initialise random number generator

!      call initran1(10)

!c initialise no. of sample plus containers

      nsample=0

!c indicates that instrument data have been read in

      iinst=0

!c indicates that beam data have been read in

      ibeam=0

!c indicates that vanadium data have been read in

      ivan=0

!c indicates that sample background data have been read in

      ibs=0

!c initialise indicators as to whether sample or container data have been read
!c in previously

!      do is=1,mcont
!            nfilesold(is)=0
!            runsold(is)=-1
!      end do

!c open input file

!c      write(6,*) 'Type name of file to read input data from'
!     read(5,'(a)') fname
! Check to see if an input file has been specified.
      if(iargc().gt.0) then
         call getarg(1,fname)
      else
         fname='gudrun_dcs.dat'
      endif
      open(15,file=fname,status='old',iostat=ierr)
      if(ierr.ne.0) then

         close(15)
         write(6,'(2(a))') 'Cannot open ',fname(1:index(fname," "))
         stop

      else
        write(6,*) 'gudrun_dcs.dat opened correctly'
      
      endif

!c Ignore empty lines or lines that begin with spaces. The smallest allowable
!c character is a blank - ASCII decimal 32. Anything less than or equal to this should be ignored.

    nchar=0
    do while(nchar.eq.0.and.ierr.eq.0)
        read(15,'(a)',iostat=ierr) line
        nchar=len_trim(line)
    end do
    if(ierr.ne.0) then
        write(6,'(a,1x,i4)') 'gudrun_dcs> ierr = ',ierr
        stop
    else
        nchar=len_trim(line)
    end if
!    write(6,*) line(1:nchar)

!c The first valid line read must contain the strings which signify the spacing between
!c values (spc2), the spacing between the last value and the following comment (spc5), and the 
!c pathseparator character. These values are delineated by apostrophes as the number of 
!c spaces is important for these values. Hence a special reading method is required.
!c
!c Find the first and second occurrence of a '
!c
    index1=index(line(1:nchar),'''')
    index2=index1+index(line(index1+1:nchar),'''')
    spc2=' '
    lenspc2=1
    index1=index2+index(line(index2+1:nchar),'''')
    index2=index1+index(line(index1+1:nchar),'''')
    spc5=line(index1+1:index2-1)
    lenspc5=index2-index1-1
    index1=index2+index(line(index2+1:nchar),'''')
    index2=index1+index(line(index1+1:nchar),'''')
    pathseparator=line(index1+1:index2-1)
    lenpathseparator=index2-index1-1
    write(6,*) spc2(1:lenspc2),';',spc5(1:lenspc5),';',pathseparator(1:lenpathseparator)

    call getaline(15)
    ifoundline=nwords.gt.0
    if(ifoundline) ifoundline=line(ncf(1):ncl(1)).ne.'END'

    do while (ifoundline)

!c setup the test string. Currently the name or numbers after 'SAMPLE' or 'CONTAINER' are not used

        if((line(ncf(1):ncl(1)).eq.'CONTAINER') &
        .or.(line(ncf(1):ncl(1)).eq.'SAMPLE' &
        .and.line(ncf(1):ncl(2)).ne.'SAMPLE BACKGROUND')) then
            teststring=line(ncf(1):ncl(1))
            lenteststring=ncl(1)-ncf(1)+1
        else
            teststring=line(ncf(1):ncl(nwords))
            lenteststring=ncl(nwords)-ncf(1)+1
        endif

!c INSTRUMENT section

        if(teststring(1:lenteststring).eq.'INSTRUMENT') then

            write(6,601) teststring(1:lenteststring)
601         format(/' gudrun_dcs> Got to: ',a)

            call getaline(15)
            inst=line(ncf(1):ncl(1))
            write(6,*) inst(1:index(inst,' ')-1)
            call getaline(15)
            directinp=line(ncf(1):ncl(nwords))
            lendirectinp=ncl(nwords)-ncf(1)+1
            write(6,999) directinp(1:lendirectinp)
999         format(1x,a)
            call getaline(15)
            direct=line(ncf(1):ncl(nwords))
            lendirect=ncl(nwords)-ncf(1)+1
            write(6,999) direct(1:lendirect)

!c Read a line where the Raw file extension is given.

            call getaline(15)
            ext=line(1:3)


!c Read name of detector calibration file. Use * if using default calibration

            call getaline(15)
            fnamed=line(ncf(1):ncl(nwords))
            write(6,999) fnamed(1:len_trim(fnamed))
!c
!c read User Table Column number for phi values
!c
            read(15,*) phiutno
            write(6,*) phiutno
!c
!c Read name of detector groups file. Use * to use default (groups_def.dat)
!c
            call getaline(15)
            fname_grps=line(ncf(1):ncl(nwords))
            ilenfname_grps=ncl(nwords)-ncf(1)+1
            write(6,999) fname_grps(1:ilenfname_grps)
!c
!c Read name of file containing detector, module and data acquisition dead times
!c Use * to use no deadtime corrections. A file called 'deadtime.cor' will be written
!c during the run which contains a list of the module numbers and the deadtimes that have
!c been used in the program.
!c
            call getaline(15)
            deadtimefile=line(ncf(1):ncl(nwords))
            write(6,999) deadtimefile(1:index(deadtimefile,'  ')-1)
!c
!c Read number of spectra of incident and transmitted beam monitor, 
!c spectrum numbers,
!c and fraction of monitor count for time independent background
!c
            call getaline(15)
            ninmon=nwords
            write(6,*) ninmon
            if(nwords.gt.0) then
                call reallocate1d_i(incid_mon,nwords)
                read(line(ncf(1):ncl(nwords)),*) (incid_mon(i),i=1,nwords)
            end if
            write(6,*) incid_mon

!c
!c Read the input wavelength range for monitor. Use 0 0 to divide channel by channel.
!c
            read(15,*) monwave1,monwave2

!c Get the transmission monitor spectrum numbers

            call getaline(15)
            ntrmon=nwords
            if(nwords.gt.0) then
                call reallocate1d_i(trans_mon,nwords)
                read(line(ncf(1):ncl(nwords)),*) (trans_mon(i),i=1,nwords)
            end if

!c Remove any negative spectrum numbers from this list.

            i=0
            do while(i.lt.ntrmon)
               i=i+1
               if(trans_mon(i).lt.0) then
                  ntrmon=ntrmon-1
                  j=i
                  do while(j.lt.ntrmon)
                     trans_mon(j)=trans_mon(j+1)
                     j=j+1
                  end do
               endif
            end do
            write(6,*) ntrmon
            write(6,*) trans_mon

!c Read the quiet count constants for incident and transmission monitors

            read(15,*) fraci
            read(15,*) fract
!c
!c Read first and last channel numbers to check for spikes (0,0 to use all)
!c
            read(15,*) nchfir,nchlas
!c
!c Read flagsp, which decides whether to exclude spikes and steps, and no. of std devs
!c used to decide on spike or not
!c
            read(15,*) spikedev
            flagsp=0
            if(spikedev.gt.0.0) then
                  flagsp=1
            else
                  flagsp=0
            endif
            spikedev=spikedev*spikedev
!c
!c Read range of incident wavelengths to be used and step size for corrections
!c
            read(15,*) wavemin,wavemax,wavestep
!c
!c Number of smoothings for monitor and vanadium
!c
            read(15,*) nsmoom
!c
!c Qmin, qmax, and qstep for final merged data. These values will be overridden
!c for individual groups.
!c Negative value of qstep means logarithmic binning above qmin. qmin is the smallest q value written
!c to the merged data file
!c
            read(15,*) qmin,qmax,qstep
            !Set the default qmax for groups
            qmax_for_groups=qmax
!c
!c set the q ranges of particular groups. If the Q-range for a particular
!c group is not defined here, it will be left at its default value
!c
!c Use a Qmin or Qmax set to 0 to get the default value
!c Use a Qmax below or equal to Qmin to ignore this group in the final merge
!c
!c Use 0 0 0 0 to terminate this input, or to indicate no groups have special 
!c Q-ranges
!c
            ivar=0
            i=1
            ngroup_in=0
            oldsize=0
            newsize=oldsize+5
            do while (i.gt.0)
                  read(15,*) i,aval,bval,cval
                  write(6,*) i,aval,bval,cval
                  if(i.gt.0) then
                      ngroup_in=ngroup_in+1
                      if(ngroup_in.gt.oldsize) then
                          newsize=oldsize+5
                          call reallocate1d_i(igroup_in,newsize)
                          call reallocate1d_r(qmingrp_in,newsize)
                          call reallocate1d_r(qmaxgrp_in,newsize)
                          call reallocate1d_r(background_factor_in,newsize)
                          oldsize=newsize
                      end if
                      igroup_in(ngroup_in)=i
                      qmingrp_in(ngroup_in)=aval
                      qmaxgrp_in(ngroup_in)=bval
                      background_factor_in(ngroup_in)=cval
                      write(6,*) igroup_in(ngroup_in),qmingrp_in(ngroup_in),qmaxgrp_in(ngroup_in),background_factor_in(ngroup_in)
                  endif
            end do
!c
!c Integer to define whether 1 group or ngroup groups are output at the end
!c
            ngrpwrt=1
!c
!c number of spectra to process at one time
!c
            nproc=mproc

!c the final merge of the groups is accomplished by comparing the scattering
!c level of each group with the expected scattering level (calculated from the
!c input sample parameters). If a group is outside the expected level by 
!c +-acceptance*(expected level), it is excluded from the merge. 
!c An acceptance of 1.0 means that all groups will be accepted.
!c
            read(15,*) acceptance
            read(15,*) mergepwr
!c
!c Decide whether to subtract a background from each group prior to merge or 
!c not. 1 means yes subtract a background prior to merge, 0 means no.
!c 
!c nweighterr = 1 means to use error bars in merge, 0 means uniform weights
!c
            read(15,*) nbacksub
            read(15,*) nweighterr
!c
!c make sure nweighterr is either 1 or 0. Default to 0 if not 1.
!c
            if(nweighterr.gt.2.or.nweighterr.lt.0) nweighterr=0
!c
!c read the incident flight path. If zero thisc RAW data file.
!c
            read(15,*) incidentfp

!c For reactor instrument need to define the incident wavelength

            if(inst(1:(index(inst,' ')-1)).eq.'D4C') then
               read(15,*) incidentwl,zeroangleoffset
            endif

!c Spectrum number for diagnostic files: -ve, 0 or >nspec means no diagnostic
!c files will be written. Normally this should be 0, unless malfunctioning of the
!c programme is suspected.

            read(15,*) nspecwrt

!c Get the name of the file which carries the neutron scattering lengths

            call getaline(15)
            neutronparfilename=line(ncf(1):ncl(nwords))
            write(6,999) line(ncf(1):ncl(nwords))

!c Get the integer which indicates what units for the final output.
!c 1 = Q [1/A], 2 = d-space [A], 3 = wavelength [A], 4 = energy [meV], 5 = TOF [musec]

            read(15,*) outputunitstype
            write(6,*) outputunitstype
            if(outputunitstype.lt.1.or.outputunitstype.gt.5) &
     outputunitstype=1

!c Get the integer that decides whether to subtract the wavelength data prior to merge

            read(15,*) nbacksubwavelength

!c Get the default folder for calibration, groups and other files to exist if they are not specified
!c on the input line.

            call getaline(15)
            gudrunstartupfolder=line(ncf(1):ncl(nwords))
            ilengudrunstartupfolder=ncl(nwords)-ncf(1)+1

!c Get the folder containing the GUI startup file and where instrument specific files are
!c stored

            call getaline(15)
            gudrunstartupfilefolder=line(ncf(1):ncl(nwords))
            ilengudrunstartupfilefolder=ncl(nwords)-ncf(1)+1

!c Get the power to be used to increment the step size

            read(15,*) qpower

!c This should be zero if qstep > 0

            if(qstep.gt.0.0) then            

               qpower=0.0

            endif

!c Read whether to use hard group edges (0) or soft group edges (1)

            isoftgroupedges=0
            read(15,*) isoftgroupedges
            softgroupedges=isoftgroupedges.ne.0 !N.B. to get soft group edges you need isoftgroupedges .ne. 0

            !c For nexus files we need to read in the list of paths needed to find particular data.
!c The filename will be called 'instrument'_nexus.txt

            if(ext.eq.'nxs'.or.ext.eq.'NXS') then
!Get the name of the NeXus information file
                call getaline(15)
                fname=line(ncf(1):ncl(nwords))
                call get_pathnames_nxs(fname)
            endif

!c get a list of run number factors. This is used to multiply the raw data by a factor BEFORE
!c it is aggregated with other runs. The factors are stored in a file called runfactor_list.dat.
!c If a particular run is not listed then the factor for that run is assumed to be 1.0. If the file
!c runfactor_list.dat does not exist, then the factors for all runs are assumed to be 1.0.

            call get_runfactor_list

!c set indicator

            iinst=1
            ibeam=0
            ivan=0
            ibs=0
            nsample=0

!c BEAM section

        else if(teststring(1:lenteststring).eq.'BEAM') then

            write(6,601) teststring(1:lenteststring)

!c
!c Read and generate beam parameters for this experiment
!c
!c Sample (and Vanadium) geometry

            call getaline(15)
            if(line(ncf(1):ncl(nwords)).eq.'FLATPLATE') then
                  vgeom=2
                  sgeom=2
            else
                  vgeom=1
                  sgeom=1
            endif

            call init_beam(15)

!c Get the sample dependent background factor

            call getaline(15)
            sampledependentbackgroundfactor=0.0
            if(nwords.gt.0) then
               read(line(ncf(1):ncl(1)),*) sampledependentbackgroundfactor
            endif
            subtractsampledependentbackground=sampledependentbackgroundfactor.ne.0.0

!c Get absorption coefficient for the shielding - this is assumed to affect the
!c multiple scattering (and sample dependent background if present).
 
            call getaline(15)
            amushield=0.0
            if(nwords.gt.0) then
               read(line(ncf(1):ncl(1)),*) amushield
            endif

            ibeam=1
            ivan=0
            ibs=0
            nsample=0

!c NORMALISATION section

      else if(teststring(1:lenteststring).eq.'VANADIUM' &
      .or.teststring(1:lenteststring).eq.'NORMALISATION') then

            write(6,601) teststring(1:lenteststring)

!c cannot define vanadium until instrument has been defined
!c
            if(iinst.ne.1.or.ibeam.ne.1) then
                  write(6,*) 'gudrun_dcs> Cannot define VANADIUM before INSTRUMENT and BEAM'
                  stop
            endif
!c
!c Number of vanadium runs to add, and period number for vanadium
!c
            read(15,*) nfilev,nperv
            call reallocate1d_c(runv,len(runv),nfilev)
!c
!c Vanadium run filenames
!c
            do i=1,nfilev
               call getaline(15)
               runv(i)=line(ncf(1):ncl(nwords))
            end do
!c
!c Number of vanadium background runs to add, and period number for vanadium
!c background
!c
            read(15,*) nfilebv,nperbv
            call reallocate1d_c(runbv,len(runbv),nfilebv)
!c
!c Vanadium background run numbers
!c
            do i=1,nfilebv
               call getaline(15)
               runbv(i)=line(ncf(1):ncl(nwords))
            end do
!c
!c Initialise run parameters
!c
            call get_run_par(runv,nfilev)
            write(6,*) 'gudrun_dcs> Got run parameters'
!c
!c read the default detector, module and data acquisition deadtimes, 
!c and the detector and module deadtimes for individual modules
!c
            call get_deadtimes()
            write(6,*) 'gudrun_dcs> Got deadtimes'

!c for each of the vanadium runs get the time channel boundaries or check 
!c that the number of time channels is the same as the first run

!c Initially must set the number of time channels to zero

            nchanv=0
            nchans=0
            do i=1,nfilev
                  call get_tcbv(runv(i))
            end do

!c for each of the vanadium background runs check that the
!c number of time channels is the same as the vanadium

            do i=1,nfilebv
                  call get_tcbv(runbv(i))
            end do
            write(6,*) 'gudrun_dcs> Got time channel boundaries'
!c
!c Initialise calibration parameters
!c
            call get_calibration(runv,nfilev)
            write(6,*) 'gudrun_dcs> Got calibration'
 
!c Initialise the spectral groups, the bad spectrum index, and the spike counter. Always look
!c for spec.bad and spike.dat (ignore=0, generated by purge_det).
!c
            call get_groups(0,wavemin,wavemax,outputunitstype)
            write(6,*) 'gudrun_dcs> Got groups'
!c            write(6,*) (ibad(i),i=1,20)

!c Get the list of individual group amushields from amutshield.dat - any groups not specified are
!c given the value of amushield. This file can only be read once the groups have been defined.
 
!c If amushield is zero, then this zero value is used for all groups.

            call get_amutshield(amushield)
!c
!c get the vanadium and vanadium background monitors
!c
            nspecproc=ninmon
            call reallocate1d_i(specproc,nspecproc)
            do is=1,nspecproc
                  specproc(is)=incid_mon(is)
            end do
!c
!c vanadium monitor
!c
            call get_mon(1,1)
!c
!c form the vanadium monitor sums
!c
!            call get_mon_period_sum(1)
!c
!c smooth vanadium monitor
!c
            call form_smoo_mon(1,nsmoom,wavemin,wavemax)
!           write(6,*) 'Smoothed vanadium monitor'
!c
!c vanadium background monitor
!c
            call get_mon(2,1)
!c
!c form the vanadium background monitor sums
!c
!            call get_mon_period_sum(2)
!c
!c smooth vanadium background monitor
!c
            call form_smoo_mon(2,nsmoom,wavemin,wavemax)
!            write(6,*) 'Smoothed vanadium background monitor'
!c
!c now the transmission monitor counts if present
!c
            if(ntrmon.gt.0) then
                nspecproc=ntrmon
                call reallocate1d_i(specproc,nspecproc)
                do is=1,nspecproc
                    specproc(is)=trans_mon(is)
                end do


!c
!c vanadium monitor
!c
                call get_mon(1,2)
!
! vanadium background monitor
!c
                call get_mon(2,2)
                write(6,*) 'gudrun_dcs> Got vanadium tranmission monitors'
            endif
!c
!c decide whether to force the calculation of the corrections (forcecalc=1)
!c otherwise they will be read from a file, unless the file does not exist
!c
            read(15,*) forcecalc
!c
!c Read vanadium sample parameters
!c
            call get_sam_par(15,0,runv(1),nperv,wavemin,wavemax,wavestep)
            write(6,*) 'gudrun_dcs> Got vanadium parameters'

            write(6,*) 'gudrun_dcs> NORMALISATION geometry = ',vgeom
!c Get vanadium transmission
            call get_trans(1,wavemin,wavemax,wavestep)
!c
!c if the transmission monitor is to be used to estimate the total cross section
!c of the sample, then do it now
!c
            if(vtmon.ne.0) then
                  call get_trans_cs(1)
                  write(6,*) 'gudrun_dcs> Got vanadium c/s from transmission monitor'
            endif
!c
!c Get the vanadium multiple scattering
!c
            call get_mul_corr(1,runv(1),nperv,forcecalc)
            write(6,*) 'gudrun_dcs> Got vanadium total scattering correction'
!c
!c Get the vanadium Placzek correction
!c
            call get_placzek(1,runv(1),nperv,forcecalc)
            write(6,*) 'gudrun_dcs> Got vanadium Placzek correction'
!c
!c Vanadium calibration .DCS filename. (* means it will not attempt to read a
!c vanadium .DCS file.)
!c
            call getaline(15)
            vdcsname=line(ncf(1):ncl(nwords))
            write(6,999) vdcsname(1:len_trim(vdcsname))
            vsmear=0.0
            vrmin=0.0
!c
!c get the vanadium differential c/s file if required
!c
            if(index(vdcsname,'*').le.0) call vdcs_read()
            write(6,*) 'gudrun_dcs> Got vanadium DCS data'
!c Get the lowest acceptable value for the smoothed vanadium - below this value and the detector is rejected.
            call getaline(15)
            vanlimit=0.01
            if(nwords.gt.0) then
               read(line(ncf(1):ncl(1)),*) vanlimit
            end if
!c Get the degree of smoothing on the vanadium
            call getaline(15)
            smoovlimit=1.0
            if(nwords.gt.0) then
               read(line(ncf(1):ncl(1)),*) smoovlimit
            end if
! Get vanadium signal to background acceptance ratio
            call getaline(15)
            vanbackfraction=0.0
            if(nwords.gt.0) then
               read(line(ncf(1):ncl(1)),*) vanbackfraction
            end if
!c Get the vanadium number of smooths for each spectrum if they exist
            call get_nsmoov(smoovlimit)
!      write(6,*) 'gudrun_dcs> 784'
!c
!c D: initialize the merge arrays
!c
            call init_merge(qmin,qstep,qmax,qpower,ngrpwrt)
!      write(6,*) 'gudrun_dcs> 789'
!c Initialise vanadium dcs to be used in vanadium normalisation
            call vdcs_bragg()
!      write(6,*) 'gudrun_dcs> 792'

!
!c get module sums for vanadium and vanadium background - these will be used in the deadtime
!c correction if requested
!c
            call get_module_sum(1)
            write(6,*) 'gudrun_dcs> Got module sums for vanadium'
            call get_module_sum(2)
            write(6,*) 'gudrun_dcs> Got module sums for vanadium background'

! Open a file to output the vanadium smoothing parameters

            nchanout=31
            ndettested=0
            ndetrejected=0
            open(nchanout,file='vansmo.par', status='unknown',iostat=ierr)
!
! Process the vanadium data
!
            ispec=0
            mproc=nspec
            call reallocate1d_i(specproc,mproc)
!            do while (ispec.le.nspec)
!Get the valid spectra
                nspecproc=0
                do while(nspecproc.lt.mproc.and.ispec.lt.nspec)
                    ispec=ispec+1
                    if(ibad(ispec).eq.0.and.igrp(ispec).gt.0.and.igrp(ispec).le.ngroup) then
                        nspecproc=nspecproc+1
                        specproc(nspecproc)=ispec
                    end if
                end do
                nspecf=specproc(1)
                nspecl=specproc(nspecproc)
                write(6,*) 'gudrun_dcs> Normalisation: first spec., last spec., nspecproc ',nspecf,nspecl,nspecproc
!c
!c Form sums of vanadium counts
!c
                call get_sum(1,nspecwrt)
                write(6,*) 'gudrun_dcs> Got vanadium counts'
!c
!c Divide vanadium data by smoothed monitor
!c
                call divide_by_mon(1,nspecwrt)
                write(6,*) 'gudrun_dcs> Divided vanadium by monitor'
!c
!c Form sums of vanadium background counts
!c 
                call get_sum(2,nspecwrt)
                write(6,*) 'gudrun_dcs> Got vanadium background counts'
!c
!c Divide vanadium background data by smoothed monitor
!c
                call divide_by_mon(2,nspecwrt)
                write(6,*) 'gudrun_dcs> Divided vanadium background by monitor'
!c
!c Subtract vanadium background from vanadium
!c
                call subtract_bak(1,wavemin,wavemax,nspecwrt)
                write(6,*) 'gudrun_dcs> Subtracted vanadium background from vanadium'
!c
!c Do vanadium  total scattering correction
!c
                call do_mul_corr(1,nspecwrt)
                write(6,*) 'gudrun_dcs> Renormalised vanadium data to total scattering'
!c
!c ... and smooth result
!c
                call form_smoo_van(smoovlimit,nspecwrt,vanlimit,nchanout)
                write(6,*) 'gudrun_dcs> Finished smoothing vanadium'
!                stop
!                ispec=nspecl+1
!            end do
            call save_nsmoov(smoovlimit)
            !Write the number of detectors tested and number rejected
            if(ndettested.gt.0) then
                write(nchanout,567) ndetrejected,ndettested,real(ndetrejected)/real(ndettested)
                write(6,567) ndetrejected,ndettested,real(ndetrejected)/real(ndettested)
567             format(' gudrun_dcs> Number of spectra rejected' &
                ,' by vanadium smooth = ',i5 &
                ,/' gudrun_dcs> Number of spectra tested by vanadium smooth = ' &
                ,i5,/' gudrun_dcs> Rejection ratio = ',f10.5)
            end if
            write(6,*) 'gudrun_dcs> ',nspecprocgood,' good spectra left'
            close(nchanout)
!c
!c signal that the vanadium has been done
!c
            ivan=1
            ibs=0
            nsample=0

!c SAMPLE BACKGROUND section

        else if(teststring(1:lenteststring).eq.'SAMPLE BACKGROUND') then

            write(6,601) teststring(1:lenteststring)

!c
!c Number of sample background runs to add
!c
            read(15,*) nfilebs,nperbs
            call reallocate1d_c(runbs,len(runbs),nfilebs)
!c
!c Sample background run filenames
!c
            do i=1,nfilebs
                call getaline(15)
                runbs(i)=line(ncf(1):ncl(nwords))
            end do
!c
!c for each of these runs get the time channel boundaries and check that the
!c number of time channels is the same as the other samples
!c
            do i=1,nfilebs
                  call get_tcbs(runbs(i))
            end do
            write(6,*) 'gudrun_dcs> Got sample background time channel boundaries'
!c A: get the sample background incident monitor counts and smooth if
!c requested
!c
            nspecproc=ninmon
            call reallocate1d_i(specproc,nspecproc)
            do is=1,nspecproc
                  specproc(is)=incid_mon(is)
            end do
!c
!c sample background monitor
!c
            call get_mon(3,1)
            call form_smoo_mon(3,nsmoom,wavemin,wavemax)
            write(6,*) 'gudrun_dcs> Smoothed sample background monitor'
!c
!c now the transmission monitor counts if present
!c
            if(ntrmon.gt.0) then
                nspecproc=ntrmon
                call reallocate1d_i(specproc,nspecproc)
                do is=1,nspecproc
                    specproc(is)=trans_mon(is)
                end do
                call get_mon(3,2)
            end if
!Get module sums for sample background
            call get_module_sum(3)
            write(6,*) 'gudrun_dcs> Got module sums for sample background'
!Get the sample background sums and normalise to monitor for all the good spectra
            !The good spectra are set up by the vanadium smooth
            igood=0
            mproc=nspecprocgood
            call reallocate1d_i(specproc,mproc)
!            do while (igood.lt.nspecprocgood)
                nspecproc=0
                do while(nspecproc.lt.mproc.and.igood.lt.nspecprocgood)
                    igood=igood+1
                    nspecproc=nspecproc+1
                    specproc(nspecproc)=specprocgood(igood)
                end do
                nspecf=specproc(1)
                nspecl=specproc(nspecproc)
                write(6,*) 'gudrun_dcs> Sample Background: first spec., last spec., nspecproc ',nspecf,nspecl,nspecproc
!c
!c Form sums of sample background counts
!c 
                call get_sum(3,nspecwrt)
                write(6,*) 'gudrun_dcs> Got sample background counts'
!c
!c Divide sample background data by smoothed monitor
!c
                call divide_by_mon(3,nspecwrt)
                write(6,*) 'Normalised sample background to monitor'
!            end do
!c
!c set indicator
!c
            ibs=1
            nsample=0
            
!c SAMPLE section

      else if(teststring(1:lenteststring).eq.'SAMPLE') then

            write(6,601) teststring(1:lenteststring)

!c
!c initialise sample counter
!c
            is=1
            call reallocate1d_i(nfiles,is)
            call reallocate1d_i(npers,is)
            call reallocate1d_r(sampleenvscatfrac,is)
            call reallocate1d_r(sampleenvamut,is)
!c
!c Initialise the number of cos theta values read for DCS files to be used in m.s. calculation. A value > 0
!c indicates at least one of the sample or containers has a DCS file defined. This will affect how the
!c the multiple scattering calculation is performed.

            ncval=0
!c
!c Set up run numbers for sample
!c
            read(15,*) nfiles(is),npers(is)
            call reallocate2d_c(runs,len(runs),nfiles(is),is)
!c
!c Sample run filenames
!c
            do i=1,nfiles(is)
               call getaline(15)
               runs(i,is)=line(ncf(1):ncl(nwords))
            end do
            write(6,*) 'gudrun_dcs> Got sample run files'
!c
!c get title of this run
!c
            call get_title(runs(1,1))
            write(6,*) 'gudrun_dcs> Got sample title'
!c
!c for each of these runs get the time channel boundaries and check that the
!c number of time channels is the same as the other samples
!c
            do i=1,nfiles(is)
                  call get_tcbs(runs(i,is))
            end do
            write(6,*) 'gudrun_dcs> Got sample time channel boundaries'
!c
!c decide whether to force the calculation of the corrections (forcecalc=1)
!c else they will be read from a file, unless the file does not exist
!c
            read(15,*) forcecalc
!c
!c read sample parameters
!c
            call get_sam_par(15,1,runs(1,is),npers(is),wavemin,wavemax,wavestep)
            write(6,*) 'gudrun_dcs> Got sample parameters'
!c
!c Read an additional Delta-Q (for broadening the
!c differential cross section), r-min (for subtracting the background), and 
!c Qwindow for the window function to be applied to the data prior to Fourier
!c transfomr. If Qwindow is zero no window function is applied prior to F.T.
!c
            read(15,*) qshell
            read(15,*) rshell
            read(15,*) qwindow
!c
!c save current number of sample plus containers
!c
            nsample=is
!c
!c get the minimum and maximum wavelengths for each
!c resonance, 1 pair per line, if any
!c
            nreson=0
            aval=1
            bval=1
            do while (aval.gt.0.0.and.bval.gt.0.0)
                  read(15,*) aval,bval
                  if(aval.gt.0.0.and.bval.gt.0) then
                     nreson=nreson+1
                     call reallocate1d_r(wavresf,nreson)
                     call reallocate1d_r(wavresl,nreson)
                     wavresf(nreson)=aval
                     wavresl(nreson)=bval
                  endif
            end do
            write(6,*) 'gudrun_dcs> Got sample resonances'
!c
!c write out resonances if present
!c
            if(nreson.gt.0) then
               write(6,*)
               write(6,*) 'No. of resonances specified: ',nreson
               write(6,*) 'Resonance wavelength ranges:-'
               do i=1,nreson
                  write(6,*) wavresf(i),wavresl(i)
               end do
               write(6,*)
            endif

!c Read the stretched exponential parameters if present

            itest=0
            nexpon=0
            do while (itest.eq.0.and.nexpon.lt.5)
               call getaline(15)
               itest=index(line(ncf(1):ncl(nwords)),'*')
               if(itest.eq.0) then
                  nexpon=nexpon+1
                  call reallocate1d_r(exponamp,nexpon)
                  call reallocate1d_r(expondecay,nexpon)
                  read(line(ncf(1):ncl(nwords)),*) exponamp(nexpon),expondecay(nexpon)
               endif
            end do
            write(6,*) 'gudrun_dcs> Got sample exponential parameters'

!c Read the normalisation factor: this is to compensate if the sample and normalisation sample (vanadium normally) do
!c not occupy the same area of beam. Normally this will only apply to flat plate samples.

            read(15,*) normalisationfactor

!c Get name of file containing self dcs as a function of wavelength

            call getaline(15)
            fnameselfe=line(ncf(1):ncl(nwords))
            ilenfnameselfe=ncl(nwords)-ncf(1)+1
            write(6,*) ilenfnameselfe,line(ncf(1):ncl(nwords))

!c Read the type of normalisation required on the final merged DCS data. 0 = nothing, 1 = <f>^2, 2 = <f^2>

            read(15,*) normalisationtype

!c Read the maximum radius for the final transformation

            read(15,*) rmax

!c Read whether output units are barns/atom/sr or cm**-1/sr

            read(15,*) itest
            dofr = .false.
            if (itest.eq.1) dofr = .true.

!c Set broadening power on g(r): e.g. 0 = constant, 0.5 = sqrt(r) (recommended), 1 = r

            read(15,*) grbroad

!c Set step in radius for final g(r)

            read(15,*) finalgrstep
            logrbinning = .false.
            if (finalgrstep.lt.0.0) logrbinning = .true.
            finalgrstep=abs(finalgrstep)

!c Read the data if present and required

            ifoundselfe=.false.
            if(ilenfnameselfe.gt.0.and.nbacksubwavelength.gt.0) then
               if(index(fnameselfe(1:ilenfnameselfe),'*').le.0) then

                  ifoundselfe=.true.
                  call dcs_read(fnameselfe,ilenfnameselfe)

               endif
            endif
            write(6,*) 'gudrun_dcs> Read sample dcs file'

!c Ignore the next line which is only used by the GUI

            call getaline(15)

!c Get the sample environment factors used to compensate for different attenuation and scattering in different containers

            call getaline(15)
            if(nwords.gt.1) then
               read(line(ncf(1):ncl(2)),*) sampleenvscatfrac(is),sampleenvamut(is)
            else
               sampleenvscatfrac(is)=1.0
               sampleenvamut(is)=0.0
            endif

!c CONTAINER section

      else if(teststring(1:lenteststring).eq.'CONTAINER') then

            write(6,601) teststring(1:lenteststring)

!c
!c increment sample counter
!c
            is=is+1
            call reallocate1d_i(nfiles,is)
            call reallocate1d_i(npers,is)
            call reallocate1d_r(sampleenvscatfrac,is)
            call reallocate1d_r(sampleenvamut,is)
!c
!1c read number of container runs to add and period number for container
!c
            read(15,*) nfiles(is),npers(is)
            !When allocating the space for this container, make sure not to REDUCE the current first dimension of runs!
            call reallocate2d_c(runs,len(runs),max(nfiles(is),size(runs,1)),is)
!c
!c Container run filenames
!c
            do i=1,nfiles(is)
               call getaline(15)
               runs(i,is)=line(ncf(1):ncl(nwords))
            end do
!c
!c for each of these runs get the time channel boundaries and check that the
!c number of time channels is the same as the other samples
!c
            do i=1,nfiles(is)
                  call get_tcbs(runs(i,is))
            end do
!c
!c read container parameters
!c
            call get_sam_par(15,2,runs(1,is),npers(is),wavemin,wavemax,wavestep)

!c Get the sample environment factors used to compensate for different attenuation and scattering in different containers

            call getaline(15)
            if(nwords.gt.1) then
               read(line(ncf(1):ncl(2)),*) sampleenvscatfrac(is),sampleenvamut(is)
            else
               sampleenvscatfrac(is)=1.0
               sampleenvamut(is)=0.0
            endif
!c
!c save current number of sample plus containers
!c
            nsample=is
!c
!c O.K. begin processing the data
!c
        else if(teststring(1:lenteststring).eq.'GO'.or.teststring(1:lenteststring).eq.'CROSS_SECTION') then
!c
!c first check that all the required parameters have been defined
!c
            if(iinst.eq.0) then
                write(6,*) 'Instrument parameters need to be defined'
                call endProgram(smoovlimit,amushield)
            endif
            if(ibeam.eq.0) then
                write(6,*) 'Beam parameters need to be defined'
                call endProgram(smoovlimit,amushield)
            endif
            if(ivan.eq.0) then
                write(6,*) 'Vanadium parameters need to be defined'
                call endProgram(smoovlimit,amushield)
            endif
            if(ibs.eq.0) then
                write(6,*) 'Sample background needs to be defined'
                call endProgram(smoovlimit,amushield)
            endif
            if(nsample.lt.1) then
                write(6,*) 'No. of sample plus containers less than 1'
                call endProgram(smoovlimit,amushield)
            else
                write(6,*) 'Number of sample plus containers = ',nsample
            endif
!c
!c D: initialize the merge arrays - again!
!c
            call init_merge(qmin,qstep,qmax,qpower,ngrpwrt)
!c Initialise the vanadium and sample ratio stores
            call reallocate1d_r(vancor,nspec)
            call reallocate1d_r(vancore,nspec)
            call reallocate2d_r(samcor,nspec,ncont)
            call reallocate2d_r(samcore,nspec,ncont)
            do is=1,nspec
               vancor(is)=0.0
               vancore(is)=0.0
               do i=1,ncont
                  samcor(is,i)=0.0
                  samcore(is,i)=0.0
               end do
            end do
!c
!c now process the data
!c
!c A: get the sample and container incident monitor counts and smooth if
!c requested
!c
            nspecproc=ninmon
            do is=1,nspecproc
                  specproc(is)=incid_mon(is)
            end do
!c
!c sample and container monitors
!c
            do is=1,nsample
                is3=is+3
!c
!c form sample or container monitor sums
!c
                call get_mon(is3,1)
                call form_smoo_mon(is3,nsmoom,wavemin,wavemax)
                write(6,*) 'gudrun_dcs> Smoothed sample/container monitor ',is
!c
!c get module sums for this container
!c
                call get_module_sum(is+3)
                write(6,*) 'gudrun_dcs> Got module sums for sample/container ',is
            end do
!c
!c B: get the transmission monitor counts if present
!c
            if(ntrmon.gt.0) then
                nspecproc=ntrmon
                do is=1,nspecproc
                    specproc(is)=trans_mon(is)
                end do
                do is=1,nsample
                    is3=is+3
                         call get_mon(is3,2)
                end do
                write(6,*) 'gudrun_dcs> Got sample transmission monitor data'
            endif
!c
!c if the transmission monitor is to be used to estimate the total cross section
!c of the sample, then do it now.  This is done starting with the outermost
!c container as requested since for cylinders the result for the inner containers
!c will depend on the values assumed for the inner containers
!c
            do is=1,nsample
                i=nsample-is+1
                is3=i+3
                call get_trans(is3,wavemin,wavemax,wavestep)
                write(6,*) 'gudrun_dcs> Got transmission for sample ',i
                if(stmon(i).ne.0) then
                    call get_trans_cs(is3)
                    write(6,*) 'gudrun_dcs> Got c/s from transmission monitor for sample ',i
                endif
            end do
            if(teststring(1:lenteststring).ne.'CROSS_SECTION') then
!c
!c C: read in or calculate the sample and container corrections
!c
                do is=1,nsample
                    if(is.eq.1) then
                        write(6,*) 'gudrun_dcs> For sample'
                    else
                        write(6,*) 'gudrun_dcs> For container ',is-1
                    endif
!c
!c shifted counter
!c
                    is3=is+3
!c
!c Get Multiple scattering data for the sample
!c
                    call get_mul_corr(is3,runs(1,is),npers(is),forcecalc)
                    write(6,*) 'gudrun_dcs> Got sample multiple scattering correction'
                end do
!c
!c get Placzek correction for sample
!c
                call get_Placzek(4,runs(1,1),npers(1),forcecalc)
                write(6,*) 'gudrun_dcs> Got sample Placzek correction'
!c
!c Read in the abs correction data
!c
                call get_abs_corr(runs(1,1),npers(1),forcecalc)
                write(6,*) 'gudrun_dcs> Got absorption corr'
!c
!c Now correct the sample and containers selecting each spectrum in turn, ignoring bad detectors if present. 
!c
                igood=0
                mproc=nspecprocgood
                call reallocate1d_i(specproc,mproc)
                do while (igood.lt.nspecprocgood)
                    nspecproc=0
                    do while(nspecproc.lt.mproc.and.igood.lt.nspecprocgood)
                        igood=igood+1
                        nspecproc=nspecproc+1
                        specproc(nspecproc)=specprocgood(igood)
                    end do
                    nspecf=specproc(1)
                    nspecl=specproc(nspecproc)
                    write(6,*) 'gudrun_dcs> Sample and containers: first spec., last spec., nspecproc ',nspecf,nspecl,nspecproc
!c
!c For each sample and container, subtract background, and normalise to the 
!c vanadium and subtract the multiple scattering
!c
                    do is=1,nsample
                        if(is.eq.1) then
                            write(6,*) 'gudrun_dcs> For sample'
                        else
                            write(6,*) 'gudrun_dcs> For container ',is-1
                        endif
!c
!c shifted counter
!c
                        is3=is+3
!c
!c Form sums of sample counts
!c
                        call get_sum(is3,nspecwrt)
                        write(6,*) 'gudrun_dcs> Got sample or container counts'
!c
!c Divide sample data by smoothed monitor
!c
                        call divide_by_mon(is3,nspecwrt)
                        write(6,*) 'gudrun_dcs> Normalised to monitor'
!c
!c Subtract sample background from sample
!c
                        call subtract_bak(is3,wavemin,wavemax,nspecwrt)
                        write(6,*) 'gudrun_dcs> Subtracted sample background'
!c
!c Normalise sample data by dividing by Vanadium
!c
                        call divide_by_van(is3,nspecwrt,normalisationfactor)
                        write(6,*) 'gudrun_dcs> Divided by vanadium'
!c
!c Do the multiple scattering correction for the sample
!c
                        call do_mul_corr(is3,nspecwrt)
                        write(6,*) 'gudrun_dcs> Done sample multiple scattering correction'
                    end do
!c
!c Do the absorption correction
!c
                    call do_abs_corr(nspecwrt)
                    write(6,*) 'gudrun_dcs> Done absorption correction'
!c
!c Do Placzek correction if needed (STEMP>0)
!c
                    if(stemp(1).gt.0.0) then
                        call do_Placzek_corr(nspecwrt)
                        write(6,*) 'gudrun_dcs> Done Placzek correction'
                    endif
!c
!c Subtract self term if present
!c
                    if(ifoundselfe.and.outputunitstype.ne.4) then
                        call do_self_corr(outputunitstype,nspecwrt)
                        write(6,*) 'gudrun_dcs> Done self subtraction'
                    endif
!c
!c Merge normalised data onto a common Q-scale
!c
                    call merge_det(wavemin,wavemax,ngrpwrt,nspecwrt)
                    write(6,*) 'gudrun_dcs> Done merge'
!c
!c go back and do the next batch of spectra if specproc(nspecproc).lt.nspec
!c
                end do
!C
!C Finish the merge
!C
                call finish_merge(acceptance,mergepwr,nbacksub,ngrpwrt,qshell,rshell,qwindow,rmax)
                write(6,*) 'gudrun_dcs> Finished merging data for sample ',runs(1,1)(1:len_trim(runs(1,1)))

!c Write out the spectrum ratios

                fname=runv(1)
                call change_ext(fname,'vanrat')
                call reallocate1d_r(wavebound,nspec)
                do i=1,nspec
                    wavebound(i)=i
                end do
                call w_diag_file(fname,nspec,wavebound,vancor,vancore)
                do is=1,nsample
                    fname=runs(1,is)
                    call change_ext(fname,'samrat')
                    do i=1,nspec
                        wavebound(i)=i
                        vancor(i)=samcor(i,is)
                        vancore(i)=samcore(i,is)
                    end do
                    call w_diag_file(fname,nspec,wavebound,vancor,vancore)
                end do

            end if
        end if
            
!c Get the next line

            call getaline(15)
            ifoundline=nwords.gt.0
            if(ifoundline) ifoundline=line(ncf(1):ncl(1)).ne.'END'

    end do

    call endProgram(smoovlimit,amushield)
      
    end 
    
    subroutine endProgram(smoovlimit,amushield)
        
        use ran1_routines
        use normalisation_routines
        use amutshield_routines


        real smoovlimit,amushield
      
        close(15)
        call saveran1(10)

!c Save the vanadium smooth parameters for next time

        call save_nsmoov(smoovlimit)

!c Save the shielding attenuation factors if non-zero

        call save_amutshield(amushield)

        call exit()

    end subroutine endProgram
     
