      subroutine getBeamParameters()

      include 'dimension.inc'
      include 'beam.inc'
      include 'sam_par.inc'
      include 'inputfilestrings.inc'
      logical found1,found2,found3,found4

      theta2theta=.false.
      beamcompensation=.false.
      detectorcompensation=.false.

c      character*256 line
c
c Read and generate beam parameters for this experiment
c
      open(15,file='BeamParameters.txt',status='old',iostat=ierr)
      if (ierr.eq.0) then


c Ignore empty lines or lines that begin with spaces. The smallest allowable
c character is a blank - ASCII decimal 32. Anything less than or equal to this should be ignored.

         line=' '
         m=0
         do while (test1(1).le.32)
            read(15,'(a)') line
            m=m+1
         end do

c number of characters in line

         nchar=Len(line)
c         write(6,*) nchar,line

c number of c

c The first valid line read must contain the strings which signify the spacing between
c values (spc2), the spacing between the last value and the following comment (spc5), and the
c pathseparator character. These values are delineated by apostrophes as the number of
c spaces is important for these values. Hence a special reading method is required.
c
c Find the first and second occurrence of a '
c
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
         write(6,*) spc2(1:lenspc2),';',spc5(1:lenspc5),';'
     1,pathseparator(1:lenpathseparator)


         call getaline(15,ierr)
         if(line(ncf(1):ncl(nwords)).eq.'FLATPLATE') then
            vgeom=2
            sgeom=2
         else
            vgeom=1
            sgeom=1
         endif
         call init_beam(15)

c Continue to read lines, looking for one containing the text "Theta-2Theta scanning?"

         found1=.false.
         found2=.false.
         found3=.false.
         found4=.false.
         do while (ierr.eq.0.and.(.not.found1.or.
     *.not.found2.or..not.found3.or..not.found4))

            call getaline(15,ierr)
            if(nwords.gt.0) then

                itest=index(line, "Theta-theta scanning?")
                if (itest.gt.0) then

                   found1=.true.
                   read(line,*) itheta2theta
                   theta2theta=itheta2theta.eq.1

                endif
                itest=index(line, "Fixed slits?")
                if (itest.gt.0) then

                   found2=.true.
                   read(line,*) ibeamcompensation
                   beamcompensation=ibeamcompensation.eq.0

                endif

            endif
         end do
      end if
      write(6,*) theta2theta,beamcompensation
      close(15)
      return
      end

      subroutine getaline(unit,ierr)

      include 'inputfilestrings.inc'
      integer*4 unit

c Procedure to read and parse character strings, using the spc2 string for between
c words and spc5 to indicate end of line.

c Find the next valid line

      read(unit,'(a)',iostat=ierr) line
      nwords=0
      do while (nwords.eq.0.and.ierr.eq.0)

c Parse line into words

         call parse(line,nwords,ncf,ncl,nchartext,nwordsdim
     1,spc2,spc5,lenspc2,lenspc5,ierrline)
         if(nwords.eq.0) read(unit,'(a)',iostat=ierr) line
      end do

      return
      end

      subroutine parse(line,nwords,ncf,ncl,nchartext,nwordsdim
     1,spc2,spc5,lenspc2,lenspc5,ierrline)

c parses the text in sometext to separate words with the assumed delimiter spc2
c and spc5 is the string between the last value and any comments in the line.
c ncf and ncl contain the first and last character numbers of each word
c ncharword is the maximum number of characters in a word, nworddim is dimension of the
c word array, nspcword is the number of spaces between words needed to decide that the
c end of text has been reached.

      character*256 line      !text to be parsed
      character*20 spc2,spc5      !Spacing and end of line character strings.
      integer*4 lenspc2,lenspc5      !length of spacing and end-of-line character strings
      integer*4 nwords,ncf(nwordsdim),ncl(nwordsdim),nchartext,ierrline
      integer*4 i

c step through the characters and identify words (separated by spc2 strings). Extra spaces
c are ignored unless the number of spaces exceeds nspcword, in which case the search
c for words terminates

      nwords=0
      ierrline=0

c Get number of characters up to end of last valid variable

      nchartext=len_trim(line)
      if(nchartext.lt.2) then
c If there is nothing on the line at all will report an error
         ierrline=1
         return
      endif

c Search for words in the current string

      i=1
      do while(nwords.lt.nwordsdim.and.i.lt.nchartext)

c Ignore empty space before a word

            do while(line(i:i).eq.' '.and.i.lt.nchartext)
            i=i+1
         end do

         if(i.lt.nchartext) then

            index1=index(line(i:nchartext),spc2(1:lenspc2))-1
            if(index1.gt.0) then
               nwords=nwords+1
               ncf(nwords)=i
               ncl(nwords)=i+index1-1
               i=ncl(nwords)+1
            else
c If index1.lt.1 it can only mean this is the last word of the set
               nwords=nwords+1
               ncf(nwords)=i
               ncl(nwords)=nchartext-1
               i=ncl(nwords)+1
            endif
         endif
      end do

      return
      end

***********************************************************************************
*
*      init_beam.FOR
*
*      A K Soper, March 2001
*
*      Reads the beam profile and other parameters from unit nin
*       to be used for the calculation of attenuation and multiple
*      scattering corrections.
*
***********************************************************************************
      subroutine init_beam(nin)
      implicit none
      include 'dimension.inc'
      include 'beam.inc'
      include 'inputfilestrings.inc'
      integer nin            !unit number to read data from
      integer i            !do loop counter
      real am            !maximum value of beam profile
      real lowang,highang      !lowest and highest scattering angles for corrs.
      real range            !range of scattering angles to be used.
      real pi            !pi
C
C READ NO. OF PROFILE VALUES AND VALUES
C
      READ(NIN,*) NPROF
      if(nprof.gt.mprof) nprof=mprof
      READ(NIN,*) (PROFIL(I),I=1,NPROF)
C
C NORMALIZE PROFILE VALUES TO MAXIMUM VALUE
C
      AM=0.
            DO 10 I=1,NPROF
              AM=AMAX1(AM,PROFIL(I))
10          CONTINUE
            DO 11 I=1,NPROF
              PROFIL(I)=PROFIL(I)/AM
11          CONTINUE
c
C INPUT INTEGRATION PARAMETERS
C
      READ(NIN,*) stepa,stepm,nslice
C
C INPUT NO. OF degrees between corrections
C
      READ(NIN,*) NDEG
C
C READ position of edges of incident beam, and scattered beam relative to
c centre of sample.
C
      READ(NIN,*) B,A,HBDOWN,HBUP
      read(nin,*) B1,A1,HSBDOWN,HSBUP
c
c define profile step
c
      PRSTEP=(A-B)/(NPROF-1)
      return
      end

**********************************************************************************
*
*      set_corr_ang.FOR
*
*      A K Soper, February 2003
*
*      Generates theta and phi values for specified sample geometry.
*
*      For cylinders the z-axis runs along the axis of the cylinder and is assumed
*       to be at right angles to the incident beam.
*
*      For flat plates the z-axis is assumed to be perpendicular to the surface of
*      slab, with the positive axis on the side of the transmitted beam
*
***********************************************************************************
      subroutine set_corr_ang(ngeom,angfirst,anglast,azidet
     *,sang,ang,azi,angrot,nang,nazi,ttheta,nangles)
      IMPLICIT NONE
      include 'dimension.inc'
      include 'beam.inc'
c
c internal variables
c
      integer ngeom            !1 = cylinder geometry, 2 = flat plate
      integer nang,nazi,nangles            !no. of theta and phi values for corrections
      integer i,ib,j,index,il,k      !internal indices
      integer is,id            !internal indices
      real sang             !Sample rotation angle (flat plate only)
      real ang(mcorrang)      !theta values for corrections
      real azi(mcorrang)      !azimuthal angle values for corrections
      real angrot(mcorrang)  !sample rotation angles
      real angmax,angmin      !maximum and minimum theta values
      real azimax,azimin      !maximum and minimum azi values
      real angfind,azifind      !theta and phi values corresponding to specified detector
      real angst,azist      !angular step sizes
      real tthetat,phit,ttheta(mcorrang)
      real angfirst,anglast,angstep,azidet
c
c step through good detectors and determine the maximum and minimum values
c of theta and phi for the respective geometries
c
      angmax=0.0
      angmin=180.0
      azimax=0.0
      azimin=180.0
      tthetat=angfirst
c      angstep=anglast-angfirst
      phit=azidet
      nangles=0
      do while(tthetat.le.anglast*1.0001.and.nangles.lt.mcorrang)

c get the corresponding scattering and azimuthal angles for this detector

         nangles=nangles+1
         if(theta2theta) then
            angrot(nangles)=0.5*tthetat+sang
         else
            angrot(nangles)=sang
         endif
         call get_corr_ang(ngeom,angrot(nangles),tthetat,phit
     1,angfind,azifind)

c         write(6,*) tthetat,phit,angfind,azifind
         ang(nangles)=angfind
         azi(nangles)=azifind
         ttheta(nangles)=tthetat
         tthetat=tthetat+real(ndeg)

      end do
c For cylindrical geometry the scattering angles will appear in azi,
c For flat plate geometry the scattering angles will appear in ang

      if(ngeom.eq.1) then
         nang=1
         nazi=nangles
      else
         nang=nangles
         nazi=1
      endif

      tthetat=angfirst
c      write(6,100) nang,nazi,(tthetat+(i-1)*real(ndeg)
c     *,angrot(i),ang(i),azi(i),i=1,nangles)
c100   format(2(1x,i5)/(4(1x,f10.5)))

c      stop
      return
      end
***********************************************************************************
*
*      get_corr_ang.FOR
*
*      A K Soper, February 2003
*
*      Generates theta and phi values for specified sample geometry.
*
*      For cylinders the z-axis runs along the axis of the cylinder and is assumed
*       to be at right angles to the incident beam.
*
*      For flat plates the z-axis is assumed to be perpendicular to the surface of
*      slab, with the positive axis on the side of the transmitted beam
*
***********************************************************************************
      subroutine get_corr_ang(ngeom,sangd,tdetd,pdetd,tfind,pfind)
c
c A K Soper, February 2003
c
c converts specified detector coordinates (z-axis in downstream beam direction), x-axis
c in horizontal plane left of z-axis, y-axis vertically upwards) to coordinates appropriate
c to cylinders (z axis vertically upwards, x axis along downstream beam direction, y-axis
c in horizontal plane to left of beam direction when facing in same direction as beam) or
c flat plate (z-axis perpendicular to plane of slab on downstream side, y-axis perpendicular
c this is vertical plane, x-axis perpendicular to z-axis in horizontal plane.
c
      implicit none
      integer ngeom            !Sample geometry
      real sangd            !sample angle (deg)
      real sangr            !sample angle (rads)
      real tdetd,pdetd      !Detector coordinates in instrument coordinate axes (deg)
       real tdetr,pdetr      !Detector coordinates in instrument coordinate axes (rad)
      real tfind,pfind      !Detector coordinates in sample coordinates.
      real piconv            !converts angles to radians
      real costs,xs,ys      !temporary values
      real x,y,z,xprime,yprime,csangr,ssangr
c      write(6,*) ngeom,sangd,tdetd,pdetd
      piconv=4*atan(1.0)/180.0
c
c convert to radians
c
      sangr=sangd*piconv
      tdetr=tdetd*piconv
      pdetr=pdetd*piconv
c
c Determine sample geometry
c
      if(ngeom.eq.1) then
c
c cylindrical geometry - corrections will be symmetric about tfind=90
c
         costs=abs(sin(tdetr)*sin(pdetr))
         if(abs(costs).gt.1.0) costs=sign(1.0,costs)
         tfind=acos(costs)/piconv
c
c only need absolute value of pfind, since corrections will be symmetrical
c about pfind=0
c
         xs=sin(tdetr)*cos(pdetr)
         ys=cos(tdetr)
         pfind=abs(atan2(xs,ys)/piconv)

      else if(ngeom.eq.2) then
c
c flat plate geometry, sang specifies the rotation about the vertical detector y axis.
c Only the angle between the detector and the normal to the sample face is needed
c so there will be no azimuthal dependence in all cases. This angle will however depend
c on the detector azimuthal angle if sang is not 0 or 180
c
         z=cos(tdetr)
         x=sin(tdetr)*cos(pdetr)
         y=sin(tdetr)*sin(pdetr)
         csangr=cos(sangr)
         ssangr=sin(sangr)
         costs=z*csangr+x*ssangr
         if(abs(costs).gt.1.0) costs=sign(1.0,costs)
         tfind=acos(costs)/piconv
c
c Get the phi value in the rotated coordinate geometry. Remember slab has been rotated about the instrument
c y-axis (normally the vertical axis).

         xprime=x*csangr+z*ssangr
         yprime=y
         pfind=atan2(yprime,xprime)/piconv

c Adjust any values .le. -180

         do while (pfind.le.-180.0)

            pfind=pfind+360.0

         end do

      else
            write(6,*) 'Specified geometry ',ngeom,' is not possible'
            stop
      endif
      return
      end

C
C      *********************
C      *                   *
C      *  FUNCTION PROBE   *
C      *                   *
C      *********************
C
      FUNCTION PROBE(X,PROFIL,NPROF,PRSTEP,A1,B1)
      DIMENSION PROFIL(NPROF)
      IF(X.GT.A1.OR.X.LT.B1) then
            probe = 0.0
      else
            DIFF=A1-X
            APOS=DIFF/PRSTEP
            NPOS=INT(APOS)
            DIFF=APOS-FLOAT(NPOS)
            NPOS1=NPOS+1
            NPOS2=NPOS+2
            IF(DIFF.ne.0.) then
               PROBE=PROFIL(NPOS1)+(PROFIL(NPOS2)-PROFIL(NPOS1))*DIFF
            else
               PROBE=PROFIL(NPOS1)
            endif
      endif
      RETURN
      END

      function diskfootprint(rad,beamheight,beamwidth,sectheta)

c Calculates the area of intersection of disk of radius rad with rectangular beam of height beamheight and width beamwidth
c incident at an angle sectheta to the normal to the disk.

      real rad,beamheight,beamwidth,sectheta
      real height,width,t,t1,t2,alpha,theta,fullarea1,fullarea2
     1,fullarea
      real diskfootprint

      height=min(rad,beamheight)
      width=abs(beamwidth)
      t=width*abs(sectheta)

c Intersection of top of beam with disk

      t1=sqrt(rad*rad-height*height)

c Maximum value of t

      t2=rad

c Get the area of the intersection in the event that sample is fully exposed to beam

      alpha=acos(t1/rad)
      fullarea1=0.5*height*t1
      fullarea2=0.5*rad*rad*alpha
      fullarea=fullarea1+fullarea2

      if(t.ge.t2) then

         diskfootprint=1.0
         return

      else if (t.le.t1) then

         diskfootprint=t*height/fullarea
         return

      else

         theta=acos(t/rad)
         wt=sqrt(rad*rad-t*t)
         diskfootprint=(fullarea1+0.5*t*wt+0.5*rad*rad*(alpha-theta))
     1/fullarea
      write(6,*) sectheta,t,wt,alpha,theta,height,t1,rad,fullarea1
     1,fullarea
         return

      endif

      end



