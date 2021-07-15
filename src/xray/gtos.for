c Program to remove the single atom scattering using a 3D top hat function
c Copyright Alan Soper April 2005

      program gtos_test

      implicit none

      include 'dimension.inc'
      include 'run_par.inc'
      include 'calibration.inc'
      include 'groups.inc'
      include 'merge_arrays.inc'

      character*80 fname,fnameout,formatout
      character*80 baddetfname     !bad detector file full filename
      character*80 groupsfilefname !groups file full filename
      character*20 uname            !username
      character*20 sttime            !start time and date

      real qdata(mq),dcsdata(mq),dcserr(mq),qmax
      real qscale(mq),deltaq,qstep,dcsbin(mq),errbin(mq)
      real qshell,smoothdat(mq),dcssave(mq)
      real pi,rho,rstep,rtemp(mq),gr1temp(mq)
      real gr2temp(mq),grlimit,sqlimit,rshell,leveltest,minfac
      real sq1temp(mq),sq2temp(mq),sqstart(mq),grstart(mq)
      real lentemp(mgroup),tthtemp(mgroup),phitemp(mgroup)
      real factor,rmax,lorchwindow,rpower
      real sumhighq,rstep1,qnext,qfirst,qlast
      real q,q1,q2,errdat(mq),factor2,r1,r2,rfac
      real q13,q23,sumq3,sumiq3,term,dwamp,dwarg,exparg
      real, dimension(:), allocatable     :: atdensities
      real, dimension(:,:), allocatable   :: rvaldata,gofrdata
      real, dimension(:,:), allocatable   :: qvaldata,sofqdata

      character*2000 sometext,text
      integer*4 ierr,ndata,iq,dataformat,nbroad,nr,nfiles
      integer*4 xcode,ycode            !used by Genie load command
      integer*4 nsumhighq,narg,iqsave,idata,iarg,i,iter,niter

      integer, parameter :: nwordsdim=100
      character*20 spc2,spc5      !Spacing and end of line character strings.
      integer*4 lenspc2,lenspc5      !length of spacing and end-of-line character strings
      integer*4 nwords,ncf(nwordsdim),ncl(nwordsdim),nchartext,ierrline
      integer*4 idataset,ndatasets,iline,nlines,iword

c Set up some values not used in this program

      info=' '
      uname=' '
      sttime=' '
      baddetfname=' '
      groupsfilefname=' '
      xcode=-1
      ycode=6
      lenin=10
      lentemp(1)=0
      tthtemp(1)=0.0
      phitemp(1)=0.0
      ngroupt=1

c Values for parsing strings

      spc2=' '
      spc5='     '
      lenspc2=1
      lenspc5=5

      write(6,100)
100      format(/'gtos - Fourier transform g(r) to S(Q)'
     */' October 2018'/)

c Get input from any supplied arguments

            narg=iargc()
      write(6,*) narg

c First argument is the filename to be processed (no spaces in filename are allowed)

           iarg=0
      ierr=0
      if(iarg.lt.narg.and.ierr.eq.0) then
         iarg=iarg+1
         call getarg(iarg,text)
c      write(6,*) text
         read(text,'(a)',iostat=ierr) fname
      else
         write(6,129)
      endif
129      format(/'gtos> Not enough arguments supplied.')

c Next is whether data are histogram (1) or point (2) data.

      if(iarg.lt.narg.and.ierr.eq.0) then
         iarg=iarg+1
         call getarg(iarg,text)
         read(text,*,iostat=ierr) dataformat
      else
         write(6,129)
      endif

c Next is a factor for the data.

      if(iarg.lt.narg.and.ierr.eq.0) then
         iarg=iarg+1
         call getarg(iarg,text)
         read(text,*,iostat=ierr) factor
      else
         write(6,129)
      endif

c Next is Qstep

      if(iarg.lt.narg.and.ierr.eq.0) then
         iarg=iarg+1
         call getarg(iarg,text)
         read(text,*,iostat=ierr) rstep
      else
         write(6,129)
      endif

c Next is Qmax

      if(iarg.lt.narg.and.ierr.eq.0) then
         iarg=iarg+1
         call getarg(iarg,text)
         read(text,*,iostat=ierr) rmax
      else
         write(6,129)
      endif

c Read the data to be processed

      write(6,101) fname
101      format(/'gtos> Name of the file to be processed: ',a50)

c Open the file and read it

      open(10,file=fname,status='old')

c Ignore lines with # in them

      sometext='#'
      do while(index(sometext,'#').gt.0)
         text=sometext
         read(10,'(a)') sometext
      end do
      write(6,'(a)') text(1:50)

c Determine the number of files to be transformed, and get their densities

      call parse(text,nwords,ncf,ncl,nchartext,nwordsdim
     1,spc2,spc5,lenspc2,lenspc5,ierrline)

c Numbor of files is 1 less than the number of words since the first word is a #

      ndatasets=nwords-1
      write(6,'("gtos> ",i5," sets of data found in file ",a)') 
     *ndatasets
     *,fname
      if(allocated(atdensities)) deallocate(atdensities
     1,rvaldata,gofrdata,qvaldata,sofqdata)
      allocate(atdensities(ndatasets))
! Retrieve the atomic number densities
      iword=1
      do while(iword.lt.nwords)
         iword=iword+1
         read(text(ncf(iword):ncl(iword)),*) atdensities(iword-1)
      end do

c Determine the number of lines with valid data

      iline=0
      ierr=0
      do while (ierr.eq.0)
         call parse(sometext,nwords,ncf,ncl,nchartext,nwordsdim
     1,spc2,spc5,lenspc2,lenspc5,ierrline)
         if(nwords.gt.0) then
            iline=iline+1
         end if
         read(10,'(a)',iostat=ierr) sometext
      end do
      nlines=iline
      write(6,'("gtos> ",i5," lines of data found in file ",a)') nlines
     *,fname
      allocate(rvaldata(nlines,ndatasets),gofrdata(nlines,ndatasets))

c Now read the data. It is assumed that the last line containing data is terminated
c with a new line character. Otherwise the last line will be lost.

      ierr=0
      ndata=0
      qmax=0
      rewind(10)

c Ignore lines with # in them

      sometext='#'
      do while(index(sometext,'#').gt.0)
         read(10,'(a)',iostat=ierr) sometext
      end do
      iline=0
      ierr=0
      do while (ierr.eq.0)
         call parse(sometext,nwords,ncf,ncl,nchartext,nwordsdim
     1,spc2,spc5,lenspc2,lenspc5,ierrline)
         if(nwords.gt.0) then
            iline=iline+1
            read(sometext,*,iostat=ierr) (rvaldata(iline,idataset)
     1,gofrdata(iline,idataset),idataset=1,ndatasets)
         end if
         read(10,'(a)',iostat=ierr) sometext
      end do
      write(6,'("gtos> ",i5," lines of data read from file ",a)') iline
     *,fname
      close(10)

c Set up the Q-scale on which to bin this data

      if (nlines.gt.1) then
         deltaq=rvaldata(2,1)-rvaldata(1,1)
      else
         write (6,'(a)') 'gtos> Insufficient data to Fourier transform!'
         stop
      end if
      write(6,104) deltaq
104      format(/'gtos> Specify the width of the r bins: ',f10.5)

c Determine the total number of bins required and set up the Q-scale
c This will be histogram format, starting from zero

      qmax=rvaldata(nlines,1)
      nq=nint(qmax/deltaq)+1
      if(nq.gt.mq) nq=mq
      write(6,102) qmax
102   format(/'gtos> Maximum r found is: ',f10.5)
      allocate(qvaldata(nq,ndatasets),sofqdata(nq,ndatasets))
      
c Supply a factor to multiply these data by

      write(6,109)
109      format(/'gtos> Specify a factor for these data: '
     *,$)
      factor2=factor*factor

c Step through the datasets and Fourier transform each one.

      do idataset=1,ndatasets
      rho=atdensities(idataset)
      ndata=0
      do while(ndata.lt.nlines.and.ndata.lt.mq)

         ndata=ndata+1
         qdata(ndata)=rvaldata(ndata,idataset)
         dcsdata(ndata)=gofrdata(ndata,idataset)
         qmax=max(qmax,qdata(ndata))

      end do

c Save the first and last input Q values

      qfirst=qdata(1)
      qlast=qmax

      idata=0
      do while(idata.lt.ndata)

         idata=idata+1
         dcsdata(idata)=dcsdata(idata)*factor
         dcserr(idata)=0.0
 
      end do

c If data format is point, we need to convert to histogram

      if(dataformat.eq.2) then


c Each bin boundary will be placed midway between the existing bins

         q1=qdata(1)
         do iq=2,ndata
            q2=qdata(iq)
            qnext=0.5*(q1+q2)
            qdata(iq)=qnext
            q1=q2
         end do
c At this point qdata(1) has not been modified
         qdata(1)=2*qdata(1)-qdata(2)
c Will need one more boundary (containing zero data)
         if(ndata.lt.mq) then
            ndata=ndata+1
            qdata(ndata)=2*q2-qdata(ndata-1)
            dcsdata(ndata)=0.0
            dcserr(ndata)=0.0
         endif
      endif

      fnameout=fname(1:index(fname,'.')-1)//'.hist'
      call w_diag_file(fnameout,ndata,qdata,dcsdata,dcserr)

      qscale(1)=0.0
      qstep=deltaq*0.5
      do iq=2,nq
         qscale(iq)=qscale(iq-1)+qstep
         qstep=deltaq
      end do

c Rebin the data onto this scale and write it out again

      call rebinq(ndata,qdata,dcsdata,nq,qscale
     1,dcsbin,mq,mq,2,1)
      call rebinq(ndata,qdata,dcserr,nq,qscale
     1,errbin,mq,mq,3,0)
      !Infill points below Q=qdata(1)
      q1=qscale(1)
      i=1
! Fill in data at low Q
      do while (qscale(i).lt.qdata(1).and.i.lt.nq)
         dcsbin(i)=dcsdata(1)
         errbin(i)=dcserr(1)
         i=i+1
      end do
      fnameout=fname(1:index(fname,'.')-1)//'.qbin'
      call w_diag_file(fnameout,nq,qscale,dcsbin,errbin)

c If qshell is -ve then need to sum data for Q > -qshell to determine high Q limit

         !Find the constant best suited to give g(0)=-grlimit
         sumhighq=1.0

c Form the average high Q limit and subtract it.

         do iq=1,nq-1
            smoothdat(iq)=sumhighq
            dcsdata(iq)=smoothdat(iq)
            if(dcsbin(iq).ne.0.0) then
               dcssave(iq)=dcsbin(iq)-dcsdata(iq)
            else
               dcssave(iq)=0.0
            endif
         end do
         nbroad=0
         qshell=0.0

      fnameout=fname(1:index(fname,'.')-1)//'.qsmooth'
      call w_diag_file(fnameout,nq,qscale,dcsdata,errbin)

      fnameout=fname(1:index(fname,'.')-1)//'.qsub'
      call w_diag_file(fnameout,nq,qscale,dcssave,errbin)

c set the r step value to calculate Fourier Transforms

      pi=3.141592653
      rstep1=pi/qscale(nq)

c set up a radius scale for F.T.

      do iq=1,nq
            rtemp(iq)=rstep*(iq-1)
      end do

c Fourier transform to Q-space

      call stogtos(1,rho,nq,qscale,dcssave,nq,rtemp,gr1temp)

c Write diagnostic file

      fnameout=fname(1:index(fname,'.')-1)//'.sq1'
      call w_diag_file(fnameout,nq,rtemp,gr1temp,gr1temp)
      do iline=1,nq
         qvaldata(iline,idataset)=rtemp(iline)
         sofqdata(iline,idataset)=gr1temp(iline)
      end do

c Do next dataset
  
      end do

      fnameout=fname(1:index(fname,'.')-1)//'.sofq'

c Open the file and write it

      open(10,file=fnameout,status='unknown')

      write(formatout,332) '(a1,',ndatasets,'(1x,f10.6))'
332   format(a,i3.3,a)
      write(10,formatout) '#',(atdensities(idataset)
     1,idataset=1,ndatasets)
      write(formatout,332) '(',ndatasets,'(1x,f8.3,1x,f10.5))'
      do iline=1,nq
         write(10,formatout) (qvaldata(iline,idataset)
     1,sofqdata(iline,idataset),idataset=1,ndatasets)
      end do

      stop
      end

      subroutine tophat3d(nsmoo,nchan,detcount,errcount)
c
c does a square-wave smoothing with width proportional to total number of points in array
c input in detcount, smoothed result in errcount
c
      real*4 detcount(*),errcount(*)
c
c determine the first and last non-zero bins of input array
c
      nfirst=0
      nlast=nchan
      do ic=1,nchan
         errcount(ic)=0.0
         if(nfirst.eq.0.and.detcount(ic).ne.0.0)
     1 nfirst=ic
         if(detcount(ic).ne.0.0) nlast=ic
      end do
c
c smooth the input array with a top hat function in 3D
c
      n=nsmoo
      rn=n
      rn2=2.0*rn+1
      factv=1.0/(rn2*rn2*rn2)
      do m1=nfirst,nlast
         errcount(m1)=0.0
         m=m1-1
         rm=m
         rm2=2.0*rm   
         rmn1=rm2*rm2-rn2*rn2+1.0
         factt=factv/rm2
         mmn=m-n
         mpn=m+n
c
c if mmn .lt. 1, then need to include the cases when the whole volume is 
c included
c
         if(mmn.lt.1) then
            mmax=iabs(mmn)
c
c volume of smallest sphere
c
            l=1
            if(l.lt.nfirst) then 
               add=detcount(nfirst)
            else
               add=detcount(l)
            endif
            errcount(m1)=errcount(m1)+add*factv
               if(mmax.gt.0) then
               do l=1,mmax
                  lref=l+1
                  rl=l
                  rl2=12.0*rl*rl+1.0
                  fact=2.0*rl2*factv
                  if(lref.le.nfirst) then
                     add=detcount(nfirst)
                  else if(lref.ge.nlast) then
                     add=detcount(nlast)
                  else
                     add=detcount(lref)
                  end if
                  errcount(m1)=errcount(m1)+add*fact
               end do
               end if
               mmin=mmax+1
         else
               mmin=mmn            
         endif
         mmax=mpn
         do l=mmin,mmax
            lref=l+1
            rl=l
            rl2=rl*rl
            fact=rm2*(12.0*rl2+1.0)-3.0*rl*(4.0*rl2+rmn1)
            fact=fact*factt
            if(lref.lt.nfirst) then
               add=detcount(nfirst)
            else if(lref.gt.nlast) then
               add=detcount(nlast)
            else
               add=detcount(lref)
            end if
            errcount(m1)=errcount(m1)+fact*add
         end do
      end do
      return
      end

      subroutine tophat3dmod(nchan,delta,rmin,cons,xin
     1,yin,yout,expamp,expdecay,expstretch)
c
c does a square-wave smoothing with width proportional to total number of points in array
c input in detcount, smoothed result in errcount
c
      real*4 xin(*),yin(*),yout(*)

c Determine the first channel at or after rmin

      imin=0
      ic=0
      do while (ic.lt.nchan.and.imin.eq.0)

         ic=ic+1
         if(xin(ic).ge.rmin) imin=ic

      end do

      do ic=1,nchan
         r=xin(ic)
         yout(ic)=0.0
         if(r.lt.rmin) then
            yout(ic)=yin(ic)+cons
c            yout(ic)=yin(ic)-yin(imin)
         endif
           delr=delta*r
         if(delr.gt.0.0) then
            delr3=3.0/(delr*delr*delr)
            pofr=delr3*(sin(delr)-delr*cos(delr))
            if(r.lt.rmin) then
c             yout(ic)=yout(ic)-cons*pofr
            else
               yout(ic)=-yin(ic)*pofr/(1.0-pofr)
               expon=expdecay*r**expstretch
               if(expon.gt.0.and.expon.lt.30.0) then
                  yout(ic)=yout(ic)+expamp*exp(-expon)
               endif
!               yout(ic)=min(yout(ic),(yin(ic)+cons)) !Ensure the correction does not drive g(r) unphysically negative
            endif
         endif
      end do
      return
      end

      subroutine w_diag_file(fname,nx,x,y,e)
      character*80 fname
      integer*4 nx
      real*4 x(*),y(*),e(*)
      real*4 rat,err
      open(10,file=fname,status='unknown')
      write(10,100) '#',fname
100      format(a1,1x,a)
      write(10,100) '#'
      write(10,100) '#'
      write(10,100) '#'
      do ic=1,nx
            rat=y(ic)
            if(e(ic).gt.0.0) then
               err=sqrt(e(ic))
            else
               err=0.0
            endif
            write(10,101) x(ic),rat,err
      end do
      close (10)
101      format(1x,3(1x,e13.6))
      return
      end

      subroutine w_int(fname,uname,sttime,baddetfname
     2,groupsfilefname,xcode,ycode,lentemp,tthtemp,phitemp)
      parameter (nwrmax=100)
      include 'dimension.inc'
      include 'run_par.inc'
      include 'calibration.inc'
      include 'groups.inc'
      include 'merge_arrays.inc'
      integer*4 xcode,ycode            !used by Genie load command
      integer*4 nskip      !how many lines to skip
      real*4 lentemp(*),tthtemp(*),phitemp(*)
      character*60 fname      !name of file to write data to
      character*80 baddetfname     !bad detector file full filename
      character*80 groupsfilefname !groups file full filename
      character*20 uname            !username
      character*20 sttime            !start time and date
      character*32 formatout
      character*2 blk
c
c After the initial 3 lines nskip is the number of lines of info that need to be skipped
c to get to the actual information
c
        nskip=10
c
c set up number of output files to write to
c
      nt=(ngroupt+(nwrmax-1))/nwrmax
      ilen=index(fname,' ')-1
      iwr1=1
      do it=1,nt
c
c number of elements to write for each x value
c
            iwr=ngroupt-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
            iwr2=iwr1+iwr-1
c
c setup output format specification
c
            nblock=2*iwr+1
            write(formatout,332) '(a1,',nblock,'(1x,e14.7))'
332      format(a4,i3.3,a11)
            write(6,*) formatout
            lrec=nblock*15+1
            if(lrec.lt.82) lrec=82
            write(blk,331) it
331      format(i2.2)
            fname=fname(1:ilen)//blk(1:2)
            write(6,101) nq-nqfirst+1,ngroupt,fname
101            format(1x,'Writing ',i6,' points in ',i3,' groups to ',a)
            open(10,file=fname,status='unknown',form='formatted'
     *,recl=lrec)
                  write(10,99) fname
                  write(10,99) info
                  write(10,100) ngroupt,nq
                  write(10,100) nskip
                  write(10,99) uname
                     write(10,99) sttime
                  write(10,99) baddetfname
                  write(10,99) groupsfilefname
                  write(10,100) xcode
                  write(10,100) ycode
                  write(10,102) lenin
99      format('# ',a)
100      format('#',2(1x,i5))
102      format('#',6(1x,e14.7))
                  xo=0.0
c
c write out group secondary flight paths
c
                  write(10,formatout) '#',xo
     *,(lentemp(ir),xo,ir=iwr1,iwr2)
c103      format(a1,037(1x,e14.7))
c
c write out group scattering angles
c
                  write(10,formatout) '#',xo
     *,(tthtemp(ir),xo,ir=iwr1,iwr2)
c
c write out group azimuthal angles
c
                  write(10,formatout) '#',xo
     *,(phitemp(ir),xo,ir=iwr1,iwr2)
c
c write the data for this range of groups
c
                  do iq=nqfirst,nq
                        write(10,formatout) ' ',qbound(iq)
     *,(aggsweights(iq,ir),aggeweights(iq,ir),ir=iwr1,iwr2)
                  end do
                  iwr1=iwr2+1
            close(10)
      end do
      return
      end

      subroutine rebinq(nold,xold,yold,nnew,xnew,ynew,mold,mnew
     *,iopt,iopt1)
c
c Performs a rebin from nold channels with x-values xold, to nnew channels
c with x-values xnew. Partly overlapping channels are left empty at the 
c beginning and end of the new array.
c
c The supplied x-values are assumed to be BIN BOUNDARIES throughout 
c Thus there nold values of xold, but only nold-1 values of yold. Same with nnew
c xnew and ynew.
c 
c
c In this version, where the new bin boundary occurs between two old bin 
c boundaries, linear interpolation is made between the neighbouring values,
c rather than simple partitioning between the two.
c
c 09/05/2008 Corrected for binning errors at beginning and end of range
c
      integer*4 nold            !number of old values
      integer*4 nnew            !number of new values
      integer*4 iopt            !1 = histograms,2 = ratios, 3 = ratio std. dev.
      integer*4 iopt1            !0=no printout, 1 = printout
      integer*4 oldref,newref      !intermediate counters
      integer*4 oldmin,oldmax !range of old values to be interpolated over
      real*4 xold(mold)      !old x- boundary values 
      real*4 yold(mold)      !old y-values (nold values)
      real*4 xnew(mnew)      !new x- boundary values 
      real*4 ynew(mnew)      !new y-values (nnew values)
      real*4 lbold,ubold       !temporary store for lower and upper boundaries of current bin
      real*4 lold,uold,lnew,unew,llim,ulim      !intermediate limit values
      real*4 binfrac            !accumulates fractional contribution to bin
      real*4 binwidth
      real*4 bincentre
      real*4 oldbinwidth,newbinwidth
      real*4 fact,fact1,grad,add,cons
      if(iopt1.eq.1) then
            open(25,file='junk.dat',status='unknown')
      endif
c
c zero all new elements
c
      do i=1,nnew
            ynew(i)=0.0
      end do

c If the number of x values is less than 3, then routine will not work. 
c Therefore simply find the nearest new bin which is closest to this old value

      if(nold.lt.3) then

            xoldl=xold(1)
            xoldu=xold(2)
c
c centre of current old bin
c

            xoldbin=0.5*(xoldu+xoldl)
            xnewl=xnew(1)
            inewref=1
            xdiffmin=0.0
            do i=2,nnew
                  xnewu=xnew(i)

c Centre of new bin

                  xnewbin=0.5*(xnewu+xnewl)
                  xdiff=abs(xnewbin-xoldbin)
                  if(i.eq.2.or.xdiff.lt.xdiffmin) then
                        xdiffmin=xdiff
                        inewref=i-1
                  endif
                  xnewl=xnewu
            end do

c Only increment the bin if the difference is less than the width of the bin

            if(xdiffmin.le.0.5*(xnew(inewref+1)-xnew(inewref))) then
                  ynew(inewref)=yold(1)
            endif
            return

      endif

c Search the old data for zeros. The first non-zero value determines the
c first value to be interpolated. The last non-zero value determines the 
c last value to be interpolated.

      oldmin=0
      oldmax=0
      do i=1,nold-1
            if(yold(i).ne.0.0.and.oldmin.eq.0) oldmin=i
            if(yold(i).ne.0.0) oldmax=i
      end do

c Save the minimum and maximum non-zero x values in the old data

      oldrangemin=xold(oldmin)
      oldrangemax=xold(oldmax+1)

      if(iopt1.eq.1) write(25,*) oldmin,oldmax,oldrangemin,oldrangemax
      oldref=oldmin
c
c lower boundary of the current old bin
c
      lold=xold(oldref)
c
c upper limit on current old bin - this is set to the centre of the next old bin
c
      uold=0.5*(xold(oldref+1)+xold(oldref+2))
c
c define the boundaries of the first new bin
c
      newref=1
      lnew=xnew(newref)
      unew=xnew(newref+1)
      do while(uold.le.lnew.or.lold.gt.lnew)
c
c Ensure upper limit of current old bin is ABOVE lower limit of current
c new bin
c
         if(uold.le.lnew) then
            oldref=oldref+1
            lold=uold
            uold=0.5*(xold(oldref+1)+xold(oldref+2))
c
c Ensure lower limit of current old bin is BELOW or EQUAL to lower limit of 
c current new bin
c
         else if(lold.gt.lnew) then
            newref=newref+1
            lnew=unew
            unew=xnew(newref+1)
         endif
      end do
      binfrac=0.0
300      continue
      grad=2*(yold(oldref+1)-yold(oldref))
     $/(xold(oldref+2)-xold(oldref))
      cons=yold(oldref)
c
c form the limits to be integrated over
c
      ulim=min(unew,uold)
      llim=max(lnew,lold)
c
c set up bin width and centre
c
      binwidth=ulim-llim
      bincentre=0.5*(ulim+llim)-0.5*(xold(oldref)+xold(oldref+1))
      oldbinwidth=uold-lold
      newbinwidth=unew-lnew
c
c Form the binning factor for the current old bin and add to new bin
c      
      fact=binwidth/newbinwidth
      fact1=oldbinwidth/newbinwidth
      binfrac=binfrac+fact
c
c integrate the straight line between yold(oldref) and yold(oldref+1) 
c between the lower and upper boundaries
c
      add=fact*((grad*bincentre)+cons)
      if(iopt.eq.1) then
            ynew(newref)=ynew(newref)+add
      else if(iopt.eq.2) then
            ynew(newref)=ynew(newref)+add
          else
            ynew(newref)=ynew(newref)+add*fact1
          endif
      if(iopt1.eq.1) then
            write(25,*) lold,uold,lnew,unew,llim,ulim
            write(25,*) binwidth,bincentre,cons,grad
            write(25,*) add,fact,binfrac,ynew(newref)
            write(25,*) oldref,oldmax,newref,nnew
      endif
c
c if ulim.eq.unew, time to move on to the next new bin
c
      if(ulim.ge.unew) then
c
c for iopt=2 or 3 divide by bin fraction unless the bin is less than 90% filled
c
            if(iopt.gt.1) then
                  if(binfrac.gt.0.9) then
                        if(iopt.eq.2) then
                              ynew(newref)=ynew(newref)/binfrac
                        else if(iopt.eq.3) then
                              ynew(newref)=ynew(newref)
     */(binfrac*binfrac)
                        endif
                  else
                        ynew(newref)=0.0
                  endif
            endif
            newref=newref+1
c
c exit when all new bins are exhausted
c
            if(newref.ge.nnew) then
                  if(iopt1.eq.1) close(25)
                  return
            endif
            binfrac=0.0
            lnew=unew
            unew=xnew(newref+1)
            go to 300
      else
c
c exit when all old bins are exhausted
c
            if(oldref.ge.oldmax-1) then

c Ensure interpolation proceeds all the way to the last bin

               if(unew.lt.oldrangemax) then

                  uold=oldrangemax

               else
c
c zero the last new bin, in case it is only part filled
c
                  ynew(newref)=0.0
                  if(iopt1.eq.1) close(25)
                  return

               endif

            else
c
c else move on to the next old bin
c
               oldref=oldref+1
               lold=uold
               uold=0.5*(xold(oldref+1)+xold(oldref+2))

            endif
            go to 300
      endif
      end

      subroutine stogtos(ntype,dens,lptin,xin,yin,lptout,xout,yout)
      real*4 xin(*),yin(*),xout(*),yout(*)
      real*8 j1qr,q1,q2,q,r1,r2,r,sum
      real*8 j1q1r,j1q2r
      pi=3.141592653
c
c ntype=1 means gtos, otherwise stog
c
      if(ntype.eq.1) then
         trfac=4.0*pi*dens
      else
         trfac=0.5/pi/pi/dens
      endif
      yin(lptin)=0.0
      r1=xout(1)
      do i=1,lptout-1
         r2=xout(i+1)
         r=0.5*(r1+r2)
         rfac=trfac/(r*r*r)
         sum=0.
         q1=xin(1)
         j1q1r=j1qr(q1,r)
         do j=1,lptin-1
            q2=xin(j+1)
            qdel=q2-q1
            qr=0.5*r*(q1+q2)
            j1q2r=j1qr(q2,r)
c Form integral of (Q sin Qr)/r from Q1 to Q2
            sum=sum+yin(j)*(j1q2r-j1q1r)
            q1=q2
            j1q1r=j1q2r
         end do
         yout(i)=rfac*sum
         r1=r2
      end do
      yout(lptout)=0.0
      return
      end

      subroutine stogtoslorch(ntype,dens,lptin,xin,yin
     1,lptout,xout,yout,qwindow,broad)

c Integrate Q sin Qr over range for each bin

      real*4 broad   !Determines how the broadening changes with r
      real*4 xin(*),yin(*),xout(*),yout(*)
      real*8 j1qr,q1,q2,q,r1,r2,r,sum
      real*8 j1q1r,j1q2r
      real*4 qwindow,bwindow,bwidth

      pi=4.0*atan(1.0)
c
c ntype=1 means gtos, otherwise stog
c
      if(ntype.eq.1) then
         trfac=4.0*pi*dens
      else
         trfac=0.5/pi/pi/dens
      endif
      r1=xout(1)
      do i=1,lptout-1
         r2=xout(i+1)
         r=0.5*(r1+r2)
         rfac=trfac/(r*r*r)
         bwindow = qwindow*(1.0+r**broad)
         sum=0.
         q1=xin(1)
         j1q1r=j1qr(q1,r)
         do j=1,lptin-1
            q2=xin(j+1)
            j1q2r=j1qr(q2,r)
            bwidth=0.5*(q1+q2)*bwindow
c Form integral of (Q sin Qr)/r from Q1 to Q2
            sum=sum+yin(j)*(j1q2r-j1q1r)
     1*pofr(bwidth)
            q1=q2
            j1q1r=j1q2r
          end do
         yout(i)=rfac*sum
         r1=r2
       end do
      yout(lptout)=0.0
      return
      end

      function j1qr(q,r)

      real*8 j1qr,q,r,qr

      qr=q*r
      j1qr = (sin(qr)-qr*cos(qr))

      return
      end

      function pofr(delr)

      if (delr.le.0.0) then
         pofr=1.0
      else
         delr3=3.0/(delr*delr*delr)
         pofr=delr3*(sin(delr)-delr*cos(delr))
      endif

      return
      end
