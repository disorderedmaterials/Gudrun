C     Last change:  GC   15 Jun 2005   10:00 pm
c----6---1---------2---------3---------4---------5---------6---------7--
c----6---1---------2---------3---------4---------5---------6---------7--
c
c      RRRRR   EEEEEE   GGGGG  DDDD    4         CCCCC
c      R    R  E       G       D   D   4        C     
c      R    R  E      G        D    D  4       C      
c      RRRRR   EEEE   G   GGG  D    D  4   4   C      
c      R  R    E      G     G  D    D  444444  C      
c      R   R   E       G    G  D   D       4    C     
c      R    R  EEEEEE   GGGG   DDDD        4     CCCCC
c
c----6---1---------2---------3---------4---------5---------6---------7--
c----6---1---------2---------3---------4---------5---------6---------7--
c
c  Compiling in d4lnx: g77 -ofile.exe -fno-silent regd4c31.for
      program regd4c
c
c  NM_NUM: Maximum number of numors
c  NM_CHA: Maximum number of channels in spectra
c  NM_ARG: Maximum number of arguments in command line
c
      parameter (NM_NUM=1000) 
      parameter (NM_CHA=1800)
      parameter (NM_ARG=21)
c
c  angle: 
c	vector of NM_CHA elements containing the scattering angle
c	(2theta in degrees) of the final spectra
c  spect:
c  	matrix 10 x NM_CHA. It represents 10 vectors containing
c	the number of counts for each of 9 detector banks. The
c	tenth vector contains the whole summed spectrum, i.e.,
c	the 9 detectors together
c  sigma:
c  	matrix 10 x NM_CHA. It represents 10 vectors containing
c	the error for each of 9 detector banks. The tenth vector 
c	contains the whole summed error bars, i.e., the 9 detectors 
c	together.
c
      real *8 angle(NM_CHA),spect(10,NM_CHA),sigma(10,NM_CHA)
      real *8 cont(10,NM_CHA)
c
c  raw:
c	3D matrix containing all the raw data NM_NUM x 9 x 64
c	They are organized for each numor in detectors (9) and
c	cells (64)
c  releff:
c	Matrix 9x64 containing the relative efficiency for each
c	detector and cell. This is a real number for which the
c	measurement of each cell is divided for.
c  dec:
c	9 elemnts vector containing the small correction that must
c	be applied for the angular position of each detector. For
c	the 1st cell of detector d the angular position is
c		pos1(d)=2theta0+(d-1)*15+dec(d)
c
      real *8 releff(9,64),raw(NM_NUM,9,64),dec(9),xin(64),yin(64)
      real *8 monit(NM_NUM),mtime(NM_NUM),cell1(NM_NUM),norm(NM_NUM)
c
c  angbin
c       3-elements vector containing the variables of the angular
c       binning.
c               angbin(1): initial angle
c               angbin(2): final angle
c               angbin(3): channel width
c
      real *8 angbin(3),wlength,zero,tmoni,ttime,param(85)
      real *8 xou(NM_CHA),you(NM_CHA),norma
      integer *4 h,i,j,k,l,m,n,ier,iunit,ndat,klong,tipo
      integer *4 iargc,narg,ntnum,ntchan,nlong
      integer *2 ipar(600)
      character *40 arg(NM_ARG),rep,regfile,regori,getfile
      character *6 numor(NM_NUM)
      character *100 comen(25)
      character *80 linea(15)
      character *16 parlbl(85)

      open (7,file="regd4c.log",status="unknown")
c
c  The function iargc gives the number of arguments in the command
c  line. 
c
      narg=iargc()
  20  continue
      if (narg.le.0) then
c----6---1---------2---------3---------4---------5---------6---------7--
      write (6,*)' Syntax: '
      write (6,*)'        regd4c [-h -9 -o filout -e fileff -d fildec -i
     # filign -b fillis -c filcor -r datadir'
      write (6,*)'               -a iangle fangle delta -z zeroangle -l 
     #lambda -n norm -g fileget] numorini '
      write (6,*)'               [numorfin]'
      write (6,*)
      write (6,*)'----------------'
      write (6,*)' Default values'
      write (6,*)'    -o numorfin.reg    -e effd4c.eff    -d d4cdec.dec'
      write (6,*)'    -a 0 140 0.125     -z 0             -l 0'
      write (6,*)'    -n 1000000 '
      write (6,*)'    -r /hosts/d4/users/data/ '
      write (6,*)'    numorfin: numorini'
      write (6,*)'----------------'
      write (6,*)
      stop
c----6---1---------2---------3---------4---------5---------6---------7--
      elseif ((narg+1).gt.NM_ARG) then
        write(6,*)
        write(6,*)"-----------------------------------------------"
        write(6,*)" ERROR - ERROR - ERROR - ERROR - ERROR - ERROR "
        write(6,*)"  Number of arguments greater than the maximum "
        write(6,*)"  allowed. Max = ",NM_ARG-1
        write(6,*)"-----------------------------------------------"
        write(6,*)
        stop
      endif
c
c  The subroutine getarg puts the argument i at the i+1 position of the
c  vector arg, which must have a dimension of narg+1, because the 0th
c  element is the name of the executable itself. This is the reason for
c  the error condition just above.
c
      do i=0,narg
        call getarg(i,arg(i+1))
      enddo
      
      write (7,*)" >>>>> Calling INPAR "
      call inpar(narg,arg,nlong,regori,releff,dec,ncor,angbin,
     # zero,wlength,ntnum,numor,norma,rep,getfile,comen,ier)
c      write (6,*)'ier',ier
c      write (6,*)' a la salida de inpar, nnormon ',nnormon
      if (ier.lt.0) then
        narg=0
        goto 20
      endif

      if (ntnum.gt.NM_NUM) then
        write(6,*)
        write(6,*)"-----------------------------------------------"
        write(6,*)" ERROR - ERROR - ERROR - ERROR - ERROR - ERROR "
        write(6,*)"  Total number of numors greater than the maximum "
        write(6,*)"  allowed. Max = ",NM_NUM
        write(6,*)"-----------------------------------------------"
        write(6,*)
        stop
      endif

      ntchan=(angbin(2)-angbin(1))/angbin(3)+1
c *** Initialisation
      do i=1,ntchan
        do j=1,10
           cont(j,i)=0.d0        ! counter
           spect(j,i)=0.d0       ! spectra
           sigma(j,i)=0.d0       ! exptal errors
        enddo
      enddo
      
      tmoni=0.0d0 ! Total monitor counts
      ttime=0.0d0 ! Total measuring time

      do n=1,ntnum
        do iunit=6,7
          write (iunit,*)" Reading numor ",numor(n)
c	  write (6,*)'rep ',rep
        enddo
        write (7,*)" >>>>> Calling READNUM "
        call readnum(n,numor,rep,param,raw,comen,linea,parlbl)
c        write (7,*)" >>>>> Exiting READNUM "

        if (getfile.ne.' ') then
        write (7,*)" >>>>> Calling GETALL "
	write (6,*)"getfile: ",getfile
          if (n.eq.1) then
            wlength=-1.0d0
            open (12,file=getfile,status='OLD')
            read (12,*)regori
            write (6,'(a)')regori
	    
	    write (6,'(a,a)')"Output get: ",regori
            read (12,*)tipo
	    write (6,*)'tipo',tipo
            if (tipo.eq.6) read (12,*)ipar(1),(ipar(i),i=2,ipar(1))
            close (12)
            open (11,file=regori,status='UNKNOWN')
            write (11,'(a,a)')'# ',regori
          endif
          call getall(tipo,n,numor,param,linea,raw,parlbl,ipar)
          goto 30
        endif

        write (7,*)" >>>>> Calling DEADTIME "
        call deadtime(n,param(1),raw,param(2))
	monit(n)=param(1)
        mtime(n)=param(2)
c  norm(n) is the normalisation variable. Here it is possible to
c  choose the normalisation method:
c     param(1): monitor     param(2): time     param(40): detector
c     param(31): det1       param(32): det2    param(33): det3
c     param(34): det4       param(35): det5    param(36): det6
c     param(37): det7       param(38): det8    param(39): det9
	norm(n)=param(1)
        if (norma.lt.0.d0) norm(n)=param(2)
        tmoni=tmoni+monit(n)
        ttime=ttime+mtime(n)
c   cell1 is the angular position of first cell of first detector
c  Here it could be possible to change the position of the detector
c  between the required (param(12)) and real (param(11)) angular
c  position.
        cell1(n)=param(11)-zero    ! Zero-angle correction: angle=raw-zero
   30 continue
      enddo
      if (wlength.lt.(0.0d0)) then
        close (11)
        write (6,*)'----> ',regori
        goto 15
      endif

      write (comen(12),'(a45,2(g20.9))')
     #"# Integrated monitor (counts) and time (sec): ",tmoni,ttime
c
c  Creating the spectra. In the first if-block the case of several
c  numors (more than 1) is considered. In the second one, the special
c  case of 1 numor is treated. In this case, there is no binning and
c  the real position of each cell is used as angle.
c
      if (ntnum.ge.2) then
        do n=1,ntnum
          do i=1,9
            write (7,'(a,i2)')' ----- Binning detector ',i
            do j=1,64
	      xin(j)=cell1(n)+(i-1)*15.0d0+dec(i)+(j-1)*0.125d0
	      yin(j)=0.0d0
              if (releff(i,j).gt.0.d0) yin(j)=raw(n,i,j)/releff(i,j)
              call eqbin(64,xin,yin,0.125d0,angbin,angle,you,ndat,ier)
	    enddo
	    do k=1,ndat
	      if (you(k).ne.0.d0) then
	        spect(i,k)=spect(i,k)+you(k)
	        cont(i,k)=cont(i,k)+norm(n)
	        spect(10,k)=spect(10,k)+you(k)
	        cont(10,k)=cont(10,k)+norm(n)
	      endif
	    enddo
	  enddo
        enddo
      else
        do i=1,9
          write (7,'(a,i2)')' ----- Positioning detector ',i
	  ndat=64*9  ! ndat=576
	  do j=1,64
	    xin(j)=cell1(1)+(i-1)*15.0d0+dec(i)+(j-1)*0.125d0
	    yin(j)=0.0d0
            if (releff(i,j).gt.0.d0) yin(j)=raw(1,i,j)/releff(i,j)
	  enddo
	  do j=1,64
	    k=64*(i-1)+j
     	    angle(k)=xin(j)
	    spect(i,k)=spect(i,k)+yin(j)
	    cont(i,k)=cont(i,k)+norm(1)
	    spect(10,k)=spect(10,k)+yin(j)
            cont(10,k)=cont(10,k)+norm(1)
	  enddo
	enddo
      endif	
        
c Normalisation and error evaluation
      do i=1,10
        do k=1,ndat
	  if (spect(i,k).ne.0.d0) then
	    sigma(i,k)=spect(i,k)/cont(i,k)*dabs(norma)
     #	             *dsqrt(1.0d0/spect(i,k)+1.0d0/cont(i,k))
	    spect(i,k)=spect(i,k)/cont(i,k)*dabs(norma)
	  endif
	enddo
      enddo
      
      regfile=regori
      klong=longit(regori)
c
c  Writing the output file (regfile). Again two cases must be considered.
c  If it is only 1 numor, then there is no binning and the real angular
c  position for each cell is taken into account and a fourth column is
c  written representing the cell which has counted each position,
c  according with this formula: col4=100*ndet+ncell
c
      write (7,*)" ----- Writing output file "
      n=10
      if (nlong.eq.1) n=1
      do i=10,n,-1
        if (i.ne.10) write (regfile(klong:klong),'(i1)')i
        open (1,file=regfile,status="unknown")
        write (6,*)'----> ',regfile
        do j=1,19
          write (1,'(a100)')comen(j)
        enddo
	h=0
	m=1
        do k=1,ndat
	  if (ntnum.eq.1) then
            h=h+1
            write (1,11)angle(k),spect(i,k),sigma(i,k),m*100+h
            l=64*(k/64)-k
	    if (l.eq.0) then
	      m=m+1
	      h=0
	      write (1,*)
	    endif
	  else
	    if (spect(i,k).ne.0.d0) 
     #        write (1,10)angle(k),spect(i,k),sigma(i,k)
          endif
	enddo
        close (1)
      enddo

      regfile=regori
      if (wlength.ne.0.d0) then
        write (7,*)" ----- Writing output file "
        n=10
        if (nlong.eq.1) n=1
        do i=10,n,-1
          if (i.ne.10) then
            write (regfile(klong:klong),'(i1)')i
            regfile(klong+1:klong+2)=".q"
          else
            regfile(klong:klong+1)=".q"
          endif
          open (1,file=regfile,status="unknown")
          write (6,*)'----> ',regfile
          do j=1,19
            write (1,'(a100)')comen(j)
          enddo
	  h=0
	  m=1
          do k=1,ndat
	    if (ntnum.eq.1) then
              h=h+1
              write (1,11)atoq(angle(k),wlength),spect(i,k),sigma(i,k)
     #                    ,m*100+h
              l=64*(k/64)-k
	      if (l.eq.0) then
	        m=m+1
	        h=0
	        write (1,*)
	      endif
	    else
	      if (spect(i,k).ne.0.d0)
     #          write (1,10)atoq(angle(k),wlength),spect(i,k),sigma(i,k)
            endif
	  enddo
          close (1)
        enddo
      endif

      write (7,*)
      write (7,*)" ++++   Heading   ++++"
      do i=1,19
         write(7,'(a100)') comen(i)
      enddo


   10 format (1x,f7.3,1x,g18.9,1x,g11.4)
   11 format (1x,f7.3,1x,g18.9,1x,g11.4,1x,i4)
   15 continue
      close(7)
      write (*,*)
      write (*,*)"See regd4c.log to have information about the run"
      write (*,*)
      stop
      end
c----6---1---------2---------3---------4---------5---------6---------7--
c
c      EEEEEE  N     N  DDDD    
c      E       NN    N  D   D   
c      E       N N   N  D    D  
c      EEEE    N  N  N  D    D  
c      E       N   N N  D    D  
c      E       N    NN  D   D   
c      EEEEEE  N     N  DDDD    
c
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c       GGGGG  EEEEEE  TTTTTTT    AAA    L       L      
c      G       E          T      A   A   L       L      
c     G        E          T     A     A  L       L      
c     G   GGG  EEEE       T     AAAAAAA  L       L      
c     G     G  E          T     A     A  L       L      
c      G    G  E          T     A     A  L       L      
c       GGGG   EEEEEE     T     A     A  LLLLLL  LLLLLL 
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine getall(tipo,n,numor,param,linea,raw,parlbl,ipar)
      parameter (NM_NUM=1000)
      integer *4  i,j,k,n,tipo
      integer *2  ipar(600)
      real *8 param(85),raw(NM_NUM,9,64)
      real *8 elapsed(NM_NUM),dia(NM_NUM),hora,minu,segu,cnum
      character *80 linea(15)
      character *100 comen(25)
      character *6 numor(NM_NUM)
      character *40 getfile
      character *16 parlbl(85)
      character *40 ti

      k=0
      do j=1,80
        if (linea(13)(j:j).eq.":") then
	  k=k+1
	  if (k.eq.2) goto 10
	endif
      end do
   10 continue
      ti=linea(13)(j-2:j-1)
      hora=cnum(ti)
      ti=linea(13)(j+1:j+2)
      minu=cnum(ti)
      ti=linea(13)(j+4:j+5)
      segu=cnum(ti)
      elapsed(n)=hora+minu/60.d0+segu/3600.d0
      if ((n.ne.1).and.(elapsed(n).lt.elapsed(n-1)))
     # elapsed(n)=elapsed(n)+24.d0
      write (6,*)linea(13)
      write (6,*)hora,minu,segu
      write (6,*)elapsed(n)

      if (tipo.eq.1) then       ! Tipo=1 heading
        write (11,'(a)')"----------------------------------------------"
	write (11,'(a,a8)')"Numor ",numor(n)
        write (11,'(a)')linea(6)
        do i=9,13
          write (11,'(a)') linea(i)
        enddo
      elseif (tipo.eq.2) then   ! Tipo=2 numor,etime,time,moni,tdet
        if (n.eq.1)
     #   write (11,'(a8,4a16)')
     #   '#  Numor','ElapsedTime',parlbl(2),parlbl(1),parlbl(40)
        write (11,'(a8,4(1pe16.8))')
     #   numor(n),elapsed(n),param(2),param(1),param(40)
      elseif (tipo.eq.3) then   ! Tipo=3 numor,etime,Tsam,Treg,Tset
        if (n.eq.1)
     #   write (11,'(a8,4a16)')
     #   '#  Numor','ElapsedTime',parlbl(41),parlbl(42),parlbl(43)
        write (11,'(a8,4(1pe16.8))')
     #   numor(n),elapsed(n),param(41),param(42),param(43)
      elseif (tipo.eq.4) then   ! Tipo=4 numor,etime,time,det1,...,det9,tdet
        if (n.eq.1)
     #   write (11,'(a8,12a16)')
     #   '#  Numor','ElapsedTime',parlbl(1)
     #    ,parlbl(31),parlbl(32),parlbl(33)
     #    ,parlbl(34),parlbl(35),parlbl(36)
     #    ,parlbl(37),parlbl(38),parlbl(39),parlbl(40)
        write (11,'(a8,12(1pe16.8))')
     #   numor(n),elapsed(n),param(1),param(31),param(32),param(33)
     #   ,param(34),param(35),param(36),param(37),param(38)
     #   ,param(39),param(40)
      elseif (tipo.eq.5) then
c  Tipo=5 numor,etime,lambda,2tread,2treal,BSC,BSD,Vol1,Vol2
        if (n.eq.1)
     #   write (11,'(a8,12a16)')
     #   '#  Numor','ElapsedTime',parlbl(5),parlbl(11),parlbl(12)
     #    ,parlbl(24),parlbl(25),parlbl(22),parlbl(23)
        write (11,'(a8,12(1pe16.8))')
     #   numor(n),elapsed(n),param(5),param(11),param(12)
     #   ,param(24),param(25),param(22),param(23)
      elseif (tipo.eq.6) then
c  Continuar aquí con la rutina getall
        if (n.eq.1)
     #   write (11,'(a8,20a16)')
     #   '#  Numor','ElapsedTime',(parlbl(ipar(j)),j=2,ipar(1))
        write (11,'(a8,20(1pe16.8))')
     #   numor(n),elapsed(n),(param(ipar(j)),j=2,ipar(1))
      endif

      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--


c----6---1---------2---------3---------4---------5---------6---------7--
c
c      EEEEEE    QQ     BBBBBB   III  N     N  
c      E        Q  Q    B     B   I   NN    N  
c      E       Q    Q   B     B   I   N N   N  
c      EEEE    Q    Q   BBBBBB    I   N  N  N  
c      E       Q  Q Q   B     B   I   N   N N  
c      E        Q  QQ   B     B   I   N    NN  
c      EEEEEE    QQ  Q  BBBBBB   III  N     N  
c
c----6---1---------2---------3---------4---------5---------6---------7--
c  This subroutine makes a binning of the input vector yin using xin
c      subroutine eqbin(nin,xin,yin,derr,dxin,
c     # xou1,xou2,dxou,xou,you,sig,nou,ier)
      subroutine eqbin(nin,xin,yin,dxin,bin,xou,you,nou,ier)
      parameter (NM_CHA=1800)
      real *8 xin(64),yin(64),dxin
      real *8 xou(NM_CHA),you(NM_CHA),hits(NM_CHA)
      real *8 bin(3),f,x1,x2,fy
      integer *4 nou,nin,i,j,i1,i2,ier
C----6---1---------2---------3---------4---------5---------6---------7-2
      ier=0
      if (bin(3).le.0.d0) then
         ier=-1
         write (7,*)ier,' --> Binning width is zero or negative'
         return
      endif
C *****  Number of equal bins of the output vectors (xou,you,sig)
      nou=int((bin(2)-bin(1))/bin(3))+1
      if ((nou.lt.1).or.(nou.gt.NM_CHA)) then
         ier=-2
         write (7,*)ier,' --> Dimension problem in output binning'
         return
      endif
      if ((nin.lt.1).or.(nin.gt.64)) then
         ier=-3
         write (7,*)ier,' --> Dimension problem in input data'
         return
      endif
      
C *****  Initialization of output vectors and hit counter
      do i=1,nou
         xou(i)=bin(1)+(i-1)*bin(3)
         you(i)=0.d0
         hits(i)=0.d0
      enddo

C *****  Here the subroutine starts
C        A cycle over all the input data is done, checking its 
C        contribution to each output bin
C
      if (dxin.lt.0.d0) then
        ier=-4
        write (7,*)ier,' --> Input binning width is zero or negative'
        return
      endif

C----6---1---------2---------3---------4---------5---------6---------7-2

      f=dxin/bin(3)
      do j=1,nin
C  --- These are the limits of the j-th channel
C      x1,x2 are the limits in abcise coordinates whereas
C      i1,i2 are the corresponding limits in position
         x1=xin(j)-dxin/2.d0
         x2=xin(j)+dxin/2.d0
         i1=1+int((x1-bin(1)+bin(3)/2.d0)/bin(3))
         i2=1+int((x2-bin(1)+bin(3)/2.d0)/bin(3))

C  --- Now the different possibilities must be considered
C
C  --- If i1 is equal to i2, this channel is completely included in
C      one bin
C
        if (i1.eq.i2) then
           call addone(i1,you,hits,f,yin(j),1.d0)
        else
c	   write (7,*)'i1 es distinto de i2',i1,i2

           fy=(xou(i1)+bin(3)/2.d0-x1)/dxin
           call addone(i1,you,hits,f,yin(j),fy)

           fy=(x2-xou(i2)+bin(3)/2.d0)/dxin
           call addone(i2,you,hits,f,yin(j),fy)
	   
           if (i2.gt.(i1+1)) then
             fy=bin(3)/dxin
             do i=i1+1,i2-1
               call addone(i,you,hits,f,yin(j),fy)
             enddo
           endif
        endif
      enddo
      
C----6---1---------2---------3---------4---------5---------6---------7-2
C ****  Normalisation
      do i=1,nou
         if (hits(i).le.0.d0) then
            you(i)=0.d0
         else
c	    write (7,*)'*hits',i,hits(i),xou(i)
c	    write (7,*)' -----> ',you(i),sig(i)
            you(i)=you(i)/hits(i)
         endif
      enddo

      return
      end
C----6---1---------2---------3---------4---------5---------6---------7-2

c----6---1---------2---------3---------4---------5---------6---------7--
c
c       AAA    DDDD    DDDD      OOO    N     N  EEEEEE 
c      A   A   D   D   D   D    O   O   NN    N  E      
c     A     A  D    D  D    D  O     O  N N   N  E      
c     AAAAAAA  D    D  D    D  O     O  N  N  N  EEEE   
c     A     A  D    D  D    D  O     O  N   N N  E      
c     A     A  D   D   D   D    O   O   N    NN  E      
c     A     A  DDDD    DDDD      OOO    N     N  EEEEEE 
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine addone(i,you,hits,f,y,fy)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (NM_CHA=1800)
      real *8 you(NM_CHA),hits(NM_CHA)
      real *8 f,y,fy
      integer *4 i
C----6---1---------2---------3---------4---------5---------6---------7-2
C ****  you(i) is the bin to where counts must be added
C
C ****  hits(i) is a counter which counts the 'number' of times that the
C       i-th bin has been modified. This number is required for 
C       normalization
C
C ****  f is the fraction of the output bin that the input width 
C       represents
C
C ****  y is the histogram value of the input data
C
C ****  fy is the fraction of the input histogram which contributes to
C       the given bin
C----6---1---------2---------3---------4---------5---------6---------7-2
      if (y.gt.0.d0) then
        you(i)=you(i)+fy*y
        hits(i)=hits(i)+fy*f
      endif
      return
      end
C----6---1---------2---------3---------4---------5---------6---------7-2

c----6---1---------2---------3---------4---------5---------6---------7--
c
c      III  N     N  PPPPP     AAA    RRRRR   
c       I   NN    N  P    P   A   A   R    R  
c       I   N N   N  P    P  AA   AA  R    R  
c       I   N  N  N  PPPPP   AAAAAAA  RRRRR   
c       I   N   N N  P       AA   AA  R  R     
c       I   N    NN  P       AA   AA  R   R    
c      III  N     N  P       AA   AA  R    R  
c
c----6---1---------2---------3---------4---------5---------6---------7--
c  Subroutine inpar
c    Reads parameters from unix command line
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine inpar(narg,arg,nlong,regfile,releff,dec,ncor,angbin,
     # zeroangle,xlam,ntnumor,numor,norma,repert,getfile,comen,ier)

      parameter (NM_ARG=21)
      parameter (NM_NUM=1000) 
      real *8 releff(9,64),dec(9),xlam,zeroangle,angbin(3),cnum,norma
      integer *4 ncor(9,64),ntnumor,narg,nfin,nbeg
      integer *4 numori,numorf,nign,nextra,numex
      integer *4 h,i,j,k,l,m,n,ier,ini
      integer *4 nini(NM_ARG),ignn(NM_NUM),nefftyp
      character *40 arg(NM_ARG)
      CHARACTER *24 the_date
      character *100 comen(25),linput
      character *6 numor(NM_NUM),nkilled
      character *40 regfile,efffile,corfile,decfile,ignfile,repert
      character *40 lisfile,getfile


      ier=-1
      nlong=0
      xlam=0.d0
      zeroangle=0.d0
      norma=1000000.d0
      angbin(1)=0.d0
      angbin(2)=140.d0
      angbin(3)=0.125d0
      regfile=' '
      efffile='effd4c.eff'
      corfile=' '
      decfile='dec.dec'
      ignfile='ign.ign'
      getfile=' '
      lisfile=' '
      repert='/net/serdon/illdata/data/d4/'
      do i=1,25
        comen(i)=" "
      enddo

c nini is a set of flags which distinguish arguments not read (nini=1) from those
c already read

      do i=1,NM_ARG  ! Here these flags are initialized to 0
        nini(i)=0
      enddo
      do i=2,narg+1 ! And here those corresponding to real arguments are switched to 1
        nini(i)=1
      enddo
      
      do i=2,narg+1
        if (arg(i)(1:2).eq.'-h') call help
      enddo
      
      do i=2,narg+1
        if (arg(i)(1:1).eq.'-') then
          nini(i)=0
          if (arg(i)(2:2).ne.'9') nini(i+1)=0
          if (arg(i)(2:2).eq.'9') nlong=1
          if (arg(i)(2:2).eq.'a') then
            angbin(1)=cnum(arg(i+1))
            angbin(2)=cnum(arg(i+2))
            angbin(3)=cnum(arg(i+3))
            nini(i+2)=0
            nini(i+3)=0
          endif
          if (arg(i)(2:2).eq.'l') xlam=cnum(arg(i+1))
          if (arg(i)(2:2).eq.'z') zeroangle=cnum(arg(i+1))
          if (arg(i)(2:2).eq.'n') then
            write (6,*)arg(i),arg(i+1)
            norma=cnum(arg(i+1))
            write (6,*)norma
          endif
c          if (arg(i)(2:2).eq.'t') nnormon=-cnum(arg(i+1))
          if (arg(i)(2:2).eq.'o') regfile=arg(i+1)
          if (arg(i)(2:2).eq.'b') lisfile=arg(i+1)
          if (arg(i)(2:2).eq.'e') efffile=arg(i+1)
          if (arg(i)(2:2).eq.'d') decfile=arg(i+1)
          if (arg(i)(2:2).eq.'i') ignfile=arg(i+1)
          if (arg(i)(2:2).eq.'c') corfile=arg(i+1)
          if (arg(i)(2:2).eq.'g') getfile=arg(i+1)
          if (arg(i)(2:2).eq.'r') repert=arg(i+1)
        endif
      enddo
      
      ini=0
      do i=2,narg+1
        ini=ini+nini(i)
      enddo
      
      if ((ini.lt.1).or.(ini.gt.2)) return
      do i=2,narg+1
        j=nini(i)*ini
        if (j.gt.0) then 
          nbeg=i
          goto 20
        endif
      enddo
   20 continue
      nfin=nbeg
      if (j.eq.2) nfin=nbeg+1

      call ignumor(ignfile,nign,ignn)
      
      numori=cnum(arg(nbeg))
      numorf=cnum(arg(nfin))

      ntotal=numorf-numori+1

      ntnumor=0
      comen(3)="# All numors in the range"
      do i=numori,numorf
        do j=1,nign
          if (i.eq.ignn(j)) then
            call numtocha(i,nkilled)
            write (7,*)' Numor ',nkilled,' not included'
            write (6,*)' Numor ',nkilled,' not included'
            write (comen(3),'(a32,a40)')
     #       "# There are ignored numors. See ",ignfile
            goto 30
          endif
        enddo
        ntnumor=ntnumor+1
        call numtocha(i,numor(ntnumor))
   30   continue
      enddo
      
      if (lisfile.ne.' ') then
        open (97,file=lisfile,status="old")
	read (97,*)nextra
	do j=1,nextra
	  read(97,*)numex
	  ntnumor=ntnumor+1
          call numtocha(numex,numor(ntnumor))
	enddo
      endif
      
      do i=2,narg+1
        if (arg(i)(1:2).eq.'-o') goto 40
      enddo
      do i=1,6
        regfile(i:i)=numor(ntnumor)(i:i)
      enddo
      regfile(7:10)='.reg'
      regfile(11:40)='                              '
   40 continue
   
      nefftyp=3
      if (efffile.eq.'1') nefftyp=1
      if (efffile.eq.'2') nefftyp=2
      call eficiencia(nefftyp,releff,efffile)

      call decal(decfile,dec)
c      if (corfile.ne.' ') call cellcor(corfile,ncor)
      
c      write (6,*)ntnumor,' ',numor(1),' ',numor(ntnumor)
c      write (6,*)'Output file: -',regfile,'-',longit(regfile)
c      write (6,*)'Decala file: ',decfile,longit(decfile)
c      write (6,*)'Effic. file: ',efffile,longit(efffile),nefftyp
c      write (6,*)'Ignore file: ',ignfile,longit(ignfile)
c      write (6,*)'Correc file: ',corfile,longit(corfile)
c      write (6,*)'Extran file: ',lisfile,longit(lisfile)
c      write (6,*)'Data direct: ',repert,longit(repert)
c      write (6,*)'nlong ',nlong
c      write (6,*)'angles ',(angbin(i),i=1,3)
c      write (6,*)'zeroan ',zeroangle
c      write (6,*)'lambda ',xlam

c
c line  1: Filename
c line  2: Numors
c line  3: Ignored numors
c line  4: Command line
c line  5: Efficiency file
c line  6: Decalage file
c line  7: Sample
c line  8: User
c line  9: Date Time
c line 10: Corrected cells
c line 11: Normalisation method
c line 12: Integrated monitor and time
c line 13: zeroangle
c line 14: Angle formula
c line 15: Decalages
c line 16: Binning
c line 17: Wavelength
c line 18: Nothing
c line 19: Columns headings
c

      write (comen(1),'(a2,a40)')"# ",regfile
      write(comen(2),'(a20,a6,a4,a6)')"# regd4c for numors ",numor(1),
     #" to ",numor(ntnumor)
      write(comen(4),'(a16,a80)')"# Command line: ",
     # linput(narg,arg)
      if (nefftyp.eq.1) then
        comen(5)="# Efficiency 1 for all cells"
      elseif (nefftyp.eq.2) then
        comen(5)="# Efficiency 1 for all cells, except 1 and 64 (-1)"
      else
        write(comen(5),'(a19,a40)')"# Efficiency file: ",efffile
      endif
      write(comen(6),'(a17,a40)')"# Decalage file: ",decfile
      call fdate(the_date)
      write (comen(9),'(a8,a24)')"# Date: ",the_date
      if (corfile(1:1).eq.' ') then
         comen(10)="# Correction by software: no"
      else
         write(comen(10),'(a20,a40)')"# Correction file: ",corfile
      endif
      if (norma.gt.0.d0) then
         write(comen(11),'(a37,f14.1)')
     #    "# Normalisation by monitor (counts): ",norma
      elseif (norma.lt.0) then
         write(comen(11),'(a29,f14.1)')"# Normalisation by time (sec): "
     #    ,-norma
      else
        write (6,*)" Error in normalisation constant"
        stop
      endif
      write(comen(13),'(a20,f14.7)')"# zeroangle (deg) = ",zeroangle
      comen(14)="# angle = angle.raw - zeroangle [+ dec]"
      write(comen(15),'(a13,9(f7.3))')"# dec (deg)= ",(dec(i),i=1,9)
      write(comen(16),'(a11,3(f9.3))')"# Binning= ",(angbin(i),i=1,3)
      write(comen(17),'(a18,f14.7)')"# Wavelength (A)= ",xlam
      comen(18)="# "
      if (ntnumor.eq.1) then
        comen(19)="# Angle(deg)   Counts         Sigma     Cell"
      else
        if (xlam.eq.0) then
          comen(19)="# Angle(deg)   Counts         Sigma "
        else
          comen(19)="#  Q(1/A)      Counts         Sigma "
        endif
      endif

      ier=0
c      write (6,*)"ier en INPAR: ",ier 
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c      RRRRRR   EEEEEE    AAA    DDDD    N     N  U     U  MM   MM
c      R     R  E        A   A   D   D   NN    N  U     U  M M M M
c      R     R  E       A     A  D    D  N N   N  U     U  M  M  M
c      RRRRRR   EEEE    AAAAAAA  D    D  N  N  N  U     U  M     M
c      R  R     E       A     A  D    D  N   N N  U     U  M     M
c      R   R    E       A     A  D   D   N    NN  U     U  M     M
c      R    R   EEEEEE  A     A  DDDD    N     N  UUUUUUU  M     M
c
c----6---1---------2---------3---------4---------5---------6---------7--
c Reads data from a numor file
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine readnum(num,numor,rep,param,cuentas,comen,linea,parlbl)
      parameter (NM_NUM=1000) 
      integer *4  i,j,det,num,ndir,ier
      real *8 param(85),cuentas(NM_NUM,9,64)
      character *80 linea(15)
      character *100 comen(25)
      character *6 numor(NM_NUM)
      character *40 rep,fichero
      character *16 parlbl(85)
      
      ier=0
  20  goto 25    
  15  continue
      ier=1
  25  continue
      if (ier.ne.0) then
        rep='/net/serdon/illdata/data-1/d4/'
	write (6,*)'Searching in: ',rep
      endif 
      ndir=longit(rep)
      fichero(1:ndir-1)=rep(1:ndir-1)
      fichero(ndir:ndir+5)=numor(num)
      fichero(ndir+6:40)='                                           '

      open(1,file=fichero,status="old",err=15)

      do i=1,13
        read (1,'(a)')linea(i)
      enddo
      do i=1,3
        read (1,'(a)')linea(15)
      enddo
      comen(7)=linea(9)(11:80)
      comen(8)=linea(10)(13:80)
      comen(7)(1:10)="# Sample: "
      comen(8)(1:8) ="# User: "

      do i=1,8
        read (1,'(5a16)')(parlbl(j+5*(i-1)),j=1,5)
      enddo

      do i=1,8
        read (1,*)(param(j+5*(i-1)),j=1,5)
      enddo

      do i=1,3
        read (1,'(a)')linea(15)
      enddo
      do i=9,12
        read (1,'(5a16)')(parlbl(j+5*(i-1)),j=1,5)
      enddo
      do i=9,12
        read (1,*)(param(j+5*(i-1)),j=1,5)
      enddo

      do i=1,3
        read (1,'(a)')linea(15)
      enddo
      do i=13,17
        read (1,'(5a16)')(parlbl(j+5*(i-1)),j=1,5)
      enddo
      do i=13,17
        read (1,*)(param(j+5*(i-1)),j=1,5)
      enddo

      do det=1,9
        do i=1,4
          read (1,'(a)')linea(15)
        enddo
        do i = 1,6
           read (1,*)(cuentas(num,det,j+10*(i-1)),j=1,10)
        enddo
        read (1,*)(cuentas(num,det,60+j),j=1,4)
        linea(15)="Detector "
        write (7,'(a9,i2,a6)')linea(15),det," read"
      enddo
      close (1)

  10  format (i6,i12,2x,1pe9.2,2x,i3,2x,1pe14.7)
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--


c----6---1---------2---------3---------4---------5---------6---------7--
c
c N     N  U     U  MM   MM  TTTTTTT   OOO      CCCCC  H     H    AAA     
c NN    N  U     U  M M M M     T     O   O    C       H     H   A   A    
c N N   N  U     U  M  M  M     T    O     O  C        H     H  A     A
c N  N  N  U     U  M     M     T    O     O  C        HHHHHHH  AAAAAAA   
c N   N N  U     U  M     M     T    O     O  C        H     H  A     A    
c N    NN  U     U  M     M     T     O   O    C       H     H  A     A    
c N     N  UUUUUUU  M     M     T      OOO      CCCCC  H     H  A     A
c
c----6---1---------2---------3---------4---------5---------6---------7--
c     Converts an integer in a 6-character string
c
c  numero: integer
c  numor: string of 6 characters
c
      subroutine numtocha(numero,numor)
      character *6 numor
      integer *4 numero,n(6),i,k
      
      k=numero
      
      n(1)=numero/100000
      numero=numero-n(1)*100000
      n(2)=numero/10000
      numero=numero-n(2)*10000
      n(3)=numero/1000
      numero=numero-n(3)*1000
      n(4)=numero/100
      numero=numero-n(4)*100
      n(5)=numero/10
      n(6)=numero-n(5)*10
      do i=1,6
        write (numor(i:i),'(i1)')n(i)
      enddo
      numero=k
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c         CCCCC  N     N  U     U  MM   MM 
c        C       NN    N  U     U  M M M M   
c       C        N N   N  U     U  M  M  M 
c       C        N  N  N  U     U  M     M   
c       C        N   N N  U     U  M     M    
c        C       N    NN  U     U  M     M    
c         CCCCC  N     N  UUUUUUU  M     M 
c
c----6---1---------2---------3---------4---------5---------6---------7--
c       This function converts a character string in a number
c----6---1---------2---------3---------4---------5---------6---------7--
      real *8 function cnum(a)
      character *40 a
      real *8 s,p(40)
      integer *2 f,nc,np,nt,i

      cnum=0.d0
      nc=1
      np=-1
      f=1
      do i=1,40
        if (a(i:i).eq.'-') then
          f=-1
          nc=2
          goto 10
        elseif (a(i:i).eq.'+') then
          nc=2
          goto 10
        elseif(a(i:i).eq.'.') then
          np=i
          goto 10
        elseif (a(i:i).eq.' ') then
          nt=i-1
          goto 20
        else
          p(i)=ichar(a(i:i))-48
        endif
 10     continue
      enddo
      return
 20   continue

      if (np.eq.-1) np=nt+1
      if (np.gt.nc) then
        if (np.gt.(nc+1)) then
          s=0.1d0
          do i=np-1,nc,-1
            s=s*10.d0
            cnum=cnum+s*p(i)
          enddo
        else
          cnum=p(nc)
        endif
      endif
      s=1.d0
      do i=np+1,nt
        s=s/10.d0
        cnum=cnum+s*p(i)
      enddo
      cnum=f*cnum

      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c      III   GGGGG  N     N  U     U  MM   MM    OOO    RRRRRR  
c       I   G       NN    N  U     U  M M M M   O   O   R     R 
c       I  G        N N   N  U     U  M  M  M  O     O  R     R 
c       I  G   GGG  N  N  N  U     U  M     M  O     O  RRRRRR  
c       I  G     G  N   N N  U     U  M     M  O     O  R  R    
c       I   G    G  N    NN  U     U  M     M   O   O   R   R   
c      III   GGGG   N     N  UUUUUUU  M     M    OOO    R    R  
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine ignumor(ignfile,nign,ignn)
      parameter (NM_NUM=1000) 
      character *40 ignfile
      character *80 linea
      integer *4 ignn(NM_NUM),nign,i,j

      nign=0
      write (7,*)"+++ Opening the file ",ignfile
      open (97,file=ignfile,status="old")
      do i=1,10000
        read(97,'(a)',end=20)linea
        if (linea(1:1).eq."#") goto 10
        do j=1,80
          if (linea(j:j).ne.' ') goto 40
        enddo
        goto 10
   40   backspace (97)
        nign=nign+1
        read (97,*)ignn(nign)
   10   continue
      enddo
   20 continue
      close (97)
      write (7,*)"+++ Closing the file ",ignfile

      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c     DDDD    EEEEEE    AAA    DDDD   TTTTTTT  III  MM   MM  EEEEEE 
c     D   D   E        A   A   D   D     T      I   M M M M  E      
c     D    D  E       AA   AA  D    D    T      I   M  M  M  E      
c     D    D  EEEE    AAAAAAA  D    D    T      I   M     M  EEEE   
c     D    D  E       AA   AA  D    D    T      I   M     M  E      
c     D   D   E       AA   AA  D   D     T      I   M     M  E      
c     DDDD    EEEEEE  AA   AA  DDDD      T     III  M     M  EEEEEE 
c
c----6---1---------2---------3---------4---------5---------6---------7--
c
c      This subroutine performs the dead-time correction on the monitor
c    and detector comptage.
c      The dead-time for the monitor is 7.2 microseg and for the detec-
c    tor is 7.17 microsec. These data are taken from the Internal ILL 
c    Report (ILL98FI15T) by H.E. Fischer, P. Palleau and D. Feltin.
c      The formula to correct the measured counting rate (nu_m) is as
c    follows:
c             nu=nu_m/(1-nu_m*tau)       where tau is the dead-time
c
c    Warning! The output is the same as the input, variable cmon (moni-
c    tor counts) and matrix counts are modified by calling this subrou-
c    tine.
c
c      Henry Fischer says that better values for monitor and detector 
c    dead-time should be 12 and 0 microsec, respectively.
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine deadtime(num,cmon,counts,temp)
      parameter (NM_NUM=1000)
      real *8 fmon,fdet,cmon,temp 
      integer *4 ncell,ndet,num
      real *8 counts(NM_NUM,9,64)
      
      fmon=12.0d-6/temp  ! Deadtime for monitor
      fdet=7.17d-6/temp  ! Deadtime for microstrips
      
      cmon=cmon/(1.0d0-fmon*cmon)
      do ndet=1,9
        do ncell=1,64
          counts(num,ndet,ncell)=counts(num,ndet,ncell)
     #              /(1.0d0-fdet*counts(num,ndet,ncell))
        enddo
      enddo

      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--
      
c----6---1---------2---------3---------4---------5---------6---------7--
c
c EEEEEE  FFFFFF III   CCCCC  III  EEEEEE  N     N    CCCCC  III   AAA    
c E       F       I   C        I   E       NN    N   C        I   A   A   
c E       F       I  C         I   E       N N   N  C         I  A     A  
c EEEE    FFFF    I  C         I   EEEE    N  N  N  C         I  AAAAAAA  
c E       F       I  C         I   E       N   N N  C         I  A     A  
c E       F       I   C        I   E       N    NN   C        I  A     A  
c EEEEEE  F      III   CCCCC  III  EEEEEE  N     N    CCCCC  III A     A  
c
c----6---1---------2---------3---------4---------5---------6---------7--
c  Subroutine eficiencia
c
c
c  type:  1- All efficiencies are 1 (default)
c         2- All 1, except the ends (-1)
c         3- Read an efficiency file (ficheff)
c  eff: matrix 9x64 containing the relative efficiencies
c  ficheff: efficiency file (only for type 3)
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine eficiencia(tipo,eff,ficheff)
      integer *4 tipo,det,cell,i,j
      real *8 eff(9,64),eff0
      character *12 ficheff
      character *80 linea

      write (7,*)"+++ By default, all the efficiencies are 1"
      do det=1,9
        do cell=1,64
          eff(det,cell)=1.d0
        enddo
      enddo
      if (tipo.eq.2) then
        write (7,*)"+++ ... but now first and last cells are killed"
        do det=1,9
          eff(det, 1)=-1.d0
          eff(det,64)=-1.d0
        enddo
      else if (tipo.eq.3) then
        write (7,*)"+++ ... but now the file ",ficheff," is used"
c  The first character of the filename should not be an space
   50   continue
        if (ficheff(1:1).eq.' ') then
          ficheff(1:11)=ficheff(2:12)
          goto 50
        endif
        open (97,file=ficheff,status="old")
        do i=1,10000
          read(97,'(a)',end=20)linea
          if (linea(1:1).eq."#") goto 10
          do j=1,80
            if (linea(j:j).ne.' ') goto 40
          enddo
          goto 10
   40     backspace (97)
          read (97,*)det,cell,eff0
          eff(det,cell)=eff0
   10     continue
        enddo
   20   continue
        close (97)
        write (7,*)"+++ Closing the file ",ficheff
      endif
      write (7,*)"*** Efficiencies"
      do cell=1,64
        write (7,30)cell,(eff(det,cell),det=1,9)
      enddo
      if (tipo.eq.1) then
        write (7,*)" All cells from all detectors are considered"
      else
        do det=1,9
	  do cell=1,64
            if (eff(det,cell).le.0) 
     #          write (6,*)" Bad cell",cell," in detector",det
          enddo
        enddo
      endif
   30 format (i3,9(g9.3))
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c         AAA   TTTTTTT   OOO      QQ   
c        A   A     T     O   O    Q  Q  
c       AA   AA    T    O     O  Q    Q  
c       AAAAAAA    T    O     O  Q    Q   
c       AA   AA    T    O     O  Q  Q Q     
c       AA   AA    T     O   O    Q  QQ    
c       AA   AA    T      OOO      QQ  Q  
c
c----6---1---------2---------3---------4---------5---------6---------7--
c       This function converts angles (in deg) into Q (in reciprocal
c     Angstroms), using a given wavelength (xlam).
c----6---1---------2---------3---------4---------5---------6---------7--
      real *8 function atoq(angle,xlam)
      real *8 xlam,pi,angle

      pi =4.0d0*datan(1.0d0)
      atoq=4.d0*pi/xlam*dsin(angle*pi/180.d0/2.0d0)

      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c      DDDD    EEEEEE    CCCCC    AAA    L      
c      D   D   E        C        A   A   L      
c      D    D  E       C        A     A  L      
c      D    D  EEEE    C        AAAAAAA  L      
c      D    D  E       C        A     A  L      
c      D   D   E        C       A     A  L      
c      DDDD    EEEEEE    CCCCC  A     A  LLLLLL 
c
c----6---1---------2---------3---------4---------5---------6---------7--
c       This subroutine reads the small corrections to the zeroangle
c     for each detector.
c
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine decal(decfile,dec)
      real *8 dec(9),deca
      integer *4 i,j,ndet
      character *40 decfile
      character *80 linea
      
      write (7,*)"+++ Opening the file ",decfile
      open (97,file=decfile,status="old")
      do i=1,10000
        read(97,'(a)',end=20)linea
        if (linea(1:1).eq."#") goto 10
        do j=1,80
          if (linea(j:j).ne.' ') goto 40
        enddo
        goto 10
   40   backspace (97)
        read (97,*)ndet,deca
        dec(ndet)=deca
   10   continue
      enddo
   20 continue
      close (97)
      write (7,*)"+++ Closing the file ",decfile

      write (7,*)"   Decalages "
      do ndet=1,9
        write (7,*)ndet,dec(ndet)
      enddo
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c       L        OOO    N     N    GGGGG  III TTTTTTT   
c       L       O   O   NN    N   G        I     T      
c       L      O     O  N N   N  G         I     T      
c       L      O     O  N  N  N  G   GGG   I     T      
c       L      O     O  N   N N  G     G   I     T       
c       L       O   O   N    NN   G    G   I     T       
c       LLLLLL   OOO    N     N    GGGG   III    T      
c
c----6---1---------2---------3---------4---------5---------6---------7--
c       This function obtains the length of a character string
c----6---1---------2---------3---------4---------5---------6---------7--
      integer *4 function longit(a)
      character *40 a
      integer *4 i
      
      do i=1,40
        if (a(i:i).eq.' ') goto 10
      enddo
      i=i-1
   10 longit=i
      
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--

c----6---1---------2---------3---------4---------5---------6---------7--
c
c      L       III  N     N  PPPPP   U     U  TTTTTTT
c      L        I   NN    N  P    P  U     U     T   
c      L        I   N N   N  P    P  U     U     T   
c      L        I   N  N  N  PPPPP   U     U     T   
c      L        I   N   N N  P       U     U     T    
c      L        I   N    NN  P       U     U     T    
c      LLLLLL  III  N     N  P       UUUUUUU     T   
c
c----6---1---------2---------3---------4---------5---------6---------7--
c       This function produces a string containing the input line
c----6---1---------2---------3---------4---------5---------6---------7--
      character *100 function linput(narg,arg)
      character *40 arg(20),argu
      integer *4 narg,ind,i,j,l

      linput=arg(1)
      ind=longit(linput)

      do i=1,narg
        argu=arg(i+1)
        l=longit(argu)
        linput(ind+1:ind+l)=argu(1:l)
        ind=ind+l
      enddo
      
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--


c----6---1---------2---------3---------4---------5---------6---------7--
c
c      H     H  EEEEEE  L       PPPPP 
c      H     H  E       L       P    P
c      H     H  E       L       P    P
c      HHHHHHH  EEEE    L       PPPPP 
c      H     H  E       L       P      
c      H     H  E       L       P      
c      H     H  EEEEEE  LLLLLL  P     
c
c----6---1---------2---------3---------4---------5---------6---------7--
c----6---1---------2---------3---------4---------5---------6---------7--
      subroutine help
      write (6,*)
      write (6,*)'------------- Instructions to use regd4c ------------'
      write (6,*)' Syntax: '
      write (6,*)'        regd4c [-h -9 -o filout -e fileff -d fildec -i
     # filign -b fillis -c filcor -r datadir'
      write (6,*)'               -a iangle fangle delta -z zeroangle -l 
     #lambda -n norm -g fileget] numorini '
      WRITE (6,*)'               [numorfin]'
      write (6,*)
      write (6,*)' numorini is a 6-digits integer representing the first
     # numor to be regrouped. It is the only'
      write (6,*)'          required argument.'
      write (6,*)' numorfin is a 6-digits integer representing the last
     #numor to be regrouped.' 
      write (6,*)
      write (6,*)' Arguments whithin [ ] are optional.'
      write (6,*)
      write (6,*)' Attention! Between switches and arguments, or between 
     # arguments, a space is always required.'
      write (6,*)'----------------'
      write (6,*)' Default values'
      write (6,*)'    -o numorfin.reg    -e effd4c.eff    -d d4cdec.dec'
      write (6,*)'    -a 0 140 0.125     -z 0             -l 0'
      write (6,*)'    -n 100000 '
      write (6,*)'    -r /hosts/d4/users/data/ '
      write (6,*)'    numorfin: numorini'
      write (6,*)'----------------'
      write (6,*)' Switches:'
      write (6,*)'  -h:   Show this help and stop'
      write (6,*)
      write (6,*)'  -9:   Produces a long output, i.e., creates one indi
     #vidual file for each detector.'
      write (6,*)'        The names of these files have the same root as
     # the oputput file, but 1,2,...,9 as last'
      write (6,*)'        character. By default, regd4c does not create 
     #these files.'
      write (6,*)
      write (6,*)'  -o:   Output filename (usually with extension .reg).
     # This character string will be used as'
      write (6,*)'        root for other output files. By default this f
     #ilename is numorfin with extension .reg'
      write (6,*)
      write (6,*)'  -e:   Efficiency file containing the relative effici
     #ency for each cell. By default regd4c'
      write (6,*)'        uses effd4c.eff.'
      write (6,*)
      write (6,*)'  -d:   File containing the angular shift for each det
     #ector. This is a small correction to the'
      write (6,*)'        supposed 7 degrees in between two adjacent det
     #ectors.'
      write (6,*)
      write (6,*)'  -i:   File containing the numors to be ignored. By d
     #efault regd4c does not ignore any numor.'
      write (6,*)
      write (6,*)'  -b:   File containing a list of extra numors to be r
     #egrouped. By default regd4c does not '
      write (6,*)'        include extra numors.'
      write (6,*)
      write (6,*)'  -c:   File containing the cells to be corrected by t
     #he charge effect, i.e., when one cell '
      write (6,*)'        count less and the cells beside have a higher 
     #comptage, but mantaining a good average.'
      write (6,*)'        By default regd4c does not correct any cell.'
      write (6,*)'        The cells that can be corrected must be in the
     # range (3,62). If cells 2 or 63 are '
      write (6,*)'        affected, they should be killed.'
      write (6,*)
      write (6,*)'  -a:   This switch is to define the angular histogram
     # for the output file. It requires three'
      write (6,*)'        arguments: first, last, and step angle. Defaul
     #t values are 0, 140 and 0.125 deg.'
      write (6,*)'        The recommended value for the step is 0.125 (o
     #ne cell).'
      write (6,*)
      write (6,*)'  -z:   This switch is to define the zeroangle for the
     # first cell of the first detector, as '
      write (6,*)'        obtained by a nickel calibration. The default 
     #value is 0.'
      write (6,*)
      write (6,*)'  -l:   This switch is to define the wavelength of inc
     #ident neutrons (in A), as obtained by a'
      write (6,*)'        nickel calibration. The default value is 0, fo
     #r which the output spectrum has angles'
      write (6,*)'        as abcise. A positive value produces spectra i
     #n Q-scale using the given wavelength.'
      write (6,*)
      write (6,*)'  -n:   This switch is to define the normalisation con
     #stant and method. The argument must be '
      write (6,*)'        an integer. If the number is the positive, it 
     #will be used for monitor normalisation.'
      write (6,*)'        If it is negative, -time is the time (in sec)
     #to which data will be normalised.'
      write (6,*)
      write (6,*)'  -r:   Directory where data are located. If data are
     #in the same directory as program, the '
      write (6,*)'        argument should be -r ./'
      write (6,*)
      write (6,*)'  -g:   Getall file containing switches defining the r
     #equired info as function of numors. '
      write (6,*)'        See for example the input file d4cget.get.'
      write (6,*)
      write (6,*)' In case of problems using this program, please visit
     #D4 web page (www.ill.fr/YelloBook/D4/home.html)'
      write (6,*)' or contact the author (Gabriel Cuello, cuello@ill.fr)
     #.'
      write (6,*)
      write (6,*)'-----------------------------------------------------'
      
c----6---1---------2---------3---------4---------5---------6---------7--
      stop
      return
      end
c----6---1---------2---------3---------4---------5---------6---------7--
c----6---1---------2---------3---------4---------5---------6---------7--
c----6---1---------2---------3---------4---------5---------6---------7--

