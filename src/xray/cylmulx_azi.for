	subroutine CYLMULTOF(nant,radt1,radt2,theight,dent,sigtl,captcs)
c
c	Original   from AKS
c
C PROGRAM MULSCA - CONTROL PROGRAM FOR THE MULTIPLE SCATTERING
C CALCULATION 
C               
c Modified March 2001 into a subroutine to be linked into GUDRUN
c
c	modified 5/02/2003 to include situation with azimuthal detectors
c
	include 'dimension.inc'
	include 'beam.inc'
	include 'mul_corr_azi.inc'
	integer*4 NT2(10000)
      real*4 SUMA(mcorrang,mcorrang,mcont,mcont)
     *,SUMB(mcorrang,mcorrang,mcont,mcont)
     *,SUMC(mcorrang,mcorrang,mcont)
     *,SUMP(mcorrang,mcorrang,mcont),ALEN1(2000),NTEST(2000)
     *,ALEN2(mcorrang,10000),OM(2000)
     *,SUMD(mcorrang,mcorrang,mcont,mcont)
     *,SUME(mcorrang,mcorrang,mcont,mcont)
     *,SUMS(mcorrang,mcorrang,mcont)
     *,SUMT(mcorrang,mcorrang,mcont)
     *,SIGSL(mcorrwav,mcont),SIGTL(mcorrwav,mcont)
     *,SIGT(mcont),MUS
	real*4 dent(mcont)		!temporary density values
	real*4 radt1(mcont),radt2(mcont)!temporary radius values
	real*4 captcs(mcont)		!capture cross sections
      COMMON/ANNBLO/ NAN,RAD1(mcont),rad2(mcont),NRAD(mcont),RST(mcont)
     *,DEN(mcont),SIGS(mcont),MUS(mcont)
      COMMON/INTBLO/ R(20),RSQ(20),NOMEGA(20),AREAEL(20)
      COMMON/CHRIS/CSUMPR(mcorrwav,mcorrang,mcorrang)
     *,CSUMB(mcorrwav,mcorrang,mcorrang)
     *,CTOSCA(mcorrwav,mcorrang,mcorrang)
     *,CTOSC(mcorrwav,mcorrang,mcorrang)
C                                                                     
C ARRAYS - PRIMARY SCATTERING:-                                       
C        SUMC = C-C                                                     
C        SUMS = C-P                                                    
C        SUMT = P-C                                                     
C        SUMP = P-P                                                    
C                                                                      
C  ARRAYS - SECONDARY SCATTERING:-                                      
C                                                                      
C        SUMA = C-C                                                     
C        SUMD = C-P                                                    
C        SUME = P-C                                                    
C        SUMB = P-P                                                    
C                                                                      
c output unit for diagnostic messages
c
	nw=6
	write(nw,*) 'Got to CylMultof'
      NBLEN1=2000                                                      
      NBLEN2=10000
C                                                                      
C SET NZ IN THE RATIO OF HEIGHT/ASTEP                                  
C                                                                       
	astep=stepm
	height=theight
c
c set beam height parameters with reference to bottom of sample rather
c than centre: this is to correspond to definitions used in cylmultof
c
	h2=0.5*height
	hdown=hbdown+h2
	hup=hbup+h2
	hsdown=hsbdown+h2
	hsup=hsbup+h2
      NZ=(HEIGHT+0.5*ASTEP)/ASTEP                                      
      IF(NZ.GT.200) NZ=200                                              
      IF(NZ.LT.1) NZ=1                                                  MUL01890
C                                                                       MUL01900
C DEFINE INTEGERS CORRESPONDING AS CLOSE AS POSSIBLE                    MUL01910
C TO THE HEIGHT OF THE INCIDENT BEAM                                    MUL01920
C                                                                       MUL01930
	IF(HUP.GT.height) HUP=height
      IF(HDOWN.LT.0.0) HDOWN=0.0
      HIGHB=HUP-HDOWN                                                   MUL01960
      ZSTEP=HEIGHT/NZ                                                   MUL01970
      MIN=(HDOWN+0.5*ZSTEP)/ZSTEP+1                                     MUL01980
      MAX=1+(HUP-0.5*ZSTEP)/ZSTEP                                       MUL01990
C                                                                       MUL02000
C                                                                       MUL02010
C DEFINE INTEGERS CORRESPONDING AS CLOSE AS POSSIBLE TO THE             MUL02020
C HEIGHT OF THE SCATTERED BEAM                                          MUL02030
C                                                                       MUL02040
      IF(HSUP.GT.HEIGHT) HSUP = HEIGHT
      IF(HSDOWN.LT.0.) HSDOWN=0.                                        MUL02060
      HIGHS=HSUP-HSDOWN                                                 MUL02070
      MINS=(HSDOWN+0.5*ZSTEP)/ZSTEP+1                                   MUL02080
      MAXS=(HSUP-0.5*ZSTEP)/ZSTEP+1                                     MUL02090
C                                                                       MUL02100
C DETERMINE HEIGHT OF OVERLAP BETWEEN INCIDENT AND SCATTERED BEAMS      MUL02110
C                                                                       MUL02120
      UP=AMIN1(HUP,HSUP)                                                MUL02130
      DOWN=AMAX1(HDOWN,HSDOWN)                                          MUL02140
      HIGHBS=UP-DOWN                                                    MUL02150
      IF(HIGHBS.LT.0.) HIGHBS=0.                                        MUL02160
c	write(6,*) a,b,a1,b1
c	write(6,*) HIGHB,HIGHS,HIGHBS,HEIGHT
C
C CALCULATE  SCATTERING C/S AT EACH WAVELENGTH ASSUMING 1/LAMBDA ABSORPTION
C
	nan=nant
	do i=1,nan
		den(i)=dent(i)
		rad1(i)=radt1(i)
		rad2(i)=radt2(i)
!	write(6,*) i,den(i),rad1(i),rad2(i)
	end do
      DO IL=1,nwavmul
		do i=1,nan
!			FAC=wavmul(IL)/1.7979
!			SIGSL(IL,I)=SIGTL(Il,I)-FAC*captCS(I)
			SIGSL(IL,I)=SIGTL(Il,I)-captCS(I)
      			SIGS(I)=SIGSL(IL,I)
      			SIGT(I)=SIGTL(IL,I)
      			MUS(I)=SIGT(I)*DEN(I)
!	write(6,*) i,sigs(i),sigt(i),mus(i)
		end do
C
C SET UP THE GEOMETRY OF THE PROBLEM AND ARRAYS ALEN1,ALEN2             MUL00480
C                                                                       MUL00490
!		write(6,*) 'Going to LEN3'
      	CALL LEN3(nw,MS,ASTEP,OM,NBLEN1)	
      	MS=NRAD(1)
c		write(6,*) 'Going to LEN4'
      	CALL LEN4(nw,OM,ALEN1,NTEST,NBLEN1)
!		write(6,*) 'Going to LEN5'
      	CALL LEN5(nw,OM,ALEN2,NT2,NBLEN1,NBLEN2)
!            write(6,*) (alen2(1,j),j=1,20)
C                                                                       MUL00530
C PERFORM THE PRIMARY AND MULTIPLE SCATTERING CALCULATION               MUL00540
C                                                                       MUL00550
!		write(6,*) 'Going to SUMMU1'
      	CALL SUMMU1(MS,NZ,MIN,MAX,MINS,MAXS,SUMA,SUMD,SUME,SUMB
     1,SUMC,SUMS,SUMT,SUMP,OM,ALEN1,ALEN2,NTEST,NT2,NBLEN1,NBLEN2)
!	write(6,*) (SUMP(1,kk,1),kk=1,nazimul)
!	do i=1,nan
!            do j=1,nan
!      write(6,*) j,i
!      write(6,*) suma(il,1,j,i),sumd(il,1,j,i),sume(il,1,j,i)
!     1,sumb(il,1,j,i)
!            end do
!      end do
C                                                                       MUL00580
C OUTPUT THE RESULTS                                                    MUL00590
C                                                                       MUL00600
!		write(6,*) 'Going to DATAOP'
!            write(6,*) height,highb,highs,highbs
      	CALL DATAOP(SUMC,SUMS,SUMT,SUMP,SUMA,SUMD,SUME,SUMB
     1,HIGHB,HIGHS,HIGHBS,IL)
	end do
!      write(6,*) nan,csumpr(1,1,1),ctosca(1,1,1)
	do k3=1,nazimul
         DO K2=1,nangmul
            DO K1=1,nwavmul
      	   onescat(k1,k2,k3)=csumpr(k1,k2,k3)
	         mulscat(k1,k2,k3)=ctosca(k1,k2,k3)
	      end do
	   end do
	end do
!      write(6,*) nan,onescat(1,1,1),mulscat(1,1,1)
 888  FORMAT(1X,I3,F10.5)
 777  FORMAT(1X,F10.4,3F11.8)
      return
	end
