      SUBROUTINE SUMROM(AAA,AREA,PSUM,NRAD,A2,R1,R2,MS)
	include 'dimension.inc'
	include 'beam.inc'
      COMMON /cylabs/ NAN,RAD1(mcont),rad2(mcont),THETA,azi,PI
     1,AMUS(mcont),AMUT(mcont),DEN(mcont),abstemp(mcont)
      DIMENSION PATH(3)
      REAL LIS(5),LSS(5),LIST,LISN,LSST,LSSN
c
c sine of detector azimuthal angle
c
	detfac=abs(sin(theta))
      A3=A2
      OMGADD=0.
      IF(A3.LT.0.) OMGADD=PI
      A3=ABS(A3)
      IF(A3.GE.R2) ALIM=R2
      IF(A3.LT.R2) ALIM=A
      IF(A3.LE.R1) ALIM=R1
      M1=MS*((ALIM-R1)/(R2-R1)+0.0001)
      M2=MS-M1
      AAA=0.
      AREA=0.
      PSUM=0.
C
C INTEGRATE THE BEAM PROFILE FROM 0 TO ALIMIT
C
      WIDTH=A3
      IF(A3.GE.R2) WIDTH=R2
      ASTEP=WIDTH/40
      P1=PROBE(0.,PROFIL,NPROF,PRSTEP,A,B)/2
      DO 10 I=1,40
      X=I*ASTEP
      X=SIGN(X,A2)
      P2=PROBE(X,PROFIL,NPROF,PRSTEP,A,B)/2
      PSUM=PSUM+P1+P2
      P1=P2
   10 CONTINUE
      PSUM=PSUM*ASTEP
      OAD=PI-azi
      IF(M1.EQ.0) GO TO 41
      N=M1
      RSTEP=(ALIM-R1)/M1
      RADD=-0.5*RSTEP+R1
   31 DO 30 M=1,N
      R=M*RSTEP+RADD
      NOMEG=PI*R/RSTEP
      OMEGST=PI/NOMEG
      OMEGAD=-0.5*OMEGST+OMGADD
      AREAY=R*RSTEP*OMEGST*AMUS(NRAD)
      ICOUNT=0
      SUM1=0.
      SUM2=0.
      ARSUM=0.
      I=1
   35 IF(I.GT.NOMEG) GO TO 36
      OMEGA=I*OMEGST+OMEGAD
      D=R*SIN(OMEGA)
      IF(ABS(D).GT.A3) GO TO 101
C
C DETERMINE A PROFILE VALUE FOR THIS 'D' VALUE
C
      PROB=PROBE(D,PROFIL,NPROF,PRSTEP,A,B)
C
C CALCULATE DISTANCE INCIDENT NEUTRON PASSES THROUGH EACH ANNULUS
C
      DO 60 J=1,NAN
      LIST=DIST(R,RAD1(J),OMEGA)
      LISN=DIST(R,RAD2(J),OMEGA)
        LIS(J)=LISN-LIST
  60  CONTINUE
C
C CALCULATE DISTANCE SCATTERED NEUTRON PASSES THROUGH EACH ANNULUS
C
      O=OMEGA+OAD
c
c check this element is visible to the detector
c (Apparently this was not done in previous versions of the programme!)
c
	D=R*sin(O)
	if(D.lt.A1.and.D.gt.B1) then
		DO 63 J=1,NAN
      			LSST=DIST(R,RAD1(J),O)
      			LSSN=DIST(R,RAD2(J),O)
        		LSS(J)=LSSN-LSST
  63  		CONTINUE
      		ANGLE=OMEGA*180/PI
C
C CALCULATE ABSORPTION FOR PATH THROUGH ALL ANNULI
C
      		PATH(1)=0
      		DO 65 II=1,NAN
      			PATH(1)=PATH(1)+AMUT(II)*(LIS(II)+LSS(II))
  65  		CONTINUE
			if(detfac.gt.1.0e-10) then
	      		SUM1=SUM1+EXP(-PATH(1)/detfac)*PROB
			endif
      		ARSUM=ARSUM+PROB
	endif
   37 I=I+1
      GO TO 35
  101 I=NOMEG-I+2
      GO TO 35
   36 AAA=AAA+SUM1*AREAY
      AREA=AREA+ARSUM*AREAY
   30 CONTINUE
   41 IF(M2.EQ.0) RETURN
      N=M2
      RSTEP=(R2-ALIM)/M2
      RADD=-0.5*RSTEP+ALIM
      M2=0
      GO TO 31
      END
