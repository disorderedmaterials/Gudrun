      SUBROUTINE ACYL(MS1,AREAS,AREAC)
	include 'dimension.inc'
	include 'beam.inc'
      COMMON /cylabs/ NAN,RAD1(mcont),rad2(mcont),THETA,azi,PI
     1,AMUS(mcont),AMUT(mcont),DEN(mcont),abstemp(mcont)
C
C     DETERMINE MAXIMUM OF THE PROFILE FUNCTION
C
      PSUM=0.
      DO 10 I=1,NPROF
      PSUM=AMAX1(PSUM,PROFIL(I))
      MS=MS1
  10  CONTINUE
      AREAS=0
      ASS=0
      DO 15 I=1,NAN
      abstemp(i)=1.0
      IF(amus(I).EQ.0.)GO TO 15
C
C  NO STEPS ARE CHOSEN SO THAT STEP WIDTH IS THE SAME FOR ALL ANNULI
C
      MS=MS1*(RAD2(i)-RAD1(i))/(RAD2(1)-RAD1(1))
	if(ms.lt.1) ms=1 
      CALL SUMROM(AAAA,AREAA,PSUM,I,A,RAD1(I),RAD2(I),MS)
      CALL SUMROM(AAAB,AREAB,PSUM,I,B,RAD1(I),RAD2(I),MS)
 
      AREAS1=SIGN(AREAA,A)-SIGN(AREAB,B)
      if(areas1.gt.0.) abstemp(i)=(SIGN(AAAA,A)-SIGN(AAAB,B))/areas1
  15  CONTINUE
      RETURN
      END


