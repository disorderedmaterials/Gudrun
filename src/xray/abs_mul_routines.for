      SUBROUTINE ACYL(MS1,AREAS,AREAC)
	include 'dimension.inc'
	include 'beam.inc'
      COMMON /cylabs/ NAN,RAD1(mcont),rad2(mcont),THETA,azi,PI
     1,AMUS(mcont),AMUT(mcont),DEN(mcont),abstemp(mcont)
     1,AMUT_OUT(mcont)
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

      SUBROUTINE SUMROM(AAA,AREA,PSUM,NRAD,A2,R1,R2,MS)
	include 'dimension.inc'
	include 'beam.inc'
      COMMON /cylabs/ NAN,RAD1(mcont),rad2(mcont),THETA,azi,PI
     1,AMUS(mcont),AMUT(mcont),DEN(mcont),abstemp(mcont)
     1,AMUT_OUT(mcont)
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
		do J=1,NAN
      		LSST=DIST(R,RAD1(J),O)
      		LSSN=DIST(R,RAD2(J),O)
            LSS(J)=LSSN-LSST
            if(detfac.gt.1.0e-5) then
                lss(j)=lss(j)/detfac
            end if
      	end do
      	ANGLE=OMEGA*180/PI
C
C CALCULATE ABSORPTION FOR PATH THROUGH ALL ANNULI
C
      	PATH(1)=0
      	do II=1,NAN
      		PATH(1)=PATH(1)+(AMUT(II)*LIS(II)+AMUT_OUT(II)*LSS(II))
            end do
    		if(path(1).lt.30.0) then
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

      FUNCTION DIST(R1,RAD,OMEGA)
      R=R1
      DIST=0.
      B=R*SIN(OMEGA)
      IF(ABS(B).GT.RAD) RETURN
      T=R*COS(OMEGA)
      C=RAD*RAD-B*B
      D=SQRT(C)
      IF(R.GT.RAD) GO TO 1
      DIST=T+D
      RETURN
    1 DIST=D*(1+SIGN(1.,T))
      RETURN
      END
C                                                                       MUL05400
C        ********************                                           MUL05410
C        *                  *                                           MUL05420
C        *  FUNCTION DIST2  *                                           MUL05430
C        *                  *                                           MUL05440
C        ********************                                           MUL05450
C                                                                       MUL05460
      FUNCTION DIST2(R1,RAD,OMEGA)                                      MUL05470
C                                                                       MUL05480
C CALCULTES THE MODULUS OF THE DISTANCE OF THE POINT (R,OMEGA)          MUL05490
C FROM THE RIGHTMOST SURFACE OF A CIRCLE OF RADIUS RAD ALONG A LINE     MUL05500
C PARALLEL TO OMEGA=0 AXIS.  IF R.GT. RAD THE VALUE RETURNED IS         MUL05510
C THE LENGTH OF THE SEGMENT WHICH THE LINE CUTS.  IF                    MUL05520
C  ABS(R*SIN(OMEGA)) IS .GT. RAD, THERE IS NO PATH THROUGH THE          MUL05530
C CIRCLE AND A VALUE OF ZERO IS RETURNED                                MUL05540
C                                                                       MUL05550
      DIST2=0                                                           MUL05560
      R=R1                                                              MUL05570
      D=R*SIN(OMEGA)                                                    MUL05580
      IF(ABS(D).GT.RAD) RETURN                                          MUL05590
      T=R*COS(OMEGA)                                                    MUL05600
      DIST2=RAD*RAD-D*D                                                 MUL05610
      DIST2=SQRT(DIST2)                                                 MUL05620
      IF(R.GT.RAD) GO TO 1                                              MUL05630
      DIST2=DIST2-T                                                     MUL05640
      RETURN                                                            MUL05650
    1 DIST2=DIST2-SIGN(DIST2,T)                                         MUL05660
      RETURN                                                            MUL05670
      END                                                               MUL05680
C                                                                       MUL05690
C      ********************                                             MUL05700
C      *                  *                                             MUL05710
C      *  FUNCTION DIST3  *                                             MUL05720
C      *                  *                                             MUL05730
C      ********************                                             MUL05740
C                                                                       MUL05750
      FUNCTION DIST3(R1,RAD,OMEGA)                                      MUL05760
C                                                                       MUL05770
C SAME COMMENTS AS FOR DIST2 EXCEPT THAT THE DISTANCE CALCULATED        MUL05780
C IS FROM THE POINT (R,OMEG) TO THE LEFT-MOST SURFACE                   MUL05790
C                                                                       MUL05800
      DIST3=0                                                           MUL05810
      R=R1                                                              MUL05820
      D=R*SIN(OMEGA)                                                    MUL05830
      IF(ABS(D).GT.RAD) RETURN                                          MUL05840
      T=R*COS(OMEGA)                                                    MUL05850
      DIST3=RAD*RAD-D*D                                                 MUL05860
      DIST3=SQRT(DIST3)                                                 MUL05870
      IF(R.GT.RAD) GO TO 1                                              MUL05880
      DIST3=DIST3+T                                                     MUL05890
      RETURN                                                            MUL05900
    1 DIST3=DIST3+SIGN(DIST3,T)                                         MUL05910
      RETURN                                                            MUL05920
      END                                                               MUL05930
C                                                                       MUL05940
C     *******************                                               MUL05950
C     *                 *                                               MUL05960
C     *  FUNCTION ARSIN *                                               MUL05970
C     *                 *                                               MUL05980
C     *******************                                               MUL05990
C                                                                       MUL06000
      FUNCTION ARSIN(ASINE,SIGCOS)                                      MUL06010
C                                                                       MUL06020
CRETURNS THE ARCSINE OF ASINE BETWEEN PI AND -PI.  SIGCOS               MUL06030
C IS THE SIGN OF COS(ARCSINE(ASINE))                                    MUL06040
C                                                                       MUL06050
      ACOS=1-ASINE*ASINE                                                MUL06060
      ACOS=SQRT(ACOS)*SIGCOS                                            MUL06070
      ARSIN=ATAN2(ASINE,ACOS)                                           MUL06080
      RETURN                                                            MUL06090
      END                                                               MUL06100
C                                                                       MUL14050
C       *******************                                             MUL14060
C       *                 *                                             MUL14070
C       *  FUNCTION DIST1 *                                             MUL14080
C       *                 *                                             MUL14090
C       *******************                                             MUL14100
C                                                                       MUL14110
      FUNCTION DIST1(R1,R2,RAD,B,ASIG)                                  MUL14120
      DIST1=0                                                           MUL14130
      IF(ABS(B).GT.RAD) RETURN                                          MUL14140
      RA=R1                                                             MUL14150
      IF(RA.GT.RAD) RA=RAD                                              MUL14160
      X=RA*RA-B*B                                                       MUL14170
	if(x.lt.0.) x=0.
      DIST1=DIST1+SQRT(X)                                               MUL14180
      RA=R2                                                             MUL14190
      IF(RA.GT.RAD) RA=RAD                                              MUL14200
      X=RA*RA-B*B                                                       MUL14210
	if(x.lt.0.) x=0.
      DIST1=DIST1+SQRT(X)*ASIG                                          MUL14220
      DIST1=ABS(DIST1)                                                  MUL14230
      RETURN                                                            MUL14240
      END                                                               MUL14250
C                                                                       MUL14260
C       *******************                                             MUL14270
C       *                 *                                             MUL14280
C       *  FUNCTION AVL2  *                                             MUL14290
C       *                 *                                             MUL14300
C       *******************                                             MUL14310
C                                                                       MUL14320
      FUNCTION AVL2(AA,BB,L0)                                           MUL14330
      REAL L0,AA,BB
      REAL*8 L1,L2,L2R,L2R2,L2R3,L2R4,L2R5                              MUL14350
      DOUBLE PRECISION A,B,AL3,AL2,AR,AR2,AR3,AR4,AR5,BR,BR2,BR3,       MUL14360
     1AL2R,AL2R2,AL2R3,AL2R4,AL2R5,TEMP,CONS,                           MUL14370
     2Z1,Z2,Z3,Z4,A1,A2,BX,THRD,C1,C2,RAT,RAT1,RAT2,RAT3,RAT4,          MUL14380
     3E1,E2,COEFF1,COEFF2,COEFF3,DIFF1,DIFF2,SUM,SUM1,SUM2              MUL14390
      SUM=0.                                                            MUL14400
      THRD=1./3.                                                        MUL14410
C                                                                       MUL14420
C CHECK THAT A IS .GE. B - IF NOT WE SWITCH VALUES                      MUL14430
C                                                                       MUL14440
      A=AA                                                              MUL14450
      B=BB                                                              MUL14460
      IF(A.GE.B) GO TO 5                                                MUL14470
      A=BB                                                              MUL14480
      B=AA                                                              MUL14490
    5 CONTINUE                                                          MUL14500
C                                                                       MUL14510
C SELECT LIMITS FOR INTEGRATION                                         MUL14520
C                                                                       MUL14530
      L1=L0+B                                                           MUL14540
      L2=L0-B                                                           MUL14550
      AL3=0                                                             MUL14560
      IF(L2.LE.0.) AL3=DABS(L2)                                         MUL14570
C                                                                       MUL14580
C WE REDEFINE THE UPPER LIMIT AS L2 AND THE LOWER LIMIT AS AL2.         MUL14590
C                                                                       MUL14600
      AL2=DABS(L2)                                                      MUL14610
      L2=L1                                                             MUL14620
      IF(L0.LE.0.) GO TO 20                                             MUL14630
      AR=A/L0                                                           MUL14640
      BR=B/L0                                                           MUL14650
      IF(BR.GT.0.01) GO TO 888                                          MUL14660
      AVL2=1                                                            MUL14670
      RETURN                                                            MUL14680
  888 CONTINUE                                                          MUL14690
      L2R=L2/L0                                                         MUL14700
      AL2R=AL2/L0                                                       MUL14710
      AR2=AR*AR                                                         MUL14720
      AR3=AR2*AR                                                        MUL14730
      AR4=AR3*AR                                                        MUL14740
      AR5=AR4*AR                                                        MUL14750
      BR2=BR*BR                                                         MUL14760
      BR3=BR2*BR                                                        MUL14770
      L2R2=L2R*L2R                                                      MUL14780
      L2R3=L2R2*L2R                                                     MUL14790
      L2R4=L2R3*L2R                                                     MUL14800
      L2R5=L2R4*L2R                                                     MUL14810
      AL2R2=AL2R*AL2R                                                   MUL14820
      AL2R3=AL2R2*AL2R                                                  MUL14830
      AL2R4=AL2R3*AL2R                                                  MUL14840
      AL2R5=AL2R4*AL2R                                                  MUL14850
      TEMP=BR2-1                                                        MUL14860
      CONS=9./(8.*BR3)                                                  MUL14870
      CONS=CONS/AR3                                                     MUL14880
      Z1=AR2*TEMP                                                       MUL14890
      Z2=THRD*(-TEMP-AR2)                                               MUL14900
      Z3=0.5*AR4                                                        MUL14910
      A1=0.5*(L2R*Z1+AR2*L2R2+Z2*L2R3-0.5*L2R4+0.2*L2R5-Z3)             MUL14920
      A2=0.5*(AL2R*Z1+AR2*AL2R2+Z2*AL2R3-0.5*AL2R4+0.2*AL2R5-Z3)        MUL14930
      BX=-THRD*(-AR3*TEMP+0.2*AR5)                                      MUL14940
      Z1=-THRD*(-TEMP+0.2*AR2)*AR                                       MUL14950
      Z2=0.5*AR                                                         MUL14960
      Z3=-0.2*AR                                                        MUL14970
      Z4=0.5*AR3                                                        MUL14980
      C1=L2R*Z4+L2R2*Z1+L2R3*Z2+L2R4*Z3                                 MUL14990
      C2=AL2R*Z4+AL2R2*Z1+AL2R3*Z2+AL2R4*Z3                             MUL15000
      RAT1=L2/A                                                         MUL15010
      RAT2=AL2/A                                                        MUL15020
      RAT3=1+RAT1                                                       MUL15030
      RAT4=1+RAT2                                                       MUL15040
      E1=A1+BX                                                          MUL15050
      E2=A2+BX                                                          MUL15060
      COEFF1=E1*DLOG(RAT3)-E2*DLOG(RAT4)                                MUL15070
      RAT3=1-RAT1                                                       MUL15080
      RAT4=1-RAT2                                                       MUL15090
      RAT3=DABS(RAT3)                                                   MUL15100
      RAT4=DABS(RAT4)                                                   MUL15110
      DIFF1=BX-A1                                                       MUL15120
      DIFF2=BX-A2                                                       MUL15130
      IF(DABS(DIFF2).LT.1.0E-05) RAT4=1.                                MUL15140
      IF(DABS(DIFF1).LT.1.0E-05) RAT3=1.                                MUL15150
      COEFF2=DIFF1*DLOG(RAT3)-DIFF2*DLOG(RAT4)                          MUL15160
      COEFF3=C1-C2                                                      MUL15170
      SUM1=CONS*(COEFF1+COEFF2+COEFF3)                                  MUL15180
      SUM=SUM+SUM1                                                      MUL15190
C                                                                       MUL15200
C COMPUTE INTEGRAL FOR ENCLOSED REGION - IF APLICABLE                   MUL15210
C                                                                       MUL15220
   20 IF(AL3.LE.0.) GO TO 10                                            MUL15230
      RAT=AL3/A                                                         MUL15240
      IF(RAT.GE.1.) GO TO 30                                            MUL15250
      RAT2=RAT*RAT                                                      MUL15260
      RAT3=RAT2*RAT                                                     MUL15270
      DIFF1=1.-RAT2                                                     MUL15280
      DIFF2=(1+RAT)/(1-RAT)                                             MUL15290
      A1=-0.5*DIFF1**2*DLOG(DIFF2)                                      MUL15300
      SUM2=(RAT3+RAT+A1)*9.*L0**2*A/(8*B**3)                            MUL15310
      GO TO 40                                                          MUL15320
   30 SUM2=2.25                                                         MUL15330
   40 SUM=SUM+SUM2                                                      MUL15340
   10 AVL2=SUM                                                          MUL15350
      RETURN                                                            MUL15360
      END                                                               MUL15370
C                                                                       MUL15380
C       **********************                                          MUL15390
C       *                    *                                          MUL15400
C       *  SUBROUTINE EQUIV  *                                          MUL15410
C       *                    *                                          MUL15420
C       **********************                                          MUL15430
C                                                                       MUL15440
      SUBROUTINE EQUIV(R1,R2,MS1,AR,ARB,ARS,ARBS)
	include 'dimension.inc'
	include 'beam.inc'
	include 'mul_corr_azi.inc'
      DIMENSION ARBS(mcorrang)
      PI=3.141592653
	piconv=pi/180.0
      PI2=2*PI                                                          MUL15490
      MS=(1.-R1/R2)*MS1                                                 MUL15500
      IF(MS.LT.1) MS=1
      RSTEP=(R2-R1)/MS                                                  MUL15510
      RADD=R1-0.5*RSTEP                                                 MUL15520
      AR=0.                                                             MUL15530
      ARB=0.                                                            MUL15540
      ARS=0.                                                            MUL15550
      DO 5 I=1,nazimul
      ARBS(I)=0.
   5  continue                                                          MUL15570
      DO 10 IR=1,MS                                                     MUL15580
      R=RADD+IR*RSTEP                                                   MUL15590
      NOM=PI2*R/RSTEP                                                   MUL15600
      OMST=PI2/NOM                                                      MUL15610
      AREL=RSTEP*R*OMST                                                 MUL15620
      DO 20 IO=1,NOM                                                    MUL15630
      OMEG=(IO-1)*OMST                                                  MUL15640
      D=R*SIN(OMEG)                                                     MUL15650
      P1=PROBE(D,PROFIL,NPROF,PRSTEP,A,B)                               MUL15660
      P2=1.                                                             MUL15670
      IF(D.GT.A1.OR.D.LT.B1) P2=0.                                      MUL15680
      AR=AR+AREL                                                        MUL15690
      ARB=ARB+P1*AREL                                                   MUL15700
      ARS=ARS+P2*AREL                                                   MUL15710
      DO 30 I=1,nazimul
	ang=azimul(i)*piconv
      D1=R*SIN(OMEG-ANG)
      P2=1.                                                             MUL15740
      IF(D1.GT.A1.OR.D1.LT.B1) P2=0.                                    MUL15750
      ARBS(I)=ARBS(I)+P1*P2*AREL                                        MUL15760
   30 CONTINUE                                                          MUL15770
   20 CONTINUE                                                          MUL15780
   10 CONTINUE                                                          MUL15790
      RETURN                                                            MUL15800
      END                                                               MUL15810

C                                                                       MUL02620
C        **********************                                         MUL02630
C        *                    *                                         MUL02640
C        *   SUBROUTINE LEN3  *                                         MUL02650
C        *                    *                                         MUL02660
C        **********************                                         MUL02670
C                                                                       MUL02680
      SUBROUTINE LEN3(nw,MS,ASTEP,OM,NBLEN1)
	include 'dimension.inc'
      DIMENSION OM(NBLEN1)                                              MUL02700
      COMMON/INTBLO/ R1(20),RSQ(20),NOMEGA(20),AREAEL(20)               MUL02710
      COMMON/ANNBLO/ NAN,RAD1(mcont),rad2(mcont),NRAD(mcont),RST(mcont)
     *,DEN(mcont),SIGS(mcont),MUS(mcont)
      REAL MUS                                                          MUL02730
C                                                                       MUL02740
C THIS SUB. SETS UP ARRAY TO BE USED IN LEN4,LEN5 AND SUMMU1            MUL02750
C NAN IS THE THE NO. OF ANNULI.
C MS IS THE NUMBER OF STEPS FROM rad1 TO RAD2 - IT IS USED TO
C DETERMINE THE SIZE OF THE STEP IN EACH ANNULUS                        MUL02780
C THE ITH ANNULUS IS DEFINED TO LIE BETWEEEN RAD1(I) AND RAD2(I)        MUL02790
C WITHE A DENSITY OF DEN(MOLS PER CC), SCATTERING C/S SIGS(I)           MUL02800
C BARNS AND TOTAL ATTENUATION C/S SIGT(I) BARNS,AND TOTAL               MUL02810
C ABSORPTION COEFFICIENT MUS(I) PER CM.  ON EXIT, NRAD(I) IS THE        MUL02820
C NUMBER OF ELEMENTAL RINGS TO BE INTEGRATED OVER,                      MUL02830
C NOMEGA(I) IS THE NUMBER OMEGA VALUES IN THE ITH RING, AND OM(K)       MUL02840
C IS THE KTH OMEGA VALUE. RSQ(I) IS R1(I)**2 WHERE R1(I)                MUL02850
C IS THE RADIUS OF THE ITH RING                                         MUL02860
C                                                                       MUL02870
C                                                                       MUL02880
C SET A NUMBER TO CHECK THE NUMBER OF ELEMENTS IN R1 DOES NOT EXCEED 20 MUL02890
C                                                                       MUL02900
      NKRAD=20                                                          MUL02910
      PI=3.141592653                                                    MUL02920
      NOM=0                                                             MUL02930
      KRAD=0                                                            MUL02940
C                                                                       MUL02950
      NOMSTA=1                                                          MUL02960
      DO 10 I=1,NAN                                                     MUL02970
C                                                                       MUL02980
C DO NOT WISH TO INTEGRATE OVER NON-SCATTERING ANNULI                   MUL02990
C                                                                       MUL03000
      IF(den(I).EQ.0) GOTO 15
C                                                                       MUL03020
C DEFINE SCATTERING C/S PER UNIT VOLUME OF ITH ANNULUS                  MUL03030
C                                                                       MUL03040
      FAC=DEN(I)*SIGS(I)/(4*PI)                                         MUL03050
C                                                                       MUL03060
C DEFINE NO. OF RADIUS STEPS IN THE ITH ANNULUS                         MUL03070
C                                                                       MUL03080
      Y=RAD1(I)
      Z=RAD2(i)-Y
      X=Z/ASTEP+0.5                                                     MUL03120
      NS=INT(X)                                                         MUL03130
      IF(NS.LT.1) NS=1                                                  MUL03140
      NRAD(I)=NS                                                        MUL03150
      RSTEP=Z/NS                                                        MUL03160
      RST(I)=RSTEP                                                      MUL03170
      RADD=-0.5*RSTEP+RAD1(I)
C                                                                       MUL03190
C SET UP NS RADIUS VALUES                                               MUL03200
C                                                                       MUL03210
      DO 20 J=1,NS                                                      MUL03220
      R=J*RSTEP+RADD                                                    MUL03230
      KRAD=KRAD+1                                                       MUL03240
C                                                                       MUL03250
C CHECK THAT KRAD IS NOT GREATER THAN NKRAD                             MUL03260
C                                                                       MUL03270
      IF(KRAD.GT.NKRAD) GO TO 50                                        MUL03280
      R1(KRAD)=R                                                        MUL03290
      RSQ(KRAD)=R**2                                                    MUL03300
C                                                                       MUL03310
C DEFINE NO. OF STEPS IN OMEGA FOR EACH RADIUS VALUE                    MUL03320
C                                                                       MUL03330
C THE NUMBER OF ELEMENTS INANY GIVEN RING IS SET TO BE AN INTEGER       MUL03340
C MULTIPLE OF THE NUMBER IN THE PREVIOUS RING                           MUL03350
C                                                                       MUL03360
      X=PI*R/ASTEP                                                      MUL03370
      ATEST=X/FLOAT(NOMSTA)+0.5                                         MUL03380
      NTEST=INT(ATEST)                                                  MUL03390
   70 NOMEG=NTEST*NOMSTA                                                MUL03400
      IF(NOMEG.LT.65) GO TO 80                                          MUL03410
      NTEST=NTEST-1                                                     MUL03420
      GO TO 70                                                          MUL03430
   80 NOMSTA=NOMEG                                                      MUL03440
      OMEGST=PI/NOMEG                                                   MUL03450
C                                                                       MUL03460
C DEFINE THE AREA ELEMENTS AT THIS RADIUS AND MULTIPLY BY               MUL03470
C THE APPROPRIATE SCATTERING C/S                                        MUL03480
C                                                                       MUL03490
      AREAEL(KRAD)=R*OMEGST*RSTEP*FAC                                   MUL03500
      NOMEG=2*NOMEG                                                     MUL03510
      NOMEGA(KRAD)=NOMEG                                                MUL03520
      DO 30 L=1,NOMEG                                                   MUL03530
      NOM=NOM+1                                                         MUL03540
C                                                                       MUL03550
C CHECK THAT NOM IS NOT GREATER THAN NBLEN1                             MUL03560
C                                                                       MUL03570
      IF(NOM.GT.NBLEN1) GO TO 51                                        MUL03580
      OM(NOM)=(L-1)*OMEGST                                              MUL03590
   30 CONTINUE                                                          MUL03600
!      write(6,*) i,nrad(i),krad,nomega(krad)
   20 CONTINUE                                                          MUL03610
      GO TO 10                                                          MUL03620
   15 NRAD(I)=0                                                         MUL03630
   10 CONTINUE                                                          MUL03640
      RETURN                                                            MUL03650
   50 WRITE(NW,60) NKRAD                                                MUL03660
      STOP                                                              MUL03670
   51 WRITE(NW,61) NBLEN1                                               MUL03680
      STOP                                                              MUL03690
   60 FORMAT(1X,29HDIMENSION SPECIFIER KRAD .GT.,I5,13H IN SUB. LEN3)   MUL03700
   61 FORMAT(1X,28HDIMENSION SPECIFIER NOM .GT.,I5,13H IN SUB. LEN3)    MUL03710
      END                                                               MUL03720
C                                                                       MUL03730
C        *********************                                          MUL03740
C        *                   *                                          MUL03750
C        *  SUBROUTINE LEN4  *                                          MUL03760
C        *                   *                                          MUL03770
C        *********************                                          MUL03780
C                                                                       MUL03790
      SUBROUTINE LEN4(nw,OM,ALEN,NTEST,NBLEN)
	include 'dimension.inc'
	include 'beam.inc'
      real*4 ALEN(NBLEN),OM(NBLEN),NTEST(NBLEN),MUS
      COMMON/ANNBLO/ NAN,RAD1(mcont),rad2(mcont),NRAD(mcont),RST(mcont)
     *,DEN(mcont),SIGS(mcont),MUS(mcont)
      COMMON/INTBLO/ R1(20),RSQ(20),NOMEGA(20),AREAEL(20)               MUL03830
C                                                                       MUL03860
C THIS SETS UP THE ARRAY ALEN CORRESPONDING TO THE EXPONENTIAL          MUL03870
C SUM EXP(-SUMMATION(MUS*L1(R,OMEGA))).  THE SUMMATION IS               MUL03880
C OVER THE ANNULI AND L1 IS THE DISTANCE OF THE                         MUL03890
C CURRENT ANNULUS THROUGH THAT ANNULUS IN A DIRECTION PARALLEL          MUL03900
C TO INCIDENT BEAM.                                                     MUL03910
C ALSO AN ARRAY NTEST IS SETUP WHICH IS 1 IF A GIVEN ELEMNT IS          MUL03920
C IN INCIDENT BEAM AND ZERO IF NOT                                      MUL03930
C                                                                       MUL03940
C CALL DIST3                                                            MUL03950
C                                                                       MUL03960
      NOM=0                                                             MUL03970
      NRA=0                                                             MUL03980
C                                                                       MUL03990
C STEP THROUGH ANNULUS VALUES                                           MUL04000
C                                                                       MUL04010
      DO 10 I=1,NAN                                                     MUL04020
      NR=NRAD(I)                                                        MUL04030
      IF(NR.EQ.0) GO TO 10                                              MUL04040
      DO 20 J=1,NR                                                      MUL04050
      NRA=NRA+1                                                         MUL04060
      R=R1(NRA)                                                         MUL04070
      NOMEG=NOMEGA(NRA)                                                 MUL04080
C                                                                       MUL04090
C STEP THROUGH OMEGA VALUES                                             MUL04100
C                                                                       MUL04110
      DO 30 K=1,NOMEG                                                   MUL04120
      NOM=NOM+1                                                         MUL04130
      OM1=OM(NOM)                                                       MUL04140
C                                                                       MUL04150
C TEST WHETHER CURRENT ELEMENT IS IN BEAM                               MUL04160
C                                                                       MUL04170
      D=R*SIN(OM1)                                                      MUL04180
      NTEST(NOM)=PROBE(D,PROFIL,NPROF,PRSTEP,A,B)                       MUL04190
C                                                                       MUL04200
C AT EACH R1,OM1 VALUE L1 IS CALCULATED THROUGH EACH ANNULUS ALONG      MUL04210
C A LINE PARALLEL TO OM0=0, IN THE DIREDTION LEFT TO RIGHT              MUL04220
C                                                                       MUL04230
      SUM=0                                                             MUL04240
      DO 40 L=1,NAN                                                     MUL04260
C                                                                       MUL04280
C CALCULATE THE DISTANCE THROUGH A CIRCLE OF RADIUS RAD(L1)             MUL04290
C                                                                       MUL04300
      XSTART=DIST3(R,RAD1(L),OM1)
      XNEXT=DIST3(R,RAD2(L),OM1)                                        MUL04310
C                                                                       MUL04320
C XSTART IS THE VALUE THAT DISTANCE THROUGH THE NEXT SMALLEST           MUL04330
C ANNULUS- HENCE IF WE SUBTRACT XSTART FORM XNEXT, THIS WILL LEAVE      MUL04340
C THE DISTANCE THROUGH THE CURRENT ANNULUS                              MUL04350
C                                                                       MUL04360
      SUM=SUM+MUS(L)*(XNEXT-XSTART)                                     MUL04370
      XSTART=XNEXT                                                      MUL04380
   40 CONTINUE                                                          MUL04390
C                                                                       MUL04400
C BEFORE OUTPUT WE MUST RECONSTRUCT SCATTERING FUNCTIONS FOR EACH       MUL04410
C ANNULUS USING THE ACTUAL AREAS FOR THE GIVEN GEOMETRY                 MUL04420
C                                                                       MUL04430
      ALEN(NOM)=EXP(-SUM)                                               MUL04440
   30 CONTINUE                                                          MUL04450
   20 CONTINUE                                                          MUL04460
   10 CONTINUE                                                          MUL04470
      RETURN                                                            MUL04480
      END                                                               MUL04490
C                                                                       MUL04500
C         *********************                                         MUL04510
C         *                   *                                         MUL04520
C         *  SUBROUTINE LEN5  *                                         MUL04530
C         *                   *                                         MUL04540
C         *********************                                         MUL04550
C                                                                       MUL04560
      SUBROUTINE LEN5(NW,OM,ALEN,NTEST,NBLEN1,NBLEN2)
	include 'dimension.inc'
	include 'beam.inc'
	include 'mul_corr_azi.inc'
      real*4 OM(NBLEN1),ALEN(mcorrang,NBLEN2)
	integer*4 NTEST(NBLEN2)
      COMMON/ANNBLO/ NAN,RAD1(mcont),rad2(mcont),NRAD(mcont),RST(mcont)
     *,DEN(mcont),SIGS(mcont),MUS(mcont)
      COMMON/INTBLO/R1(20),RSQ(20),NOMEGA(20),AREAEL(20)                MUL04610
      REAL*4 MUS
C                                                                       MUL04630
C THIS SETS UP THE ARRAY ALEN2 CORRESPONDING TO THE VALUES              MUL04640
C OF EXP(-SUMMATION(MUS*L2(R,OMEGA))), WHERE L2(R,OMEGA) IS THE         MUL04650
C DISTANCE THROUGH A GIVEN ANNULUS FROM THE SCATTERING POINT            MUL04660
C (R,OMEGA) PARALLEL TO THE THE DIRECTION OF THE                        MUL04670
C COUNTER.  THE SUMMATION IS OVER ALL THE ANNULI, AND THERE IS          MUL04680
C AN ARRAY ELEMENT FOR EACH SCATTERING POINT WITHIN THE SAMPLE          MUL04690
C                                                                       MUL04700
C CALLS DIST2                                                           MUL04710
      PI=3.141592653                                                    MUL04720
	piconv=pi/180.0
      NALEN=0                                                           MUL04730
      NOM=0                                                             MUL04740
      NRA=0                                                             MUL04750
C                                                                       MUL04760
C STEP THROUGH THE ANNULI AND RADIUS VALUES                             MUL04770
C                                                                       MUL04780
      DO 10 I=1,NAN                                                     MUL04790
      NR=NRAD(I)                                                        MUL04800
      IF(NR.EQ.0) GO TO 10                                              MUL04810
      DO 20 J=1,NR                                                      MUL04820
      NRA=NRA+1                                                         MUL04830
      R=R1(NRA)                                                         MUL04840
      NOMEG=NOMEGA(NRA)                                                 MUL04850
C                                                                       MUL04860
C STEP THROUGH THE OMEGA VALUES                                         MUL04870
C                                                                       MUL04880
      NSAVE=NOM                                                         MUL04890
      DO 50 iazi=1,nazimul
      NOM=NSAVE                                                         MUL04910
      ANGLE=azimul(iazi)*piconv
      DO 30 K=1,NOMEG                                                   MUL04930
      NOM=NOM+1                                                         MUL04940
      OM1=OM(NOM)-ANGLE                                                 MUL04950
C                                                                       MUL04960
C DECIDE WHETHER SCATTERING POINT IS WITHIN THE SCATTERED BEAM -        MUL04970
C IF SO WE SET NTEST TO UNITY, OTHERWISE IT IS ZERO.                    MUL04980
C                                                                       MUL04990
      NALEN=NALEN+1                                                     MUL05000
C                                                                       MUL05010
C CHECK NALEN IS NOT GREATER THAN NBLEN2                                MUL05020
C                                                                       MUL05030
      IF(NALEN.GT.NBLEN2) GO TO 60                                      MUL05040
      NTEST(NALEN)=1                                                    MUL05050
      D=R*SIN(OM1)                                                      MUL05060
      IF(D.GT.A1.OR.D.LT.B1) NTEST(NALEN)=0                             MUL05070
C                                                                       MUL05080
C AT EACH R,OM1 VALUE THE DISTANCE L2 IS CALCULATED THROUGH EACH        MUL05090
C ANNULUS, ALONG A LINE PARALLEL TO OMEGA=0, IN THE DIRECTION           MUL05100
C LEFT TO RIGHT                                                         MUL05110
C                                                                       MUL05120
      SUM=0.                                                            MUL05130
      DO 40 L=1,NAN                                                     MUL05150
C                                                                       MUL05170
C CALCULATE THE DISTANCE THROUGH A CIRCLE OF RADIUS RAD(L)              MUL05180
C                                                                       MUL05190
      Xstart=DIST2(R,RAD1(L),OM1)                                       MUL05200
      XNEXT=DIST2(R,RAD2(L),OM1)                                        MUL05200
C                                                                       MUL05210
C XSTART IS THE VALUE OF THAT DISTANCE THROUGH THE NEXT                 MUL05220
C SMALLEST ANNULUS. HENCE IF WE SUBTRACT XSTART FROM XNEXT,             MUL05230
C THIS WILL LEAVE THE DISTANCE THROUGH THE CURRENT ANNULUS              MUL05240
C                                                                       MUL05250
      SUM=SUM+MUS(L)*(XNEXT-XSTART)                                     MUL05260
      XSTART=XNEXT                                                      MUL05270
   40 CONTINUE                                                          MUL05280
c
c step through the detector azimuthal angles (corresponds to angmul in
c cylindrical sample coordinates)
c
	do iang=1,nangmul
	   detfac=abs(sin(angmul(iang)*piconv))
c
c perform azimuthal correction
c
	   if(detfac.gt.1.0e-10) then
	      ALEN(iang,NALEN)=EXP(-SUM/detfac)
	   else
		alen(iang,nalen)=0.0
	   endif
	end do
   41 CONTINUE                                                          MUL05300
   30 CONTINUE                                                          MUL05310
   50 CONTINUE                                                          MUL05320
   20 CONTINUE                                                          MUL05330
   10 CONTINUE                                                          MUL05340
      RETURN                                                            MUL05350
   60 WRITE(NW,61) NBLEN2                                               MUL05360
      STOP                                                              MUL05370
   61 FORMAT(1X,30HDIMENSION SPECIFIER NALEN .GT.,I5,13H IN SUB. LEN5)  MUL05380
      END                                                               MUL05390
C                                                                       MUL07450
C       ***********************                                         MUL07460
C       *                     *                                         MUL07470
C       *  SUBROUTINE SUMMU1  *                                         MUL07480
C       *                     *                                         MUL07490
C       ***********************                                         MUL07500
C                                                                       MUL07510
      SUBROUTINE SUMMU1(MS,NZ,MIN,MAX,MINS,MAXS,SUMA,SUMD,SUME,SUMB     MUL07520
     1,SUMC,SUMS,SUMT,SUMP,OM,ALEN1,ALEN2,NTEST,NT2,NBLEN1,NBLEN2)
	include 'dimension.inc'
	include 'beam.inc'
	include 'mul_corr_azi.inc'
      COMMON/ANNBLO/ NAN,RAD1(mcont),rad2(mcont),NRAD(mcont),RST(mcont)
     *,DEN(mcont),SIGS(mcont),MUS(mcont)
      COMMON/INTBLO/ R(20),RSQ(20),NOMEGA(20),AREAEL(20)                MUL07560
	external cosd,sind
	integer*4 NT2(NBLEN2)
      real*4 SUMHA(mcorrang,mcorrang,2)
     *,SUMHB(mcorrang,mcorrang,2),ZVALSQ(200)
     *,SUMA(mcorrang,mcorrang,mcont,mcont)
     *,SUMB(mcorrang,mcorrang,mcont,mcont),STORSQ(128),STORSM(128)
     2,SUMP(mcorrang,mcorrang,mcont),SUMC(mcorrang,mcorrang,mcont)
     *,PRODP(mcorrang,mcorrang),PROD(mcorrang,mcorrang)
     *,PRODS(mcorrang,mcorrang),PRODT(mcorrang,mcorrang)
     3,OM(NBLEN1),ALEN1(NBLEN1),NTEST(NBLEN1),ALEN2(mcorrang,NBLEN2)
     4,ZERO(20),EQRAD(20),AFACT(200),BFACT(200),ZMULL(200)              MUL07620
     5,AREAA(mcorrang,2),AREAB(mcorrang,2),AREAC(mcorrang)
     *,AREAP(mcorrang)
     6,SUMD(mcorrang,mcorrang,mcont,mcont)
     *,SUME(mcorrang,mcorrang,mcont,mcont)
     *,SUMS(mcorrang,mcorrang,mcont)
     *,SUMT(mcorrang,mcorrang,mcont)
     7,AREAS(mcorrang),AREAT(mcorrang),AREAD(mcorrang,2)
     *,AREAE(mcorrang,2),DFACT(200)
     8,EFACT(200),SUMHE(mcorrang,mcorrang,2),SUMHD(mcorrang,mcorrang,2)
 	real*4 tema1(mcorrang),tema2(mcorrang)
 	real*4 temb1(mcorrang),temb2(mcorrang)
 	real*4 temd1(mcorrang),temd2(mcorrang)
 	real*4 teme1(mcorrang),teme2(mcorrang)
      REAL MUS
      PI=3.141592653                                                    MUL07680
C                                                                       MUL07690
C THIS ROUTINE CALCULATES THE PRIMARY AND SECONDARY SCATTERING          MUL07700
C FOR A SERIES OF CONCENTRIC INCOHERENT SCATTERERS OF ANNULAR           MUL07710
C CROSS-SECTION WITH A VARIETY OF THICKNESSES, USING EQ (13) OF         MUL07720
C BLECH AND AVERBACH 1965 AS A STARTING POINT.  THE ELEMENTAL           MUL07730
C AREAS ARE DEFINED IN SUB LEN3.  THE GEOMETRY HAS BEEN DEFINED         MUL07740
C SO AS TO KEEP THE NUMBER OF COMPUTATIONS OF  INTER-                   MUL07750
C ELEMENTAL DISTANCES (AND ANGLES) TO A MINIMUM --SEE NOTES FOR         MUL07760
C A FULLER EXPLANATION.                                                 MUL07770
C   THE INTEGRAL IS CALCULATED FOR A SERIES OF Z' VALUES AND AT         MUL07780
C THE END THE SUMS ARE WEIGHTED BY THE NUMBER OF TIMES THEY APPEAR      MUL07790
C IN THE FINAL INTEGRAL AND ARE ADDED TOGETHER                          MUL07800
C  MOST OF THE VARIABLE NAMES HAVE EITHER BEEN INTRODUCED IN            MUL07810
C THE CALL ROUTINE OR THEIR USE SHOULD BE DEDUCTABLE FROM THE           MUL07820
C PROGRAM COMMENTS                                                      MUL07830
C  RHE BEAM LIES BETWEEN X=A AND X=B, WHERE X=0 IS AN AXIS PARALLEL     MUL07840
C TO THE DIRECTION OF THE BEAM THROUGH CENTRE OF THE SYSTEM.  THE       MUL07850
C INTEGRAL IS CALCULATED BOTH FOR WHEN THE BEAM IS PARTIALLY            MUL07860
C EMMERSED AND FOR WHEN COMPLETELY EMMERSED                             MUL07870
C                                                                       MUL07880
C   THE VALUES OF ALEN1 AND ALEN2 ARE READ FROM BINARY FILES            MUL07890
C                                                                       MUL07900
C INITIALIZE VARIOUS ARRAYS                                             MUL07910
C                                                                       MUL07920
      ALIMIT=1.0E-30                                                    MUL07930
      ZSTEP=HEIGHT/NZ                                                   MUL07940
      ZADD=-ZSTEP                                                       MUL07950
      ZSTEPS=ZSTEP*ZSTEP                                                MUL07960
C                                                                       MUL07970
C CALCULATE HEIGHT OF BEAM FOR PARTIAL EMMERSION                        MUL07980
C                                                                       MUL07990
      HIGHB=(MAX-MIN+1)*ZSTEP                                           MUL08000
      HIGHS=(MAXS-MINS+1)*ZSTEP                                         MUL08010
      HIGHCC=HEIGHT*HEIGHT                                              MUL08020
      HIGHCP=HEIGHT*HIGHS                                               MUL08030
      HIGHPC=HIGHB*HEIGHT                                               MUL08040
      HIGHPP=HIGHB*HIGHS                                                MUL08050
c
c zero accumulators
c
	do k=1,nan
	   do j=1,nan
	      do iazi=1,nazimul
 	         do iang=1,nangmul
	            SUMA(IANG,iazi,J,K)=0.0
                  SUMD(IANG,iazi,J,K)=0.0
                  SUME(IANG,iazi,J,K)=0.0
                  SUMB(IANG,iazi,J,K)=0.0
	         end do
	      end do
	   end do
	end do
    2 CONTINUE                                                          MUL08150
      DO 1 I=1,NZ                                                       MUL08160
      Z=I*ZSTEP+ZADD                                                    MUL08170
      ZVALSQ(I)=Z*Z                                                     MUL08180
      AFACT(I)=0.                                                       MUL08190
      DFACT(I)=0.                                                       MUL08200
      EFACT(I)=0.                                                       MUL08210
      BFACT(I)=0.                                                       MUL08220
    1 CONTINUE                                                          MUL08230
C                                                                       MUL08240
C CALCULATE SUMMATION FACTORS FOR COMPLETE AND PARTIAL BEAMS            MUL08250
C                                                                       MUL08260
      DO 410 IZ1=1,NZ                                                   MUL08270
      DO 411 IZ2=1,NZ                                                   MUL08280
C                                                                       MUL08290
C COMPLETE BEAM, COMPLETE SCATTERING - C-C                              MUL08300
C                                                                       MUL08310
      IZP=IZ2-IZ1                                                       MUL08320
      IZP=1+IABS(IZP)                                                   MUL08330
      AFACT(IZP)=AFACT(IZP)+1.                                          MUL08340
C                                                                       MUL08350
C COMPLETE BEAM, PARTIAL SCATTERING - C-P                               MUL08360
C                                                                       MUL08370
      IF(IZ2.LT.MINS.OR.IZ2.GT.MAXS) GO TO 412                          MUL08380
      DFACT(IZP)=DFACT(IZP)+1.                                          MUL08390
412   CONTINUE                                                          MUL08400
C                                                                       MUL08410
C PARTIAL BEAM, COMPLETE SCATTERING - P-C                               MUL08420
C                                                                       MUL08430
      IF(IZ1.LT.MIN.OR.IZ1.GT.MAX) GO TO 411                            MUL08440
      EFACT(IZP)=EFACT(IZP)+1.                                          MUL08450
C                                                                       MUL08460
C PARTIAL BEAM, PARTIAL SCATTERING - P-P                                MUL08470
C                                                                       MUL08480
      IF(IZ2.LT.MINS.OR.IZ2.GT.MAXS) GO TO 411                          MUL08490
      BFACT(IZP)=BFACT(IZP)+1.                                          MUL08500
  411 CONTINUE                                                          MUL08510
  410 CONTINUE
!      write(6,*) min,max,mins,maxs
!      write(6,*) afact(2),dfact(2),efact(2),bfact(2)
C                                                                       MUL08530
C SET UP MEAN VALUES OF 1/(L**2) FOR SINGLE ELEMENTS                    MUL08540
C                                                                       MUL08550
      NR1=1                                                             MUL08560
      DO 400 I1=1,NAN                                                   MUL08570
      RSTEP=RST(I1)                                                     MUL08580
      NRAD1=NRAD(I1)                                                    MUL08590
      IF(NRAD1.EQ.0) GO TO 400                                          MUL08600
      DO 401 J1=1,NRAD1                                                 MUL08610
      TSTEP=2*PI*R(NR1)/NOMEGA(NR1)                                     MUL08620
C                                                                       MUL08630
C GENERATE AN EQUIVALENT RADIUS BASED ON THE SURFACE AREA               MUL08640
C OF A CUBE, AND EVALUATE THE ZEROTH CONTRIBUTION                       MUL08650
C                                                                       MUL08660
      AREA=2*(RSTEP*TSTEP+RSTEP*ZSTEP+TSTEP*ZSTEP)                      MUL08670
      ZERO(NR1)=33.78/AREA                                              MUL08680
      ARAD=2.25/ZERO(NR1)                                               MUL08690
      EQRAD(NR1)=SQRT(ARAD)                                             MUL08700
      NR1=NR1+1                                                         MUL08710
  401 CONTINUE                                                          MUL08720
  400 CONTINUE                                                          MUL08730
C                                                                       MUL08740
C***********************************************************************MUL08750
C                                                                       MUL08760
C CALCULATE PRIMARY SCATTERING FOR COMPLETE AND PARTIALLY EMMERSED      MUL08770
C BEAMS                                                                 MUL08780
C                                                                       MUL08790
      NOM1=1                                                            MUL08800
      NR1=0                                                             MUL08810
      NSERCH=1                                                          MUL08820
      DO 6 I1=1,NAN                                                     MUL08830
C                                                                       MUL08840
C INITIALIZE THE AREA SUMS WHICH WILL BE USED TO NORMALIZE SUMC,        MUL08850
C SUMS, SUMT,  AND SUMB                                                 MUL08860
C                                                                       MUL08870
      DO IAZI=1,nazimul
         AREAC(IAZI)=0.
         AREAS(IAZI)=0.
         AREAT(IAZI)=0.
         AREAP(IAZI)=0.
 	   do iang=1,nangmul
	      SUMC(IANG,iazi,I1)=0.0
            SUMS(IANG,iazi,I1)=0.0
            SUMT(IANG,iazi,I1)=0.0
            SUMP(IANG,iazi,I1)=0.0
	   end do
	end do
      NRAD1=NRAD(I1)                                                    MUL08980
      IF(NRAD1.EQ.0) GOTO 6                                             MUL08990
      DO 7 J1=1,NRAD1                                                   MUL09000
      NR1=NR1+1                                                         MUL09010
      NOMEG1=NOMEGA(NR1)                                                MUL09020
      AREA1=AREAEL(NR1)                                                 MUL09030
      DO 302 IAZI=1,nazimul
	do iang=1,nangmul
         PROD(IANG,iazi)=0.
         PRODS(IANG,iazi)=0.
         PRODT(IANG,iazi)=0.
         PRODP(IANG,iazi)=0.
	end do
  302 CONTINUE                                                          MUL09090
C                                                                       MUL09100
C PICK UP THE APPROPRIATE VALUES OF ALEN1 AND ALEN2                     MUL09110
C                                                                       MUL09120
      DO 5 IAZI=1,nazimul
      NDOM1=NOM1                                                        MUL09140
      DO 8 K1=1,NOMEG1                                                  MUL09150
      RT1T=NTEST(NDOM1)                                                 MUL09160
      NT2T=NT2(NSERCH)                                                  MUL09170
      AREAC(IAZI)=AREAC(IAZI)+AREA1                                     MUL09180
      AREAS(IAZI)=AREAS(IAZI)+AREA1*NT2T                                MUL09190
      AREAT(IAZI)=AREAT(IAZI)+AREA1*RT1T                                MUL09200
      AREAP(IAZI)=AREAP(IAZI)+AREA1*RT1T*NT2T                           MUL09210
	do iang=1,nangmul
	   XPROD=ALEN1(NDOM1)*ALEN2(iang,NSERCH)
         PROD(IANG,iazi)=PROD(IANG,iazi)+XPROD
         PRODS(IANG,iazi)=PRODS(IANG,iazi)+XPROD*NT2T
         PRODT(IANG,iazi)=PRODT(IANG,iazi)+XPROD*RT1T
         PRODP(IANG,iazi)=PRODP(IANG,iazi)+XPROD*RT1T*NT2T
      end do
	NDOM1=NDOM1+1
      NSERCH=NSERCH+1                                                   MUL09280
    8 CONTINUE                                                          MUL09290
    5 CONTINUE                                                          MUL09300
c	write(6,*) ((PRODP(iang,kk),kk=1,nazimul),iang=1,nangmul)
      NOM1=NOM1+NOMEG1                                                  MUL09310
      DO IAZI=1,nazimul
         do iang=1,nangmul
	      PROD(IANG,iazi)=PROD(IANG,iazi)*AREA1
            PRODS(IANG,iazi)=PRODS(IANG,iazi)*AREA1
            PRODT(IANG,iazi)=PRODT(IANG,iazi)*AREA1
            PRODP(IANG,iazi)=PRODP(IANG,iazi)*AREA1
            SUMC(IANG,iazi,I1)=SUMC(IANG,iazi,I1)+PROD(IANG,iazi)
            SUMS(IANG,iazi,I1)=SUMS(IANG,iazi,I1)+PRODS(IANG,iazi)
            SUMT(IANG,iazi,I1)=SUMT(IANG,iazi,I1)+PRODT(IANG,iazi)
            SUMP(IANG,iazi,I1)=SUMP(IANG,iazi,I1)+PRODP(IANG,iazi)
	   end do
	end do
    7 CONTINUE
      do iang=1,nangmul
	   DO IAZI=1,nazimul
            SUMC(IANG,iazi,I1)=SUMC(IANG,iazi,I1)/AREAC(IAZI)
            IF(AREAS(IAZI).GT.ALIMIT)
     1SUMS(IANG,iazi,I1)=SUMS(IANG,iazi,I1)/AREAS(IAZI)
            IF(AREAT(IAZI).GT.ALIMIT)
     1SUMT(IANG,iazi,I1)=SUMT(IANG,iazi,I1)/AREAT(IAZI)
            IF(AREAP(IAZI).GT.ALIMIT)
     1SUMP(IANG,iazi,I1)=SUMP(IANG,iazi,I1)/AREAP(IAZI)
         end do
	end do
    6 CONTINUE                                                          MUL09540
C                                                                       MUL09550
C***********************************************************************MUL09560
C                                                                       MUL09570
C  CALCULATE SECONDARY SCATTERING                                       MUL09580
C                                                                       MUL09590
C DEFINE R1 VALUES                                                      MUL09600
C                                                                       MUL09610
      NR1=0                                                             MUL09620
      NOM1=1                                                            MUL09630
      DO 10 I1=1,NAN                                                    MUL09640
      NRAD1=NRAD(I1)                                                    MUL09650
C                                                                       MUL09660
C AVOID INTEGRATING OVER ANNULI WHICH DO NOT SCATTER                    MUL09670
C                                                                       MUL09680
      IF(NRAD1.EQ.0) GO TO 10                                           MUL09690
      NR1S=NR1                                                          MUL09700
      NOM1S=NOM1                                                        MUL09710
      NR2=NR1                                                           MUL09720
      NOM2=NOM1                                                         MUL09730
      DO 30 I2=I1,NAN                                                   MUL09740
      NRAD2=NRAD(I2)                                                    MUL09750
C                                                                       MUL09760
C DO NOT BOTHER WITH ELEMENTS WHICH DO NOT SCATTER                      MUL09770
C                                                                       MUL09780
      IF(NRAD2.EQ.0) GO TO 30                                           MUL09790
C                                                                       MUL09800
C SAVE VALUE OF NR2,NOM2                                                MUL09810
C                                                                       MUL09820
      NR2S=NR2                                                          MUL09830
      NOM2S=NOM2                                                        MUL09840
C                                                                       MUL09850
C RESET NR1                                                             MUL09860
C                                                                       MUL09870
      NR1=NR1S                                                          MUL09880
      NOM1=NOM1S                                                        MUL09890
C                                                                       MUL09900
C SET AREA PRODUCTS AND SUMH'S TO ZERO                                  MUL09910
C                                                                       MUL09920
      DO 11 J=1,2                                                       MUL09930
      DO 12 IAZI=1,nazimul
      AREAA(IAZI,J)=0.                                                  MUL09950
      AREAD(IAZI,J)=0.                                                  MUL09960
      AREAE(IAZI,J)=0.                                                  MUL09970
      AREAB(IAZI,J)=0.                                                  MUL09980
	do iang=1,nangmul
         SUMHA(IANG,iazi,J)=0.
         SUMHD(IANG,iazi,J)=0.
         SUMHE(IANG,iazi,J)=0.
         SUMHB(IANG,iazi,J)=0.
	end do
   12 CONTINUE                                                          MUL10030
   11 CONTINUE                                                          MUL10040
      DO 20 J1=1,NRAD1                                                  MUL10050
      NR1=NR1+1                                                         MUL10060
      R1=R(NR1)                                                         MUL10070
      R1SQ=RSQ(NR1)                                                     MUL10080
      R12=2*R1                                                          MUL10090
      NOMEG1=NOMEGA(NR1)                                                MUL10100
      AREA1=AREAEL(NR1)                                                 MUL10110
C                                                                       MUL10120
C SET A TEMPORARY VALUE OF NOM1                                         MUL10130
C                                                                       MUL10140
      NTEMP1=NOM1                                                       MUL10150
C                                                                       MUL10160
C DEFINE R2 VALUES  THESE AR E ALWAYS THE SAME OR GREATER THEN          MUL10170
C R1 VALUES                                                             MUL10180
C                                                                       MUL10190
      NR2=NR2S                                                          MUL10200
      NOM2=NOM2S                                                        MUL10210
C                                                                       MUL10220
C IF 2ND RING IS IN SAME ANNULUS AS FIRST THEN WE HAVE ENSURED          MUL10230
C THAT THE INITIAL R2 VAUE IS NOT LESS THAN R1 VALUE.  THIS             MUL10240
C WILL ALSO REDUCE THE NUMBER OF 2ND RING S TO BE SUMMED OVER           MUL10250
C                                                                       MUL10260
      J2STAR=1                                                          MUL10270
      IF(I1.NE.I2) GO TO 39                                             MUL10280
      NR2=NR1                                                           MUL10290
      NOM2=NOM1                                                         MUL10300
      J2STAR=J1                                                         MUL10310
   39 DO 40 J2=J2STAR,NRAD2                                             MUL10320
      R2=R(NR2)                                                         MUL10330
      R2SQ=RSQ(NR2)                                                     MUL10340
      ADD=R1SQ+R2SQ                                                     MUL10350
      NOMEG2=NOMEGA(NR2)                                                MUL10360
C                                                                       MUL10370
C TAKE RATIO OF NOMEG2 TO NOMEG1 - IF THE PROBLEM WAS SET UP CORRECTLY  MUL10380
C IN SUB LEN3, THIS SHOULD ALWAYS BE AN INTEGER.                        MUL10390
C                                                                       MUL10400
      NRAT=NOMEG2/NOMEG1                                                MUL10410
      AREA2=AREAEL(NR2)                                                 MUL10420
      ARPROD=AREA1*AREA2                                                MUL10430
C                                                                       MUL10440
C DEFINE OMEGA VALUES FOR FIRST ELEMENT                                 MUL10450
C                                                                       MUL10460
      OM1=OM(NTEMP1)                                                    MUL10470
C                                                                       MUL10480
C STEP THROUGH ELEMENTS IN 2ND RING AND CLACULATE HORIZONTAL            MUL10490
C SEPARATION OF ELEMENTS AND ABSORPTION PATH LENGTH.  SINCE THE         MUL10500
C THE GEOMETRY IS SYMMETRIC ABOUT OMEGA=0, THIS NEED ONLY BE DONE       MUL10510
C FOR HALF THE RING                                                     MUL10520
C                                                                       MUL10530
      NLAST=NOMEG2/2                                                    MUL10540
      NLAST1=NLAST+1                                                    MUL10550
C                                                                       MUL10560
C FOR THE CASES OMDIF=0 AND OMDIF=PI, THE CALCULATION OF ADDA AND       MUL10570
C SUML DOES NOT INVOLVE CALCULATION OF SINES AND COSINES                MUL10580
C                                                                       MUL10590
      RADD=R1+R2                                                        MUL10600
      RDIFF=R2-R1                                                       MUL10610
      STORSQ(1)=RDIFF*RDIFF                                             MUL10620
      STORSQ(NLAST1)=RADD*RADD                                          MUL10630
      SUMR1=0                                                           MUL10640
      SUMR2=0                                                           MUL10650
      RADST1=R1-RAD1(1)
      RADST2=R2-RAD1(1)
C                                                                       MUL10680
C THE METHOD IS TO FORM TH SUMMATION OF PRODUCTS                        MUL10690
C MUS(IAN)*L FROM THE CENTRE TO R1 (SUMR1) AND FROM THE CENTRE          MUL10700
C TO R2 (SUMR2).  FOR OMDIF=0 SUML=SUMR2-SUMR1, AND FOR OMDIF=PI        MUL10710
C SUML=SUMR1+SUMR2                                                      MUL10720
C                                                                       MUL10730
      DO 59 IAN=1,NAN                                                   MUL10740
      AMU=MUS(IAN)                                                      MUL10750
      SUMR1=SUMR1+RADST1*AMU                                            MUL10770
      RADNX1=R1-RAD2(IAN)
      IF(RADNX1.LT.0) GO TO 58                                          MUL10790
      SUMR1=SUMR1-RADNX1*AMU                                            MUL10800
      RADST1=RADNX1                                                     MUL10810
   59 CONTINUE                                                          MUL10820
   58 CONTINUE
      DO 56 IAN=1,NAN
      AMU=MUS(IAN)                                                      MUL10840
      SUMR2=SUMR2+RADST2*AMU                                            MUL10860
      RADNX2=R2-RAD2(IAN)
      IF(RADNX2.LT.0) GO TO 57                                          MUL10880
      SUMR2=SUMR2-RADNX2*AMU                                            MUL10890
      RADST2=RADNX2                                                     MUL10900
   56 CONTINUE                                                          MUL10910
C                                                                       MUL10920
C NORMALIZE THE SUMS TO THE PATH LENGTH IN EACH CASE                    MUL10930
C                                                                       MUL10940
   57 CONTINUE
      STORSM(1)=MUS(I1)                                                 MUL10950
      IF(NR1.EQ.NR2) GO TO 54                                           MUL10960
      STORSM(1)=(SUMR2-SUMR1)/RDIFF                                     MUL10970
   54 STORSM(NLAST1)=(SUMR1+SUMR2)/RADD                                 MUL10980
C                                                                       MUL10990
C SET A TEMPORARY VALUE OF NOM2. THIS IS ONE GREATER THAN THE           MUL11000
C CURRENT VALUE OF NOM2, WHICH CORRESPONDS TO OMDIF=0 AND WHICH HAS     MUL11010
C ALREADY BEEN SUMMED FOR.                                              MUL11020
C (EXCUSE ENGLISH_)                                                     MUL11030
C                                                                       MUL11040
      NTEMP2=NOM2+1                                                     MUL11050
      DO 60 K2=2,NLAST                                                  MUL11060
      OM2=OM(NTEMP2)                                                    MUL11070
      OMDIF=OM2-OM1                                                     MUL11080
C
C CALCULATE HORIZONTAL SEPARTION OF ELEMENTS                            MUL11090
C                                                                       MUL11100
      OMCOS=R2*COS(OMDIF)                                               MUL11110
      ADDA=ADD-R12*OMCOS                                                MUL11120
      XADD=SQRT(ADDA)                                                   MUL11130
C                                                                       MUL11140
C CALCULATE THE PERPENDICULAR DISTANCE OF THE LINE JOINING              MUL11150
C (R2,OMDIF) AND (R1,0)                                                 MUL11160
C FROM THE CNETRE                                                       MUL11170
C                                                                       MUL11180
      SUML=0                                                            MUL11190
      D=R1*R2*SIN(OMDIF)/XADD                                           MUL11200
      E=R1-OMCOS                                                        MUL11210
      SIGCOS=SIGN(1.,E)                                                 MUL11220
C                                                                       MUL11230
C CLACULTE THE SEPARATION OF THE POINTS - IF EITHER OR BOTH             MUL11240
C LIE OUTSIDE A CIRCLE OF RADIUS  RAD(J), THEN THE LENGHT               MUL11250
C OF THEIR SEPARAITON WITHIN THAT CIRCLE OS CALCULATED                  MUL11260
C                                                                       MUL11270
      DO 70 IAN=1,NAN                                                   MUL11290
      J=IAN+1                                                           MUL11300
C                                                                       MUL11310
C DEFINE LENGTH FO PATH BETWEEN THE TWO POINTS WHICH LIES WITHIN        MUL11320
C THE CIRCLE FO RAIUS RAD(J)                                            MUL11330
C                                                                       MUL11340
      XSTART=DIST1(R1,R2,RAD1(ian),D,SIGCOS)
      XNEXT=DIST1(R1,R2,RAD2(ian),D,SIGCOS)
C                                                                       MUL11360
C THE PATH LENGTH THROUGH ITH ANNULUS IS NOW XNEXT-XSTART  --           MUL11370
C XSTART WAS THE PATH THROUGH THE NEXT SMAALER CIRCLE                   MUL11380
C                                                                       MUL11390
      SUML=SUML+MUS(IAN)*(XNEXT-XSTART)                                 MUL11400
   70 CONTINUE                                                          MUL11420
C                                                                       MUL11430
C CONVERT SUML TO A FRACTIONAL PATH LENGTH OF XADD                      MUL11440
C                                                                       MUL11450
      SUML=SUML/XADD                                                    MUL11460
      STORSQ(K2)=ADDA                                                   MUL11470
      STORSM(K2)=SUML                                                   MUL11480
      NTEMP2=NTEMP2+1                                                   MUL11490
   60 CONTINUE                                                          MUL11500
C                                                                       MUL11510
C GENERATE REMAINING NOMEG2-NLAST ELEMENTS                              MUL11520
C                                                                       MUL11530
      NSAVE=NLAST1                                                      MUL11540
      NLAST=NLAST1+1                                                    MUL11550
      DO 80 K2=NLAST,NOMEG2                                             MUL11560
      NSAVE=NSAVE-1                                                     MUL11570
      STORSQ(K2)=STORSQ(NSAVE)                                          MUL11580
      STORSM(K2)=STORSM(NSAVE)                                          MUL11590
   80 CONTINUE                                                          MUL11600
C                                                                       MUL11610
C NOW CALCULATE  THE PRODUCTS ALEN1*ALEN2 FOR ALL                       MUL11620
C COMBINATIONS OF THE ELEMENTS IN 1ST AND 2ND RINGS - THESE ARE         MUL11630
C STORED IN NOMEG2 TEMPORARY SUMS TEMPSU (PARTIAL) AND                  MUL11640
C TEMCSU (COMPLETE), WHICH MUST INITIALLY BE SET TO ZERO                MUL11650
C                                                                       MUL11660
C                                                                       MUL11670
C RESET NTEMP2                                                          MUL11680
C                                                                       MUL11690
      NTEMP2=NOM2                                                       MUL11700
C                                                                       MUL11710
C STEP THROUGH NOMEG2 ELEMENTS FOR SECOND RING                          MUL11720
C                                                                       MUL11730
      DO 140 K2=1,NOMEG2                                                MUL11740
      STORE=STORSQ(K2)                                                  MUL11750
      STORM=STORSM(K2)                                                  MUL11760
C                                                                       MUL11770
C GENERATE THE Z DEPENDENT ARRAY ZMULL, WHICH INDEPENDENT OF            MUL11780
C K1 AND IANG                                                           MUL11790
C                                                                       MUL11800
      IF(NZ.EQ.1) GO TO 151                                             MUL11810
      DO 150 IZ2=2,NZ                                                   MUL11820
      ZSQ=STORE+ZVALSQ(IZ2)                                             MUL11830
      ZSQRT=SQRT(ZSQ)                                                   MUL11840
      ZMUL=ZSQRT*STORM                                                  MUL11850
      ZSQ=AVL2(EQRAD(NR1),EQRAD(NR2),ZSQRT)/ZSQ                         MUL11860
      ZMULL(IZ2)=EXP(-ZMUL)*ZSQ                                         MUL11870
  150 CONTINUE                                                          MUL11880
C                                                                       MUL11890
C FOR Z=0 WE USE ZERO FOR 1/L**2 FOR ELEMENTS WHOSE                     MUL11900
C SEPARATION IS ZERO.                                                   MUL11910
C                                                                       MUL11920
  151 IZ2=1                                                             MUL11930
      IF(STORE.LT.ALIMIT) GO TO 152                                     MUL11940
      ZSQRT=SQRT(STORE)                                                 MUL11950
      ZMUL=ZSQRT*STORM                                                  MUL11960
      ZSQ=AVL2(EQRAD(NR1),EQRAD(NR2),ZSQRT)/STORE                       MUL11970
      GO TO 153                                                         MUL11980
  152 ZSQ=ZERO(NR1)                                                     MUL11990
      ZMUL=MUS(I1)/SQRT(ZSQ)                                            MUL12000
  153 CONTINUE                                                          MUL12010
      ZMULL(IZ2)=EXP(-ZMUL)*ZSQ                                         MUL12020
      ZMULA=0.                                                          MUL12030
      ZMULD=0.                                                          MUL12040
      ZMULE=0.                                                          MUL12050
      ZMULB=0.                                                          MUL12060
      DO 155 IZ2=1,NZ                                                   MUL12070
      ZMULA=ZMULA+ZMULL(IZ2)*AFACT(IZ2)                                 MUL12080
      ZMULD=ZMULD+ZMULL(IZ2)*DFACT(IZ2)                                 MUL12090
      ZMULE=ZMULE+ZMULL(IZ2)*EFACT(IZ2)                                 MUL12100
      ZMULB=ZMULB+ZMULL(IZ2)*BFACT(IZ2)                                 MUL12110
  155 CONTINUE                                                          MUL12120
C                                                                       MUL12130
C STEP THROUGH THE SCATTERING ANGLES                                    MUL12140
C                                                                       MUL12150
      DO 200 IAZI=1,nazimul
C                                                                       MUL12170
C RESET NTEMP1, AND DEFINE CONSTANTS FOR REFERENCE TO ARRAY             MUL12180
C ALEN2                                                                 MUL12190
C                                                                       MUL12200
      NTEMP1=NOM1                                                       MUL12210
      NALEN1=(nazimul-1)*(NOM1-1)+(IAZI-1)*NOMEG1
      NALEN2=(nazimul-1)*(NOM2-1)+(IAZI-1)*NOMEG2
      do iang=1,nangmul
	   TEMA1(iang)=0.0
         TEMA2(iang)=0.0
         TEMD1(iang)=0.0
         TEMD2(iang)=0.0
         TEME1(iang)=0.0
         TEME2(iang)=0.0
         TEMB1(iang)=0.0
         TEMB2(iang)=0.0
	END DO
C                                                                       MUL12320
C STEP THROUGH NOMEG1 ELEMENTS IN FIRST RING                            MUL12330
C                                                                       MUL12340
      DO 130 K1=1,NOMEG1                                                MUL12350
C THE VALUE OF NTEMP2 MUST REMAIN BETWEEN THE VALUES OF NOM2            MUL12360
C AND NOM2+NOMEG2-1 , IN ORDER FOR IT TO ALWAYS REPRESENT ELEMENTS IN THMUL12370
C 2ND RING                                                              MUL12380
C                                                                       MUL12390
      NA=NTEMP2-NOM2                                                    MUL12400
      NA=NA/NOMEG2                                                      MUL12410
      NTEMP2=NTEMP2-NA*NOMEG2                                           MUL12420
C                                                                       MUL12430
C DEFINE THE REFERENCE  INTEGERS FOR ALEN2 FOR 1ST AND 2ND RINGS        MUL12440
C                                                                       MUL12450
      NSERC1=NALEN1+NTEMP1                                              MUL12460
      NSERC2=NALEN2+NTEMP2                                              MUL12470
C                                                                       MUL12480
C DEFINE PROBABILITIES OF ELEMENT BEING IN IN BEAM AND OF BEING SEEN    MUL12490
C BY DETECTOR                                                           MUL12500
C                                                                       MUL12510
      RT1T1=NTEST(NTEMP1)                                               MUL12520
      RT1T2=NTEST(NTEMP2)                                               MUL12530
      NT2T1=NT2(NSERC1)                                                 MUL12540
      NT2T2=NT2(NSERC2)                                                 MUL12550
C                                                                       MUL12560
C PRODUCTS WITH FIRST ELEMENT AS PRIMARY                                MUL12570
C                                                                       MUL12580
      do iang=1,nangmul
	   PROD1A=ALEN1(NTEMP1)*ALEN2(iang,NSERC2)
         PROD1D=PROD1A*NT2T2
         PROD1E=PROD1A*RT1T1
         PROD1B=PROD1A*NT2T2*RT1T1
C                                                                       MUL12630
C ADD TO TEMPORARY SUMS                                                 MUL12640
C                                                                       MUL12650
         TEMA1(iang)=TEMA1(iang)+PROD1A
         TEMD1(iang)=TEMD1(iang)+PROD1D
         TEME1(iang)=TEME1(iang)+PROD1E
         TEMB1(iang)=TEMB1(iang)+PROD1B
	end do
C                                                                       MUL12700
C ADD ARPROD TO APPROPRIATE AREA PRODUCT SUM                            MUL12710
C                                                                       MUL12720
      AREAA(IAZI,1)=AREAA(IAZI,1)+ARPROD                                MUL12730
      AREAD(IAZI,1)=AREAD(IAZI,1)+ARPROD*NT2T2                          MUL12740
      AREAE(IAZI,1)=AREAE(IAZI,1)+ARPROD*RT1T1                          MUL12750
      AREAB(IAZI,1)=AREAB(IAZI,1)+ARPROD*NT2T2*RT1T1                    MUL12760
C                                                                       MUL12770
C IF NR1.EQ.NR2 THERE IS NO NEED TO ADD TO SECOND TEMPORARY             MUL12780
C SUMS AS WELLSINCE THESE ADDITIONS WOULD BE INDISINGIUSHABLE           MUL12790
C FROM THOSE TO THE FIRST SUMS                                          MUL12800
C                                                                       MUL12810
  116 IF(NR1.EQ.NR2) GO TO 118                                          MUL12820
C                                                                       MUL12830
C PRODUCTS WITH SECOND ELEMENT AS PRIMARY                               MUL12840
C                                                                       MUL12850
      do iang=1,nangmul
	   PROD2A=ALEN1(NTEMP2)*ALEN2(iang,NSERC1)
         PROD2D=PROD2A*NT2T1
         PROD2E=PROD2A*RT1T2
         PROD2B=PROD2A*NT2T1*RT1T2
         TEMA2(iang)=TEMA2(iang)+PROD2A
         TEMD2(iang)=TEMD2(iang)+PROD2D
         TEME2(iang)=TEME2(iang)+PROD2E
         TEMB2(iang)=TEMB2(iang)+PROD2B
	end do
C                                                                       MUL12940
C ADD ARPROD TO APPROPRIATE AREA PRODUCTS                               MUL12950
C                                                                       MUL12960
      AREAA(IAZI,2)=AREAA(IAZI,2)+ARPROD                                MUL12970
      AREAD(IAZI,2)=AREAD(IAZI,2)+ARPROD*NT2T1                          MUL12980
      AREAE(IAZI,2)=AREAE(IAZI,2)+ARPROD*RT1T2                          MUL12990
      AREAB(IAZI,2)=AREAB(IAZI,2)+ARPROD*NT2T1*RT1T2                    MUL13000
  118 CONTINUE                                                          MUL13010
C                                                                       MUL13020
C INCREMENT NTEMP1 AND INCREMENT NTEMP2 SO THAT THE CORRESPONDING       MUL13030
C VALUE OF OM2 STARTS AT THE SAME VALUE AS OM1                          MUL13040
C                                                                       MUL13050
      NTEMP1=NTEMP1+1                                                   MUL13060
      NTEMP2=NTEMP2+NRAT                                                MUL13070
  130 CONTINUE                                                          MUL13080
C                                                                       MUL13090
C MULTIPLY THE TEMPORARY SUMS BY THE AREA PRODUCTS                      MUL13100
C AND SUM OVER Z VALUES                                                 MUL13110
C                                                                       MUL13120
      do iang=1,nangmul
	   TEMA1(IANG)=TEMA1(IANG)*ARPROD
         TEMA2(IANG)=TEMA2(IANG)*ARPROD
         TEMD1(IANG)=TEMD1(IANG)*ARPROD
         TEMD2(IANG)=TEMD2(IANG)*ARPROD
         TEME1(IANG)=TEME1(IANG)*ARPROD
         TEME2(IANG)=TEME2(IANG)*ARPROD
         TEMB1(IANG)=TEMB1(IANG)*ARPROD
         TEMB2(IANG)=TEMB2(IANG)*ARPROD
         SUMHA(IANG,IAZI,1)=SUMHA(IANG,IAZI,1)+ZMULA*TEMA1(IANG)
         SUMHD(IANG,IAZI,1)=SUMHD(IANG,IAZI,1)+ZMULD*TEMD1(IANG)
         SUMHE(IANG,IAZI,1)=SUMHE(IANG,IAZI,1)+ZMULE*TEME1(IANG)
         SUMHB(IANG,IAZI,1)=SUMHB(IANG,IAZI,1)+ZMULB*TEMB1(IANG)
         SUMHA(IANG,IAZI,2)=SUMHA(IANG,IAZI,2)+ZMULA*TEMA2(IANG)
         SUMHD(IANG,IAZI,2)=SUMHD(IANG,IAZI,2)+ZMULD*TEMD2(IANG)
         SUMHE(IANG,IAZI,2)=SUMHE(IANG,IAZI,2)+ZMULE*TEME2(IANG)
         SUMHB(IANG,IAZI,2)=SUMHB(IANG,IAZI,2)+ZMULB*TEMB2(IANG)
	end do
  200 CONTINUE                                                          MUL13290
C                                                                       MUL13300
C INCREMENT NTEMP2                                                      MUL13310
C                                                                       MUL13320
      NTEMP2=NTEMP2+1                                                   MUL13330
!      if(i1.eq.1.and.i2.le.2.and.k2.le.nomeg2) then
!            write(6,*) i1,i2,j1,j2,nr1,nr2,k2
!            write(6,*) storsq(k2),storsm(k2),sumr1,sumr2
!            write(6,*) zmula,zmuld,zmule,zmulb
!            write(6,*) TEMA1(1),temd1(1),teme1(1),temb1(1)
!            write(6,*) TEMA2(1),temd2(1),teme2(1),temb2(1)
!      end if
  140 CONTINUE                                                          MUL13340
C                                                                       MUL13350
C HAVING COMPLETED ALL INTEGRATIONS FOR THE CURRENT 2ND                 MUL13360
C ELEMENT, WE CAN LOOK AT NEXT 2ND AND INCREMENT NR2 AND NOM2           MUL13370
C                                                                       MUL13380
      NR2=NR2+1                                                         MUL13390
      NOM2=NOM2+NOMEG2                                                  MUL13400
   40 CONTINUE                                                          MUL13410
C                                                                       MUL13420
C HAVING NOW EXHAUSTED ALL VALUES OF R2 FOR THIS VALUE OF R1            MUL13430
C WE CAN INCREMENT NOM1 AND PROCEED TO A NEW VALUE OF R1                MUL13440
C                                                                       MUL13450
      NOM1=NOM1+NOMEG1                                                  MUL13460
   20 CONTINUE                                                          MUL13470
C ADD SUMH>S                                                            MUL13480
C                                                                       MUL13490
	do iang=1,nangmul
      DO 160 IAZI=1,nazimul
      SUMA1=SUMHA(IANG,IAZI,1)*ZSTEPS/HIGHCC
      SUMA2=SUMHA(IANG,IAZI,2)*ZSTEPS/HIGHCC
      SUMD1=SUMHD(IANG,IAZI,1)*ZSTEPS/HIGHCP
      SUMD2=SUMHD(IANG,IAZI,2)*ZSTEPS/HIGHCP
      SUME1=SUMHE(IANG,IAZI,1)*ZSTEPS/HIGHPC
      SUME2=SUMHE(IANG,IAZI,2)*ZSTEPS/HIGHPC
      SUMB1=SUMHB(IANG,IAZI,1)*ZSTEPS/HIGHPP
      SUMB2=SUMHB(IANG,IAZI,2)*ZSTEPS/HIGHPP
  500 FORMAT(1X,5(2X,E13.6))                                            MUL13590
      IF(I1.ne.I2) then
	   IF(AREAa(IAZI,1).ne.0.) then
            SUMA(IANG,IAZI,I1,I2)=SUMA(IANG,IAZI,I1,I2)
     1+SUMA1/AREAA(IAZI,1)
	   endif
	   if(areaa(IAZI,2).ne.0.) then
            SUMA(IANG,IAZI,I2,I1)=SUMA(IANG,IAZI,I2,I1)
     1+SUMA2/AREAA(IAZI,2)
         end if
	   IF(AREAD(IAZI,1).ne.0.) then
            SUMD(IANG,IAZI,I1,I2)=SUMD(IANG,IAZI,I1,I2)
     1+SUMD1/AREAD(IAZI,1)
         end if
	   IF(AREAD(IAZI,2).ne.0.) then
            SUMD(IANG,IAZI,I2,I1)=SUMD(IANG,IAZI,I2,I1)
     1+SUMD2/AREAD(IAZI,2)
         endif
  	   IF(AREAE(IAZI,1).ne.0.) then
     	      SUME(IANG,IAZI,I1,I2)=SUME(IANG,IAZI,I1,I2)
     1+SUME1/AREAE(IAZI,1)
         endif
	   if(AREAE(IAZI,2).ne.0.) then
	      SUME(IANG,IAZI,I2,I1)=SUME(IANG,IAZI,I2,I1)
     1+SUME2/AREAE(IAZI,2)
	   endif
         IF(AREAB(IAZI,1).ne.0.) then
            SUMB(IANG,IAZI,I1,I2)=SUMB(IANG,IAZI,I1,I2)
     1+SUMB1/AREAB(IAZI,1)
         endif
         IF(AREAB(IAZI,2).ne.0.) then
            SUMB(IANG,IAZI,I2,I1)=SUMB(IANG,IAZI,I2,I1)
     1+SUMB2/AREAB(IAZI,2)
         endif
C                                                                       MUL13820
C                                                                       MUL13830
C                                                                       MUL13840
	else
         SUMAR=AREAA(IAZI,1)+AREAA(IAZI,2)
         IF(SUMAR.NE.0.) then
            SUMA(IANG,IAZI,I1,I2)=(SUMA1+SUMA2)/sumar
 	   endif
         SUMAR=AREAD(IAZI,1)+AREAD(IAZI,2)
         IF(SUMAR.NE.0.) then
            SUMD(IANG,IAZI,I1,I2)=(SUMD1+SUMD2)/SUMAR
	   endif
         SUMAR=AREAE(IAZI,1)+AREAE(IAZI,2)
         IF(SUMAR.NE.0.) then
            SUME(IANG,IAZI,I1,I2)=(SUME1+SUME2)/SUMAR
	   endif
         SUMAR=AREAB(IAZI,1)+AREAB(IAZI,2)
         IF(SUMAR.NE.0.) then
            SUMB(IANG,IAZI,I1,I2)=(SUMB1+SUMB2)/SUMAR
	   endif
	endif
  160 CONTINUE                                                          MUL14000
	end do
!      write(6,*) i1,i2,nrad(i1),nrad(i2)
!      write(6,*) zmula,zmuld,zmule,zmulb
!      write(6,*) SUMA1,SUMd1,SUMe1,SUMb1
!      write(6,*) SUMA2,SUMd2,SUMe2,SUMb2
!      write(6,*) TEMA1(1),temd1(1),teme1(1),temb1(1)
!      write(6,*) TEMA2(1),temd2(1),teme2(1),temb2(1)
!      write(6,*) SUMHA(1,1,1),SUMHd(1,1,1),SUMHe(1,1,1),SUMHb(1,1,1)
!      write(6,*) SUMHA(1,1,2),SUMHd(1,1,2),SUMHe(1,1,2),SUMHb(1,1,2)
!      write(6,*) highcc,highcp,highpc,highpp
!      write(6,*) areaa(1,1),aread(1,1),areae(1,1),areab(1,1)
!      write(6,*) areaa(1,2),aread(1,2),areae(1,2),areab(1,2)
!      write(6,*) suma(1,4,i1,i2),sumd(1,4,i1,i2)
!     1,sume(1,4,i1,i2),sumb(1,4,i1,i2)
   30 CONTINUE                                                          MUL14010
   10 CONTINUE                                                          MUL14020
      RETURN                                                            MUL14030
      END                                                               MUL14040

C                                                                       MUL06340
C       ***********************                                         MUL06350
C       *                     *                                         MUL06360
C       *  SUBROUTINE DATAOP  *                                         MUL06370
C       *                     *                                         MUL06380
C       ***********************                                         MUL06390
C                                                                       MUL06400
      SUBROUTINE DATAOP(SUMC1,SUMS1,SUMT1,SUMP1,SUMA1,SUMD1,SUME1
     1,SUMB1,HIGHB,HIGHS,HIGHBS,K)
	include 'dimension.inc'
	include 'beam.inc'
	include 'mul_corr_azi.inc'
      DIMENSION SUMA(mcorrang),SUMB(mcorrang)
     *,SUMC(mcorrang),SUMP(mcorrang)
     1,DIVCC(mcorrang),DIVCP(mcorrang)
     *,SUMA1(mcorrang,mcorrang,mcont,mcont)
     2,SUMB1(mcorrang,mcorrang,mcont,mcont)
     *,SUMC1(mcorrang,mcorrang,mcont)
     *,SUMP1(mcorrang,mcorrang,mcont),BIGDEL(mcorrang,mcorrang)
     1,SMADEL(mcorrang,mcorrang)
     3,AREAB(mcont),AREAS(mcont),AREABS(mcorrang,mcont),ARBS(mcorrang)
     *,AREA(mcont)
     4,SUMD1(mcorrang,mcorrang,mcont,mcont)

     *,SUME1(mcorrang,mcorrang,mcont,mcont)
     *,SUMS1(mcorrang,mcorrang,mcont),SUMT1(mcorrang,mcorrang,mcont)
     1,SUMD(mcorrang)
     5,SUME(mcorrang),SUMS(mcorrang),SUMT(mcorrang)
      COMMON/ANNBLO/ NAN,RAD1(mcont),rad2(mcont),NRAD(mcont),RST(mcont)
     *,DEN(mcont),SIGS(mcont),MUS(mcont)
      COMMON/CHRIS/CSUMPR(mcorrwav,mcorrang,mcorrang)
     1,CSUMB(mcorrwav,mcorrang,mcorrang)
     *,CTOSCA(mcorrwav,mcorrang,mcorrang)
     2,CTOSC(mcorrwav,mcorrang,mcorrang)
      REAL MUS                                                          MUL06510
      PI=3.141592653                                                    MUL06520
C                                                                       MUL06530
C GENERATE THE EQUIVALENT SAMPLE AREAS IN THE CALCULATION - THIS        MUL06540
C MUST BE DONE NUMERICALLY SO AS TO INCORPORATE THE BEAM PROFILE
C                                                                       MUL06560
      DO 40 I=1,NAN                                                     MUL06570
      FACT=SIGS(I)*DEN(I)/(4*PI)                                        MUL06580
      IF(FACT.GT.0.) CALL EQUIV(RAD1(I),RAD2(I),100
     1,AREA(I),AREAB(I),AREAS(I),ARBS)
      AREA(I)=AREA(I)*FACT                                              MUL06620
      AREAB(I)=AREAB(I)*FACT                                            MUL06630
      AREAS(I)=AREAS(I)*FACT                                            MUL06640
      DO 45 KK=1,nazimul
      AREABS(KK,I)=ARBS(KK)*FACT
   45 CONTINUE
   40 CONTINUE                                                          MUL06670
c	write(6,*) (areaBS(kk,1),kk=1,nazimul)
c	write(6,*) (sumP1(1,kk,1),kk=1,nazimul)
c	write(6,*) (sumB1(1,kk,1,1),kk=1,nazimul)
C                                                                       MUL06680
C STEP THROUGH THE ANGLES
C
	DO IANG=1,NANGMUL
C                                                                       MUL06690
C A) PRIMARY SCATTERING.                                                MUL06700
C                                                                       MUL06710
      DO 50 IAZI=1,nazimul
      SUMC(IAZI)=0                                                      MUL06730
      SUMS(IAZI)=0.                                                     MUL06740
      SUMT(IAZI)=0.                                                     MUL06750
      SUMP(IAZI)=0                                                      MUL06760
      DO 60 IAN=1,NAN                                                   MUL06770
      SUMC(IAZI)=SUMC(IAZI)+SUMC1(IANG,IAZI,IAN)*HEIGHT*AREA(IAN)
      SUMS(IAZI)=SUMS(IAZI)+SUMS1(IANG,IAZI,IAN)*HIGHS*AREAS(IAN)
      SUMT(IAZI)=SUMT(IAZI)+SUMT1(IANG,IAZI,IAN)*HIGHB*AREAB(IAN)
      SUMP(IAZI)=SUMP(IAZI)+SUMP1(IANG,IAZI,IAN)*HIGHBS*AREABS(IAZI,IAN)
   60 CONTINUE                                                          MUL06820
   50 CONTINUE
C	write(6,*) (sumP(kk),kk=1,nazimul)
C                                                                       MUL06840
C B) SECONDARY SCATTERING.                                              MUL06850
C                                                                       MUL06860
      DO 70 IAZI=1,nazimul
      SUMA(IAZI)=0.                                                     MUL06880
      SUMD(IAZI)=0.                                                     MUL06890
      SUME(IAZI)=0.                                                     MUL06900
      SUMB(IAZI)=0.                                                     MUL06910
   70 CONTINUE                                                          MUL06920
C                                                                       MUL06930
C STEP THROUGH ANNULI AS PRIMARY SCATTER                                MUL06940
C                                                                       MUL06950
      DO 72 I1=1,NAN                                                    MUL06960
C                                                                       MUL06970
C STEP THROUGH ANNULI AS SECONDARY SCATTERER                            MUL06980
C                                                                       MUL06990
      DO 74 I2=1,NAN                                                    MUL07000
      ARPROA=HEIGHT*HEIGHT*AREA(I1)*AREA(I2)                            MUL07010
      ARPROD=HEIGHT*HIGHS*AREA(I1)*AREAS(I2)                            MUL07020
      ARPROE=HIGHB*HEIGHT*AREAB(I1)*AREA(I2)                            MUL07030
      ARPROB=HIGHB*HIGHS*AREAB(I1)*AREAS(I2)                            MUL07040
      DO 76 IAZI=1,nazimul
      SUMA(IAZI)=SUMA(IAZI)+ARPROA*SUMA1(IANG,IAZI,I1,I2)
      SUMD(IAZI)=SUMD(IAZI)+ARPROD*SUMD1(IANG,IAZI,I1,I2)
      SUME(IAZI)=SUME(IAZI)+ARPROE*SUME1(IANG,IAZI,I1,I2)
      SUMB(IAZI)=SUMB(IAZI)+ARPROB*SUMB1(IANG,IAZI,I1,I2)
   76 CONTINUE                                                          MUL07100
!                    write(6,*) i1,i2
!                    write(6,*) suma(1),sumd(1),sume(1),sumb(1)
   74 CONTINUE                                                          MUL07110
   72 CONTINUE                                                          MUL07120
C                                                                       MUL07130
C OUTPUT RESULTS AND FORM SOME RATIOS - DIVCP AND DIVCC                 MUL07140
C                                                                       MUL07150
      DO 7 I=1,nazimul
      DIVCC(I)=SUMA(I)/SUMC(I)                                          MUL07200
      DIVCP(I)=SUMD(I)/SUMS(I)                                          MUL07210
    7 CONTINUE                                                          MUL07220
      DO 11 IAZI=1,nazimul
      SUMPRI=SUMP(IAZI)                                                 MUL07250
C                                                                       MUL07260
C CALCULATE TOTAL MULTIPLE SCATTERING                                   MUL07270
C                                                                       MUL07280
      THIRSC=SUME(IAZI)*DIVCP(IAZI)/(1.-DIVCC(IAZI))                    MUL07290
      TOSCA=SUMB(IAZI)+THIRSC                                           MUL07300
      TOSCA1=TOSCA+SUMPRI                                               MUL07310
      CSUMPR(K,IANG,IAZI)=SUMPRI
      CSUMB(K,IANG,IAZI)=SUMB(IAZI)
      CTOSCA(K,IANG,IAZI)=TOSCA
      CTOSC(K,IANG,IAZI)=TOSCA1
   11 CONTINUE
	END DO
      RETURN                                                            MUL07350
      END                                                               MUL07440

	function f(amu,amu_out,t,sec1,sec2)
	s=t*(amu*sec1-amu_out*sec2)
	s1=amu_out*t*sec2
	if(abs(s).lt.0.01)then
c use a simple series expansion
            nterms=5
            f=1-s/real(nterms)
            do while (nterms.gt.2)
               nterms=nterms-1
               f=1.0-s*f/real(nterms)
            end do
	else
		f=(1.0-exp(-s))/s
	endif
	if(s1.gt.0.0) then
		f=f*exp(-s1)
	endif
	return
	end

      function binte1(amu,x0,delx)
      x1=x0-delx
      x2=x0+delx
      binte1=ainte1(amu,x1)+ainte1(amu,x2)-2.*ainte1(amu,x0)
      binte1=binte1/(delx*delx)
      return
      end

      function ainte1(amu,x0)
      real*8 expint,damux
      x=abs(x0)
      if(x.eq.0.0) go to 10
      amux=amu*x
      damux=amux
      sum=amux*amux*expint(damux)
      sum=sum+1.0+2.0*amux
      if(amux.lt.60.0) sum=sum+(1.0-amux)*exp(-amux)
      sum=sum*0.5/(amu*amu)
      ainte1=sum
      return
   10 continue
      ainte1=1./(amu*amu)
      return
      end

      function e1int1(b,al)
        real*8 expint,da1,dal,expintiexp
        dal=al
      if(b.le.1.) go to 10
      b1=b-1
      sum=0.0
      a1=b1*al
      a2=b*al
      if(a2.lt.40.0) then
         sum=log(b1)*exp(-a2)
      endif
      sum=sum-expintiexp(a1,a2)-expint(dal)
      e1int1=-sum
      return
   10 continue
      gam=0.5772156649
      sum=-gam-alog(al)
      a2=b*al
      sum=sum*exp(-a2)-expint(dal)
      e1int1=-sum
      return
   20 continue
	if(a2.lt.30.0) then
		sum=sum*exp(-a2)
	else
		sum=0.
	endif
      sum=sum-exp(a1-a2)/a1
      sum=sum-exp(-al)/al
      e1int1=-sum
      return
      end

      function e1int2(b,al)
      real*8 expint,da1,dal
      b1=b+1
      sum=alog(b1)
      a1=b1*al
      da1=a1
      dal=al


      if(a1.gt.100.0) go to 5
      sum=sum+expint(da1)
      a2=b*al
      sum=sum-exp(-a2)*expint(dal)
    5 continue
      e1int2=sum
      return
      end

      function aintex(al,sec,amu,del)
      amux1=0.5*del*amu*abs(sec)
      amux2=al*amu*abs(sec)
      sum=0.0
      if(amux2.le.40.0) then
         sum=exp(-amux2)
         if(amux1.le.40.0) then
            sum=sum*2.0*sinh(amux1)
         end if
      endif
      aintex=sum/(sec*amu*del)
cwrite(6,100) al,z1,z2,sec,aintex
  100 format(1x,6e13.6)
      return
      end

      function prime(sec0,sec1,amu,t)
      sum=1.
      amut=amu*t
      if(sec0.eq.sec1) go to 10
      s=amut*(sec0-sec1)
	if(abs(s).lt.0.01)then
c use a simple series expansion
            nterms=5
            f=1-s/real(nterms)
            do while (nterms.gt.2)
               nterms=nterms-1
               f=1.0-s*f/real(nterms)
            end do
            sum=f
	else
		sum=(1.0-exp(-s))/s
	endif
   10 continue
      if(sec1.lt.0.0) go to 20
      amuts=amut*sec1
      prime=0.
      if(amuts.lt.50.) prime=sum*exp(-amut*sec1)
      return
   20 continue
      prime=sum
      return
      end

      double precision function expint1(x)
      implicit double precision (a-h,o-z)
      dimension p1(5),q1(4),p2(7),q2(7),p3(6),q3(6),p4(8),q4(8)
      dimension a1(8),b1(7),a2(8),b2(7),a3(6),b3(5)

      data p1
     1/-4.34981 43832 95212e+2, +4.25696 82638 59170e+2,
     2 +2.92525 18866 92055e+2, +3.98941 53870 32107e+1,
     3 +4.29312 52343 20973e+0/
      data q1
     1/+7.53585 64359 84293e+2, +5.68052 52718 98696e+2,
     2 +1.50950 38744 25131e+2, +1.88992 88395 00297e+1/
      data p2
     1/+4.65627 10797 50957e-7, +9.99979 57705 15950e-1,
     2 +9.04161 55694 63287e+0, +2.43784 08879 13167e+1,
     3 +2.30192 55939 13335e+1, +6.90522 52278 44436e+0,
     4 +4.30967 83946 93888e-1/
      data q2
     1/+1.00000 00000 00000e+0, +1.00411 64382 90545e+1,
     2 +3.24264 21069 51381e+1, +4.12807 84189 14243e+1,
     3 +2.04494 78501 37942e+1, +3.31909 21359 33016e+0,
     4 +1.03400 13040 48740e-1/
      data p3
     1/-9.99999 99990 36009e-1, -1.96304 08535 93865e+1,
     2 -1.19557 61038 37175e+2, -2.54376 33976 88951e+2,
     3 -1.47982 19500 50448e+2, -2.39099 64453 13573e+0/
      data q3
     1/+1.00000 00000 00000e+0, +2.16304 08494 23774e+1,
     2 +1.56818 43364 53856e+2, +4.62230 27156 14783e+2,
     3 +5.30685 09610 81167e+2, +1.77600 70940 35063e+2/
      data p4
     1/-8.66937 33995 10696e+0, -5.49142 26552 10851e+2,
     2 -4.21001 61535 70699e+3, -2.49301 39345 86476e+5,
     3 -1.19623 66934 92469e+5, -2.21744 62775 88454e+7,
     4 +3.89280 42131 12014e+6, -3.91546 07380 90955e+8/
      data q4
     1/+3.41718 75000 00000e+1, -1.60708 92658 72209e+3,
     2 +3.57300 29805 85081e+4, -4.83547 43616 21635e+5,
     3 +4.28559 62461 17490e+6, -2.49033 37574 05403e+7,
     4 +8.91925 76757 56121e+7, -1.65254 29972 52109e+8/
      data a1
     1/+1.00443 10922 80779e+0, -4.32531 13287 81346e+1,
     2 +6.01217 99083 00805e+1, -3.31842 53199 72211e+1,
     3 +2.50762 81129 35598e+1, +9.30816 38566 21651e+0,
     4 -2.19010 23385 48809e+1, -2.18086 38152 07237e+0/
      data b1
     1/+5.27468 85196 29079e-1, +2.73624 11988 93281e+3,
     2 +1.43256 73812 19376e+1, +1.00367 43951 67258e+3,
     3 -6.25041 16167 18755e+0, +3.00892 64837 29152e+2,
     4 +3.93707 70185 27150e+0/
      data a2
     1/+9.99994 29607 47083e-1, -1.95022 32128 96598e+0,
     2 +1.75656 31546 96144e+0, +1.79601 68876 92516e+1,
     3 -3.23467 33030 54035e+1, -8.28561 99414 06413e+0,
     4 -1.86545 45488 33988e+1, -3.48334 65360 28526e+0/
      data b2
     1/+1.00083 86740 26391e+0, -3.43942 26689 98700e+0,
     2 +2.89516 72792 51351e+1, +7.60761 14800 77346e+2,
     3 +2.57776 38423 84399e+1, +5.72837 19383 73237e+1,
     4 +6.95000 65588 74340e+1/
      data a3
     1/+1.00000 00000 70443e+0, -3.00000 77799 35772e+0,
     2 -5.02233 17461 85109e+0, -9.14830 08216 73641e+0,
     3 -1.01047 90815 76032e+1, -2.77809 28934 43810e+1/
      data b3
     1/+1.99999 99428 26009e+0, -2.99901 18065 26193e+0,
     2 -7.18975 18395 04450e+0, +2.72761 00778 77917e+0,
     3 +1.22399 93926 82269e+2/
      data x0 /0.37250 74107 81367/

      if(x .gt. 4.0) go to 1
      if(x .gt. 1.0) go to 2
      if(x .gt. 0.0) go to 3
      if(dabs(x).le. 1.0e-20) go to 4
      y=-x
      if(x .gt. -6.0) go to 5
      v=dexp(y)/x
      if(x .gt. -12.0) go to 6
      if(x .gt. -24.0) go to 7

      expint1=v*(1.0+(a3(1)+b3(1)/(a3(2)+y+b3(2)/(a3(3)
     1+y+b3(3)/(a3(4)+y+
     1 b3(4)/(a3(5)+y+b3(5)/(a3(6)+y))))))/y)
      return

    7 expint1=v*(a2(1)+b2(1)/(a2(2)+y+b2(2)/(a2(3)
     1+y+b2(3)/(a2(4)+y+b2(4)
     1 /(a2(5)+y+b2(5)/(a2(6)+y+b2(6)/(a2(7)+y+b2(7)/(a2(8)+y))))))))
      return

    6 expint1=v*(a1(1)+b1(1)/(a1(2)+y+b1(2)/(a1(3)
     1+y+b1(3)/(a1(4)+y+b1(4)
     1 /(a1(5)+y+b1(5)/(a1(6)+y+b1(6)/(a1(7)+y+b1(7)/(a1(8)+y))))))))
      return

    5 v=0.666666666666667*y-2.0
      bp=0.
      bq=0.
      dp=p4(1)
      dq=q4(1)
      do 15 i = 2,8
      ap=bp
      bp=dp
      aq=bq
      bq=dq
      dp=p4(i)-ap+v*bp
      dq=q4(i)-aq+v*bq
   15 CONTINUE
      expint1=-dlog(y/x0)-(y-x0)*(dp-ap)/(dq-aq)
      return

    4 expint1=0.
      return

    3 expint1=-dlog(x)+
     1       (p1(1)+x*(p1(2)+x*(p1(3)+x*(p1(4)+x*p1(5)))))/
     2       (q1(1)+x*(q1(2)+x*(q1(3)+x*(q1(4)+x))))
      return

    2 y=1.0/x
      expint1=dexp(-x)*
     1(p2(1)+y*(p2(2)+y*(p2(3)+y*(p2(4)+y*(p2(5)
     1+y*(p2(6)+y*p2(7)))))))/
     2(q2(1)+y*(q2(2)+y*(q2(3)+y*(q2(4)+y*(q2(5)+y*(q2(6)+y*q2(7)))))))
      return

    1 y=1.0/x
      expint1=dexp(-x)*y*(1.0+y*
     1       (p3(1)+y*(p3(2)+y*(p3(3)+y*(p3(4)+y*(p3(5)+y*p3(6))))))/
     2       (q3(1)+y*(q3(2)+y*(q3(3)+y*(q3(4)+y*(q3(5)+y*q3(6)))))))
      return
 	end

      function expint(x)
c Requires x > 0
      real*8 x,rx,a(4),b(4),c(5),expint,sum1,sum2,sum,ratio,gam,an
     1,prod
      integer n
      data a/8.5733287401,18.0590169730,8.6347608935,0.2677737343/
      data b/9.5733223454,25.6329561486,21.0996530827,3.9584969228/
      data c/0.99999193,-0.24991055,0.05519968,-0.00976004,0.00107857/
      data gam/0.5772156649/

      x=abs(x)
      expint=0.0
      if(x.le.0.0) return

      if (x.ge.1.0) then
c Expansion from A + S p 231, 5.1.56
         rx=1.0/x
         sum1=a(4)
         sum2=b(4)
         n=4
         do while (n.gt.1)
            n=n-1
            sum1=rx*sum1+a(n)
            sum2=rx*sum2+b(n)
         end do
         sum1=1.0+sum1
         sum2=1.0+sum2
         ratio=sum1/sum2
c Form the log of arguments to prevent overflows
         sum=x+log(x)-log(ratio)
c      write(6,*) sum1,sum2,ratio,sum,exp(-sum)
         if(abs(sum).lt.50.0) then
            expint=exp(-sum)
         endif

      else

c A + S p231 5.1.53

         n=5
         sum=x*c(n)
         do while (n.gt.1)
            n=n-1
            sum=x*(c(n)+sum)
         end do
         expint=sum-gam-log(x)

      endif

      return
      end

      function expintiexp(a1,a2)
      real a1,a2
      real*8 expintiexp,da1,da2,gam,sum,t,t2,prod,an
      integer n
c calculates the value of E_i(a1)*exp(-a2) using the series expansion. Requires a1,a2 > 0.0
      da1=a1
      da2=a2
      gam=0.5772156649
      n=int(sqrt(a1))+1
      an=real(n)
      t=exp(-da2/an)
      prod=da1*t
      sum=prod/(an*an)
c      write(6,*) n,prod,t,an,sum
      t2=1
      do while (n.gt.1)
         n=n-1
         an=real(n)
         t2=t2*t
         sum=prod*(t2/an+sum)/an
c         write(6,*) n,prod,t,an,sum
      end do
      if(da2.lt.40.) then
         sum=sum+(gam+log(da1))*exp(-da2)
      endif
      expintiexp=sum
      return
      end

	function cosd(deg)
	parameter (degradconv=3.141592653/180.0)
	cosd=cos(deg*degradconv)
	return
	end
	function sind(deg)
	parameter (degradconv=3.141592653/180.0)
	sind=sin(deg*degradconv)
	return
	end

	function acosd(cosphi)
	parameter (degradconv=3.141592653/180.0)
	acosd=acos(cosphi)
	acosd=acosd/degradconv
	return
	end
	function asind(sinphi)
	parameter (degradconv=3.141592653/180.0)
	asind=asin(sinphi)
	asind=asind/degradconv
	return
	end
	function atan2d(sinphi,cosphi)
	parameter (degradconv=3.141592653/180.0)
	atan2d=atan2(sinphi,cosphi)
	atan2d=atan2d/degradconv
	return
	end

