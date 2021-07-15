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
   5  ARBS(I)=0.                                                        MUL15570
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

