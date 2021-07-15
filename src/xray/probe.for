C                                                                       MUL06110
C      *********************                                            MUL06120
C      *                   *                                            MUL06130
C      *  FUNCTION PROBE   *                                            MUL06140
C      *                   *                                            MUL06150
C      *********************                                            MUL06160
C                                                                       MUL06170
      FUNCTION PROBE(X,PROFIL,NPROF,PRSTEP,A1,B1)                       MUL06180
      DIMENSION PROFIL(NPROF)                                           MUL06190
      IF(X.GT.A1.OR.X.LT.B1) GO TO 10                                   MUL06200
      DIFF=A1-X                                                         MUL06210
      APOS=DIFF/PRSTEP                                                  MUL06220
      NPOS=INT(APOS)                                                    MUL06230
      DIFF=APOS-FLOAT(NPOS)                                             MUL06240
      NPOS1=NPOS+1                                                      MUL06250
      NPOS2=NPOS+2                                                      MUL06260
      IF(DIFF.EQ.0.) GO TO 20                                           MUL06270
      PROBE=PROFIL(NPOS1)+(PROFIL(NPOS2)-PROFIL(NPOS1))*DIFF            MUL06280
      RETURN                                                            MUL06290
   20 PROBE=PROFIL(NPOS1)                                               MUL06300
	RETURN
   10 PROBE=0.                                                          MUL06310
      RETURN                                                            MUL06320
      END                                                               MUL06330
