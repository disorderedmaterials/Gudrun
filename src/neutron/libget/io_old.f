C---------------------------------------------------------------------------
C *** CCLRC ISIS Facility GET Routines ***
C *** Original routines by Kevin Knowles, modifications by Freddie Akeroyd 
C
C $Id: io.f,v 1.36 1999/02/04 13:38:31 faa Exp $
C---------------------------------------------------------------------------
        SUBROUTINE CLOSE_DATA_FILE()
C
C  Forces the currently opened RAW file to be closed
C
C
        CHARACTER*256 FILNAM
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,
     +  NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
        LOGICAL OPENFILE
        COMMON/crpt_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,OPENFILE,NPER,
     +          DATA_HEADER
        SAVE /crpt_SPECIALS/
C
        FILNAM = ' ' 
        RETURN
        END
C-------------------------------------------------------------------------
        SUBROUTINE GETDAT(RUNID,IFSN,NOS,IDATA,LENGTH,ERRCODE)
        IMPLICIT NONE
C
C *** ERROR CODES ***
C
C    0 = all OK
C    1 = file RUNID not currently open
C    2 = file number is '00000' and so no access available
C    3 = asked for non-existent parameter
C    4 = TOO MANY SPECTRA ASKED FOR
C    5 = error in byte unpacking
C    6 = cannot understand data section
C
C
C  Gets data from a RAW file one or more spectra at a time
C
        CHARACTER*(*) RUNID
	LOGICAL CRPT_ACCESS, DAE_ACCESS
        INTEGER*4 ERRCODE,IFSN,NOS,LENGTH,IDATA(LENGTH)
C
        INTEGER*4 IWORK(256),IERR
        REAL*4 WORK(256)
        COMMON/crpt_WORK/IWORK
        SAVE /crpt_WORK/
        EQUIVALENCE (IWORK,WORK)
C temporary parameters
        INTEGER*4 IBASE,ILONG
        INTEGER I,J,ISTART,STATUS,ICOMPRESS,IBUFFER(33000)
C
C
        CHARACTER*256 FILNAM
	  CHARACTER*256 MESS
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,
     +  NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
        LOGICAL OPENFILE
        COMMON/crpt_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,OPENFILE,NPER,
     +          DATA_HEADER
        SAVE /crpt_SPECIALS/
	EXTERNAL CRPT_ACCESS, DAE_ACCESS
C

        ERRCODE=0
C
C  decide whether it's CRPT or just a file
        IF (RUNID(4:8).EQ.'00000') THEN
C               WRITE(6,6000)
C 6000          FORMAT(' Accessing of current run not available')
C               GOTO 999
           ERRCODE=2
        ENDIF
C
C  check name is valid & open file according to RUNID
        IF ((FILNAM.NE.RUNID) .OR. (FILNAM.EQ.' ')) THEN
	MESS = 'Access requested to '//RUNID//' when '
     +//FILNAM//' open'
	   CALL FERROR_ADD('GETDAT',
     1     MESS,
     2     ' ')
           ERRCODE=1
           RETURN
        ENDIF
C
C
C  read data into IDAT ....remembering there are NTC1+1 time channels and
C -----------------------------------------------------------------------
C  NSP1+1 spectra and now also NPER periods
C
        IF (IFSN.LT.0.OR.IFSN.GT.((NSP1+1)*NPER)-1) THEN
C               WRITE (6,6010) NSP1,IFSN
C 6010  FORMAT(' Sorry, only ',I6,' spectra. You asked for',I6)
                ERRCODE=4
                RETURN
        ENDIF
        IF (IFSN+NOS-1.GT.((NSP1+1)*NPER)-1) THEN
C               WRITE (6,6020) NSP1,IFSN+NOS-1
C 6020  FORMAT(' Sorry, only ',I6,' spectra. You asked for up to',I6)
                ERRCODE=4
                RETURN
        ENDIF
        IF (VER1.EQ.1) THEN
                IBASE=IFORMAT(6)
        ELSE
                IBASE=IFORMAT(7)
        ENDIF
        IF ( (IVER(7).LE.1) .OR. DAE_ACCESS(RUNID) ) THEN
C Original version of data section
                IBASE = IBASE + 1 + IFSN*(NTC1+1)
                ILONG = NOS*(NTC1+1)
C       TYPE *,' BASE and LONG',IBASE,ILONG
                CALL GETSECT(IBASE,ILONG,IDATA,49,IERR)
		CALL VAX_TO_LOCAL_INTS(IDATA, ILONG)
        ELSE
C New version of data section (may be compressed and not necessarily consecutive
C First pick up data section header, in particular compression type.
CCCCCC          CALL GETSECT(IBASE+1,32,DATA_header)
                ICOMPRESS=DATA_HEADER(1)
C The CRPT has the "compress" flag set, even though it isn't really compressed
                IF ( (ICOMPRESS .EQ. 0) .OR. CRPT_ACCESS(RUNID) ) THEN 
C uncompressed
                    IBASE = IBASE + DATA_HEADER(3) + IFSN*(NTC1+1)
                    ILONG=NOS*(NTC1+1)
                    CALL GETSECT(IBASE,ILONG,IDATA,49,IERR)
		    CALL VAX_TO_LOCAL_INTS(IDATA, ILONG)
                ELSEIF (ICOMPRESS.EQ.1) THEN 
C byte relative compression
                    J=1
                    DO I=IFSN,IFSN+NOS-1
                        ISTART=IBASE + DATA_HEADER(3) + 2*I
                        ILONG=2
                        CALL GETSECT(ISTART,ILONG,IBUFFER,49,IERR)
		    	CALL VAX_TO_LOCAL_INTS(IBUFFER, ILONG)
                        ISTART=IBASE + IBUFFER(2)
                        ILONG=IBUFFER(1)
                        CALL GETSECT(ISTART,ILONG,IBUFFER,49,IERR)
C *** We do not need to CALL VAX_TO_LOCAL_INTS() as handled in byte_rel_expn()
                        call byte_rel_expn(
     +                      IBUFFER,ILONG*4,1,IDATA(J),NTC1+1,STATUS)
C *** odd number error codes are OK ***
                        IF (MOD(STATUS,2) .EQ. 0) THEN
                          ERRCODE=5
C                         TYPE *,' ERROR - decompressing spectrum',I
                          goto 999
                        ENDIF   
                        J = J + NTC1+1
                    ENDDO
                ELSE
                    ERRCODE=6
C                    TYPE *,' ERROR - cannot understand data section'
                ENDIF
        ENDIF
C
 999    RETURN
        END
C------------------------------------------------------------------
        SUBROUTINE GETPARC(RUNID,NAME,CVALUE,LENGTH_IN,LENGTH_OUT,
     1					 ERRCODE)
        IMPLICIT NONE
C
C *** ERROR CODES ***
C
C    0 = all OK
C    1 = file RUNID not currently open
C    2 = file number is '00000' and so no access available
C    3 = asked for non-existent parameter
C    4 = other error
C
C
C  Gets character parameter(s) from RAW data file
C
        INTEGER*4 LENGTH_IN,LENGTH_OUT,ERRCODE,ILINES,I,ITEMP(33),IERR,
     1			NOTESECT
        CHARACTER*(*) RUNID
        CHARACTER*(*) NAME
        CHARACTER*(*) CVALUE(LENGTH_IN)
	INTEGER VAX_TO_LOCAL_INT
	EXTERNAL VAX_TO_LOCAL_INT
C
C
        CHARACTER*132 CTEMP
C	BYTE B(4)
        EQUIVALENCE (ITEMP,CTEMP)
C
C
        CHARACTER*256 FILNAM
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,
     +  NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32),K
	INTEGER NLINES, VER9, LLEN, OFFSET, ILLEN
        LOGICAL OPENFILE
	CHARACTER*256 MESS
        COMMON/crpt_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,OPENFILE,NPER,
     +          DATA_HEADER
        SAVE /crpt_SPECIALS/
C
      ERRCODE=0
C Unless we set it otherwise, assume returning a simple string
      LENGTH_OUT = 1
C
C  decide whether it's CRPT or just a file
        IF (RUNID(4:8).EQ.'00000') THEN
C               WRITE(6,6000)
C 6000          FORMAT(' Accessing of current run not available')
C               GOTO 999
           ERRCODE=2
        ENDIF
C
C  check name is valid & open file according to RUNID
        IF ((FILNAM .NE. RUNID) .OR. (FILNAM.EQ.' ')) THEN
           MESS = 'Access requested to '//RUNID//
     1		  ' when '//FILNAM//' open'
	   CALL FERROR_ADD('GETPARC',
     1     MESS,
     2     ' ')
           ERRCODE=1
           RETURN
        ENDIF
C
C  read variables into CVALUE
C ----------------------------
C
C  Format section
        IF (NAME.EQ.'HDR') THEN
                CALL GETSECT(1,20,ITEMP,49,IERR)
                CVALUE(1)=CTEMP
C  Run section
        ELSEIF (NAME.EQ.'TITL') THEN
                CALL GETSECT(IFORMAT(1)+2,20,ITEMP,49,IERR)
                CVALUE(1)=CTEMP
        ELSEIF (NAME.EQ.'USER') THEN
                CALL GETSECT(IFORMAT(1)+22,20,ITEMP,49,IERR)
                CVALUE(1)=CTEMP(1:20)
                CVALUE(2)=CTEMP(21:40)
                CVALUE(3)=CTEMP(41:60)
                CVALUE(4)=CTEMP(61:80)
                CALL GETSECT(IFORMAT(1)+42,20,ITEMP,49,IERR)
                CVALUE(5)=CTEMP(1:20)
                CVALUE(6)=CTEMP(21:40)
                CVALUE(7)=CTEMP(41:60)
                CVALUE(8)=CTEMP(61:80)
C  Instrument section
        ELSEIF (NAME.EQ.'NAME') THEN
                CALL GETSECT(IFORMAT(2)+1,2,ITEMP,49,IERR)
                CVALUE(1)=CTEMP(1:8)
C LOG / Notes section
        ELSEIF (NAME.EQ.'NOTE') THEN
		IF (VER1 .EQ. 1) THEN
		    NOTESECT = 7
		ELSE
		    NOTESECT = 8
		ENDIF
		CALL GETSECT(IFORMAT(NOTESECT), 1, VER9, 49, IERR)
		VER9 = VAX_TO_LOCAL_INT(VER9)
C		WRITE(6,*) 'VER9 = ', ver9
C		write(6,*) 'Dumping section'
C		do i = 1, 100
C		    call getsect(iformat(notesect)+i, 1, b, 49, ierr)
C		    write(6,*) (b(j), char(b(j)), j=1,4)
C		enddo
		IF (VER9 .EQ. 0) THEN
                    ILINES=(IFORMAT(NOTESECT+1)-IFORMAT(NOTESECT))/20
		    NLINES=ILINES
		    OFFSET=IFORMAT(NOTESECT) + 2
		    ILLEN=20	! 20*4 characters
		    LLEN=80
C		    WRITE(6,*) 'Reading ',nlines, 'lines'
		    CTEMP = ' '
		    LENGTH_OUT=MIN(NLINES,LENGTH_IN)
                    DO I=1,LENGTH_OUT
		       K=OFFSET+(I-1)*ILLEN
                       CALL GETSECT(K,ILLEN,ITEMP,49,IERR)
                       CVALUE(I)=CTEMP(1:LLEN)
		    ENDDO
		ELSEIF (VER9 .EQ. 2) THEN
		    CALL GETSECT(IFORMAT(NOTESECT)+1, 1, NLINES, 49, IERR)
		    CALL VAX_TO_LOCAL_INTS(NLINES, 1)
C		    WRITE(6,*) 'Reading ',nlines, 'lines'
		    OFFSET = IFORMAT(NOTESECT) + 2
C Each line stored as a line length + data
		    LENGTH_OUT=MIN(NLINES,LENGTH_IN)
                    DO I=1,LENGTH_OUT
			CALL GETSECT(OFFSET, 1, LLEN, 49, IERR)
		        CALL VAX_TO_LOCAL_INTS(LLEN, 1)
C		        WRITE(6,*) 'Reading line length', llen
		        ILLEN = (LLEN-1)/4 + 1 
			CTEMP = ' '
                        CALL GETSECT(OFFSET+1,ILLEN,ITEMP,49,IERR)
                        CVALUE(I)=CTEMP(1:LLEN)
C			write(6,*) ctemp
			OFFSET = OFFSET + ILLEN + 1
		    ENDDO
 		ELSE
		    CALL GETSECT(IFORMAT(NOTESECT)+1, 1, NLINES, 49, IERR)
		    CALL VAX_TO_LOCAL_INTS(NLINES, 1)
		    ILLEN=20  ! 20*4 characters per line
		    LLEN=80
C		    WRITE(6,*) 'Reading ',nlines, 'lines'
		    OFFSET=IFORMAT(NOTESECT)+2
		    CTEMP = ' '
                    LENGTH_OUT=MIN(NLINES,LENGTH_IN)
                    DO I=1,LENGTH_OUT
			    K=OFFSET+(I-1)*ILLEN
                            CALL GETSECT(K,ILLEN,ITEMP,49,IERR)
                            CVALUE(I)=CTEMP(1:LLEN)
		    ENDDO
		ENDIF
                IF (NLINES.LE.0) THEN
                    CVALUE(1)=' No notes were made'
                ENDIF
		IF (LENGTH_OUT .LT. NLINES) THEN
           	    MESS = 'Not enough space to return all of NOTES section'
	            CALL FERROR_ADD('GETPARC', MESS, ' ')
		ENDIF
C  non existant requests
        ELSE
                ERRCODE=3
     	    	MESS = 'No such CHARACTER parameter as '//NAME
	        CALL FERROR_ADD('GETPARC', 
     1      		    MESS, 
     2			    ' ')
C               TYPE *,' No such CHARACTER parameter as',NAME
        ENDIF
C
C
C
 999    RETURN
        END
C
CGET routines modified to cope with periods Feb89 !!
C-------------------------------------------------------------------------
        SUBROUTINE GETPARI(RUNID,NAME,IVALUE,LENGTH_IN,
     1			   LENGTH_OUT,ERRCODE)
        IMPLICIT NONE
C
C *** ERROR CODES ***
C
C    0 = all OK
C    1 = file RUNID not currently open
C    2 = file number is '00000' and so no access available
C    3 = asked for non-existent parameter
C    4 = other error
C
C  Gets named integer paramter(s) from a RAW data file
C  Whole sections may also be requested
C
        CHARACTER*(*) RUNID
        CHARACTER*(*) NAME
        INTEGER*4 ERRCODE,LENGTH_IN,LENGTH_OUT,LENGTH,
     1            IVALUE(LENGTH_IN),IFROM,IERR
C
C
        INTEGER*4 I, J, IWORK(256), SENUM, NOTESECT
        REAL*4 WORK(256)
        COMMON/crpt_WORK/IWORK
        SAVE /crpt_WORK/
        EQUIVALENCE (IWORK,WORK)
C
        CHARACTER*256 FILNAM
	CHARACTER*256 MESS
	INTEGER NTC,NDETY,NUSE
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,K,
     +            NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
	INTEGER*4 IJUNK, OFFSET, LLEN, ILLEN
	BYTE BJUNK(4)
        LOGICAL OPENFILE
C
        COMMON/crpt_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,OPENFILE,NPER,
     +          DATA_HEADER
        SAVE /crpt_SPECIALS/
	EQUIVALENCE (IJUNK,BJUNK)
	INTRINSIC CHAR
	LENGTH=LENGTH_IN
	LENGTH_OUT=0
        ERRCODE=0
	IERR = 0
C
C  decide whether it's CRPT or just a file
        IF (RUNID(4:8) .EQ. '00000') THEN
	   	CALL FERROR_ADD('GETPARI', 'Runid is 00000', ' ')
                ERRCODE=2
                RETURN
        ENDIF
C
C  check name is valid & open file according to RUNID
C
      	CALL OPEN_DATA_FILE(RUNID,NTC,NDETY,NUSE,ERRCODE)
	IF (ERRCODE .EQ. 1) RETURN
C
C  read variables into IVALUE
C ----------------------------
C
C  From now on just decide what has been requested and return it
        IF     (NAME.EQ.'VER1') THEN
                IVALUE(1)=VER1
		LENGTH = 1
        ELSEIF (NAME.EQ.'SFMT') THEN
                CALL GETSECT(1,31,IVALUE,49,IERR)
		LENGTH = 31
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
C
C  run section
        ELSEIF (NAME.EQ.'SRUN') THEN
                CALL GETSECT(IFORMAT(1),94,IVALUE,49,IERR)
		LENGTH = 94
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'VER2') THEN
                CALL GETSECT(IFORMAT(1),1,IVALUE,49,IERR)
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'RUN') THEN
                CALL GETSECT(IFORMAT(1)+1,1,IVALUE,49,IERR)
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF ((NAME.EQ.'RPB') .OR. (NAME.EQ.'IRPB')) THEN
                CALL GETSECT(IFORMAT(1)+62,32,IVALUE,49,IERR)
		LENGTH = 32
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
C
C  instrument section
        ELSEIF (NAME.EQ.'SINS') THEN
          IF (IVER(1).EQ.1) THEN
            CALL GETSECT(IFORMAT(2),70+(NMON*2)+(6+NEFF)*NDET,
     1                   IVALUE,49,IERR)
		LENGTH = 70+(NMON*2)+(6+NEFF)*NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
          ELSE
            CALL GETSECT(IFORMAT(2),70+(NMON*2)+(5+NEFF)*NDET,
     1                   IVALUE,49,IERR)
		LENGTH = 70+(NMON*2)+(5+NEFF)*NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
          ENDIF
        ELSEIF (NAME.EQ.'VER3') THEN
                IVALUE(1)=IVER(2)
		LENGTH = 1
        ELSEIF (NAME.EQ.'IVPB') THEN
                CALL GETSECT(IFORMAT(2)+3,64,IVALUE,49,IERR)
		LENGTH = 64
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'NDET') THEN
                IVALUE(1)=NDET
		LENGTH = 1
        ELSEIF (NAME.EQ.'NMON') THEN
                IVALUE(1)=NMON
		LENGTH = 1
        ELSEIF (NAME.EQ.'NEFF') THEN
                IVALUE(1)=NEFF
		LENGTH = 1
        ELSEIF (NAME.EQ.'NUSE') THEN
                IVALUE(1)=NEFF
		LENGTH = 1
        ELSEIF (NAME.EQ.'MDET') THEN
                CALL GETSECT(IFORMAT(2)+70,NMON,IVALUE,49,IERR)
		LENGTH = NMON
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'MONP') THEN
                CALL GETSECT(IFORMAT(2)+70+NMON,NMON,IVALUE,49,IERR)
		LENGTH = NMON
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'SPEC') THEN
                IFROM=IFORMAT(2)+70+2*NMON
                CALL GETSECT(IFROM,NDET,IVALUE,49,IERR)
		LENGTH = NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'CODE') THEN
                IF (IVER(2).NE.1) THEN
                IFROM=IFORMAT(2)+70+2*NMON+3*NDET
                CALL GETSECT(IFROM,NDET,IVALUE,49,IERR)
		LENGTH = NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
		ELSE
		IERR=1
                ENDIF
        ELSEIF (NAME.EQ.'TIMR') THEN
          IF (VER1.EQ.1) THEN
                IFROM=IFORMAT(2)+70+2*NMON+NDET
          ELSE
                IFROM=IFORMAT(4)+65+3*NDET
          ENDIF
          CALL GETSECT(IFROM,NDET,IVALUE,49,IERR)
		LENGTH = NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
C
C  sample environment section
        ELSEIF (NAME.EQ.'SSEN') THEN
          IF (IVER(3).EQ.1) THEN
                CALL GETSECT(IFORMAT(3),34+NSEP*24,IVALUE,49,IERR)
		LENGTH = 34+NSEP*24
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
          ELSE
                CALL GETSECT(IFORMAT(3),66+NSEP*32,IVALUE,49,IERR)
		LENGTH = 66+NSEP*24
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
          ENDIF
          IF (NSEP.NE.0) THEN
	   	CALL FERROR_ADD('GETPARI',
     1		'GETPAR needs adjusting to take account of SE',
     2     	' ')
             ERRCODE=4
          ENDIF
        ELSEIF (NAME.EQ.'VER4') THEN
                CALL GETSECT(IFORMAT(3),1,IVALUE,49,IERR)
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'SPB ') THEN
          IF (IVER(3).EQ.1) THEN
                CALL GETSECT(IFORMAT(3)+1,32,IVALUE,49,IERR)
		LENGTH = 32
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
          ELSE
                CALL GETSECT(IFORMAT(3)+1,64,IVALUE,49,IERR)
		LENGTH = 64
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
          ENDIF
        ELSEIF (NAME.EQ.'NSEP') THEN
                IVALUE(1)=NSEP
		LENGTH = 1
        ELSEIF (NAME(1:2) .EQ.'SE') THEN
		READ(NAME(3:4), '(I2.2)') SENUM
		IF ((IVER(3) .EQ. 1) .OR. (NSEP .LT. SENUM)) THEN
     		    MESS = 'Invalid SE block '//NAME
	   	    CALL FERROR_ADD('GETPARI',
     1		    MESS,
     2     	    ' ')
                    ERRCODE=4
		ELSE
C		    WRITE(6,'(1X,A,I2.2)') 'Reading SE',SENUM
		    K = IFORMAT(3)+34+32*SENUM
                    CALL GETSECT(K,32,IVALUE,49,IERR)
		    LENGTH = 32
		    CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
C		    DO I=1,32
C			IJUNK=IVALUE(I)
C		        WRITE(6,'(4(1X,Z2,A,A1))') (BJUNK(J), '=', CHAR(BJUNK(J)), J=1,4)
C		    ENDDO
		ENDIF
C  DAE section
        ELSEIF (NAME.EQ.'SDAE') THEN
          IF (IVER(4).EQ.1) THEN
                CALL GETSECT(IFORMAT(4),65+3*NDET,IVALUE,49,IERR)
		LENGTH = 65+3*NDET
          ELSE
                CALL GETSECT(IFORMAT(4),65+5*NDET,IVALUE,49,IERR)
		LENGTH = 65+5*NDET
          ENDIF
	  CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'VER5') THEN
                CALL GETSECT(IFORMAT(4),1,IVALUE,49,IERR)
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'DAEP') THEN
                CALL GETSECT(IFORMAT(4)+1,64,IVALUE,49,IERR)
		LENGTH = 64
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'CRAT') THEN
                CALL GETSECT(IFORMAT(4)+65,NDET,IVALUE,49,IERR)
		LENGTH = NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'MODN') THEN
                CALL GETSECT(IFORMAT(4)+65+NDET,NDET,IVALUE,49,IERR)
		LENGTH = NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'MPOS') THEN
                CALL GETSECT(IFORMAT(4)+65+2*NDET,NDET,IVALUE,49,IERR)
		LENGTH = NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'UDET') THEN
                CALL GETSECT(IFORMAT(4)+65+4*NDET,NDET,IVALUE,49,IERR)
		LENGTH = NDET
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
C
C  TCB section
        ELSEIF (NAME.EQ.'STCB') THEN
                IF (NTRG.NE.1) THEN
	   	   CALL FERROR_ADD('GETPARI',
     1		   'GETPAR needs adjusting to take account of SE',
     2     	   ' ')
                   ERRCODE=4
 		ELSE
                   CALL GETSECT(IFORMAT(5),288+NTC1+1,IVALUE,49,IERR)
		   LENGTH = 288+NTC1+1
		   CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
                ENDIF
        ELSEIF (NAME.EQ.'VER6') THEN
                CALL GETSECT(IFORMAT(5),1,IVALUE,49,IERR)
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'NTRG') THEN
                IF (NTRG.NE.1) THEN
	   	   CALL FERROR_ADD('GETPARI',
     1             'Multiple time regimes....GETPAR needs changing',
     2     	   ' ')
                   ERRCODE=4
		ELSE
                   IVALUE(1)=NTRG
		   LENGTH = 1
                ENDIF
        ELSEIF (NAME.EQ.'NFPP') THEN
                CALL GETSECT(IFORMAT(5)+2,1,IVALUE,49,IERR)
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'NPER') THEN
                IVALUE(1)=NPER
		LENGTH = 1
        ELSEIF (NAME.EQ.'PMAP') THEN
                CALL GETSECT(IFORMAT(5)+4,256,IVALUE,49,IERR)
		LENGTH = 256
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'NSP1') THEN
                IVALUE(1)=NSP1
		LENGTH = 1
        ELSEIF (NAME.EQ.'NTC1') THEN
                IVALUE(1)=NTC1
		LENGTH = 1
        ELSEIF (NAME.EQ.'TCM1') THEN
                CALL GETSECT(IFORMAT(5)+262,5,IVALUE,49,IERR)
		LENGTH = 5
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'PRE1') THEN
                CALL GETSECT(IFORMAT(5)+287,1,IVALUE,49,IERR)
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'TCB1') THEN
                CALL GETSECT(IFORMAT(5)+288,NTC1+1,IVALUE,49,IERR)
		LENGTH = NTC1+1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
        ELSEIF (NAME.EQ.'DHDR') THEN
		DO I=1, 32
		    IVALUE(I) = DATA_HEADER(I)
		ENDDO
		LENGTH = 32
	ELSEIF (NAME .EQ. 'ULEN') THEN
		IVALUE(1) = ULEN
		LENGTH = 1
C User section
        ELSEIF (NAME.EQ.'VER7') THEN
		IF (VER1 .EQ. 1) THEN
	            CALL FERROR_ADD('GETPARI', 
     1      		    'No USER section, so no VER7 parameter', 
     2			    ' ')
		ELSE
                    CALL GETSECT(IFORMAT(6),1,IVALUE,49,IERR)
		    LENGTH = 1
		    CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
		ENDIF
C data section section
        ELSEIF (NAME.EQ.'VER8') THEN
		IF (VER1 .EQ. 1) THEN
                    CALL GETSECT(IFORMAT(6),1,IVALUE,49,IERR)
		ELSE
                    CALL GETSECT(IFORMAT(7),1,IVALUE,49,IERR)
		ENDIF
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
C NOTES section
        ELSEIF (NAME.EQ.'VER9') THEN
		IF (VER1 .EQ. 1) THEN
C Don't think this version number exists
C		    IVALUE(1) = 1
                    CALL GETSECT(IFORMAT(7),1,IVALUE,49,IERR)
		ELSE
                    CALL GETSECT(IFORMAT(8),1,IVALUE,49,IERR)
		ENDIF
		LENGTH = 1
		CALL VAX_TO_LOCAL_INTS(IVALUE, LENGTH)
C Max Line length in notes section (bytes)
        ELSEIF (NAME.EQ.'NTLL') THEN
		IF (VER1 .EQ. 1) THEN
		    NOTESECT = 7
		ELSE
		    NOTESECT = 8
		ENDIF
		IF (IVER(8) .LT. 2) THEN
		    IVALUE(1) = 80
		ELSE
C Get number of lines into J
  		    CALL GETSECT(IFORMAT(NOTESECT)+1, 1, J, 49, IERR)
		    OFFSET = IFORMAT(NOTESECT) + 2
C Each line stored as a line length + data
		    IF (J .LT. 1) THEN
			IVALUE(1) = 80
		    ELSE
		        IVALUE(1) = 0
		    ENDIF
                    DO I=1,J
			CALL GETSECT(OFFSET, 1, LLEN, 49, IERR)
			ILLEN = (LLEN - 1)/4 + 1
			IVALUE(1) = MAX(IVALUE(1), ILLEN*4)
			OFFSET = OFFSET + ILLEN + 1
		    ENDDO
                ENDIF
		LENGTH = 1
        ELSEIF (NAME.EQ.'FORM') THEN
		IVALUE(1) = IFORMAT(10)
		LENGTH = 1
C Number of lines in notes section
        ELSEIF (NAME.EQ.'NTNL') THEN
		IF (VER1 .EQ. 1) THEN
		    NOTESECT = 7
		ELSE
		    NOTESECT = 8
		ENDIF
		IF (IVER(8) .EQ. 0) THEN
                    IVALUE(1)=(IFORMAT(NOTESECT+1)-IFORMAT(NOTESECT))/20
                ELSE
		    CALL GETSECT(IFORMAT(NOTESECT)+1, 1, IVALUE, 49, IERR)
		    CALL VAX_TO_LOCAL_INTS(IVALUE, 1)
		ENDIF
		IF (IVALUE(1) .LT. 1) IVALUE(1) = 1
		LENGTH = 1

C  non existant requests

        ELSE
                ERRCODE=3
		LENGTH = 0
C     	    	MESS = 'No such INTEGER parameter as '//NAME
C	        CALL FERROR_ADD('GETPARI', 
C     1      		    MESS,
C     2			    ' ')
C               TYPE *,' No such INTEGER parameter as',NAME
        ENDIF
	IF ((IERR .NE. 0) .AND. (ERRCODE .EQ. 0)) THEN
     	    MESS = 'Error in reading data from file '//RUNID
	    CALL FERROR_ADD('GETPARI', 
     1      		    MESS,
     2			    ' ')
	    ERRCODE = 4
	ENDIF
C
	LENGTH_OUT=LENGTH
        RETURN
        END
C----------------------------------------------------------------------
        SUBROUTINE GETRUN(RUNID,IARRAY,LENGTH,IERROR)
        IMPLICIT NONE
C
C IERROR = 0 : ALL OK
C IERROR = 1 : FILE OPEN ERROR
C if error in GETDAT, returns 10+(error code in getdat)
C
C  This routine returns the whole of a run file into the given array
C
        CHARACTER*(*) RUNID
        INTEGER*4 LENGTH,IERROR,ERRCODE,IERR
        INTEGER*4 IARRAY(LENGTH)
C
        LOGICAL FOUND
        INTEGER start_of_data
C
C
        CHARACTER*256 FILNAM
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,
     +  NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
        LOGICAL OPENFILE
        COMMON/crpt_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,OPENFILE,NPER,
     +          DATA_HEADER
	SAVE /crpt_SPECIALS/
C
C
C  check name is valid & open file according to RUNID
	IERROR = 0
        IF ((FILNAM.NE.RUNID) .OR. (FILNAM.EQ.' ')) THEN
                CALL OPEN_FILE(RUNID,FOUND)
                IF (.NOT.FOUND) THEN
                        IERROR=1
                        RETURN
                ENDIF
        ENDIF
C  want whole file except for log section
C get this in 2 parts: first the file up to the data section+data version number
        IF (VER1.EQ.1) THEN
            start_of_data=IFORMAT(6)
        ELSE
            start_of_data=IFORMAT(7)
        ENDIF
        CALL GETSECT(1,start_of_data,IARRAY,49,IERR)
	CALL VAX_TO_LOCAL_INTS(IARRAY, start_of_data)
C       CALL GETSECT(1,start_of_data+1+(NSP1+1)*(NTC1+1)*NPER,IARRAY)
C and now the data...
        CALL GETDAT(RUNID,0,(NSP1+1)*NPER,IARRAY(start_of_data+1),
     1              LENGTH,ERRCODE)
	IF (ERRCODE .NE. 0) THEN
	   IERROR=10+ERRCODE
	   RETURN
	ENDIF
C for compressed file, log pointer points to wrong place.
C Recalculate for data version 2.
        IF (IARRAY(start_of_data).eq.2) then
                IARRAY(29)=IFORMAT(7)+1+(NSP1+1)*(NTC1+1)*NPER
        ENDIF
C we now have uncompressed data in a Version 1 format. Set data version to 1
        IARRAY(start_of_data)=1
        RETURN
        END
C-------------------------------------------------------------------------
      SUBROUTINE OPEN_DATA_FILE(RUNID,NTC,NDET,NUSE,ERRCODE)
      IMPLICIT NONE
C
C *** INPUT PARAMETERS (UNCHANGED BY THIS ROUTINE) ***
C
C RUNID - name of data file to open
C
C *** RETURNED VALUES ***
C
C NTC - number of time channels (minimum value for NTCMAX)
C NDET - the number of detectors (minimum value for NDETMAX)
C NUSE - number of user-defined UTn tables
C ERRCODE - returned error code: 0 = all OK, a l'UNIX
C                                1 = file not found or error in opening it
C
      CHARACTER*(*) RUNID
      INTEGER*4 NTC,NDET,NUSE,ERRCODE,TRUELEN
        LOGICAL FOUND
C
C common block stuff
C
        CHARACTER*256 FILNAM
	CHARACTER*256 MESS
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET1,NMON,NEFF,
     +            NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
        LOGICAL OPENFILE
C
        COMMON/CRPT_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET1,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,
     +          OPENFILE,NPER,DATA_HEADER
        SAVE /CRPT_SPECIALS/
C
      EXTERNAL OPEN_FILE, FERROR_ADD, TRUELEN
C Open file and set up CRPT_SPECIALS
      ERRCODE=0
      FOUND=.FALSE.
      IF (FILNAM .NE. RUNID) THEN
         FILNAM = ' '
         CALL OPEN_FILE(RUNID,FOUND)
         IF (.NOT. FOUND) THEN
     	    MESS = 'File '//RUNID(:TRUELEN(RUNID))//' not found'
	    CALL FERROR_ADD('OPEN_DATA_FILE', MESS, ' ')
            ERRCODE=1
            RETURN
         ENDIF
      ENDIF
      NTC=NTC1
      NDET=NDET1
      NUSE=NEFF
      RETURN
      END
C---------------------------------------------------------------------------
        SUBROUTINE OPEN_FILE(RUNID,FOUND)
        IMPLICIT NONE
C
C  Opens a RAW data file for reading. If a different fiole is already
C  open it will be closed first.
C
        CHARACTER*(*) RUNID
        LOGICAL FOUND
C
        INTEGER*4 ITEMP(3),ISTATUS,I,TRUELEN,IERR
C
        CHARACTER*256 FILNAM
	CHARACTER*256 CFILTMP
        INTEGER*4 IFILTMP(15),VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,
     +  NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
        LOGICAL OPENFILE, DAE_ACCESS, CRPT_ACCESS
        COMMON/crpt_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,
     +          OPENFILE,NPER,DATA_HEADER
        SAVE /crpt_SPECIALS/
	INTEGER VAX_TO_LOCAL_INT, ERRCODE, FASTGET_INIT
	EXTERNAL TRUELEN, VAX_TO_LOCAL_INT, VAX_TO_LOCAL_INTS,
     +		DAE_ACCESS, CRPT_ACCESS
	EQUIVALENCE (CFILTMP,IFILTMP)
C
            
C  Check that the file exists
	I = TRUELEN(RUNID)
	IF ( DAE_ACCESS(RUNID) .OR. CRPT_ACCESS(RUNID) ) THEN
            FOUND=.TRUE.
        ELSE
            INQUIRE(FILE=RUNID(1:I),EXIST=FOUND)
	ENDIF
        IF (.NOT. FOUND) THEN
C               WRITE(6,6000)RUNID
C 6000          FORMAT(' File ',A,' not found.')
                RETURN
        ENDIF
C
C  Open new file it
C
	CFILTMP = RUNID
C *** Hack - we have FILNAM in common block, so temporarily assign it for getsect_orig.f to read
	FILNAM=RUNID
	ERRCODE = FASTGET_INIT(IFILTMP, 49)
        FILNAM=' '
	IF (ERRCODE .NE. 0) THEN
	    FOUND=.FALSE.
	    RETURN
        ENDIF
C
C  Now pick out the vital parameters for future use
C
C  version number
        CALL GETSECT(21,1,ITEMP,49,IERR)
        VER1=VAX_TO_LOCAL_INT(ITEMP(1))
C--  Format section
        CALL GETSECT(22,10,IFORMAT,49,IERR)
	CALL VAX_TO_LOCAL_INTS(IFORMAT, 10)
C--  RUN section
        CALL GETSECT(IFORMAT(1),1,IVER(1),49,IERR)
	CALL VAX_TO_LOCAL_INTS(IVER(1), 1)
C--  Instrument section
        CALL GETSECT(IFORMAT(2),1,IVER(2),49,IERR)
	CALL VAX_TO_LOCAL_INTS(IVER(2), 1)
C  NDET,NMON,NEFF
        CALL GETSECT(IFORMAT(2)+67,3,ITEMP,49,IERR)
	CALL VAX_TO_LOCAL_INTS(ITEMP, 3)
        NDET=ITEMP(1)
        NMON=ITEMP(2)
        NEFF=ITEMP(3)
C--  SE section
        CALL GETSECT(IFORMAT(3),1,IVER(3),49,IERR)
	CALL VAX_TO_LOCAL_INTS(IVER(3), 1)
C  NSEP
        IF (IVER(3).EQ.1) THEN
                CALL GETSECT(IFORMAT(3)+33,1,ITEMP,49,IERR)
        ELSE
                CALL GETSECT(IFORMAT(3)+65,1,ITEMP,49,IERR)
        ENDIF
        NSEP=VAX_TO_LOCAL_INT(ITEMP(1))
C--  DAE section
        CALL GETSECT(IFORMAT(4),1,IVER(4),49,IERR)
	CALL VAX_TO_LOCAL_INTS(IVER(4), 1)
C--  TCB section
        CALL GETSECT(IFORMAT(5),1,IVER(5),49,IERR)
	CALL VAX_TO_LOCAL_INTS(IVER(5), 1)
C  NTRG
        CALL GETSECT(IFORMAT(5)+1,1,ITEMP,49,IERR)
        NTRG=VAX_TO_LOCAL_INT(ITEMP(1))
        IF (NTRG.NE.1) THEN
C               WRITE(6,6030) NTRG
C 6030  FORMAT(' This version expects only 1 time regime. NTRG=',I9)
		write(6,*) 'NTRG problem',NTRG
                FOUND=.FALSE.
                RETURN
        ENDIF
C  NPER
        CALL GETSECT(IFORMAT(5)+3,1,ITEMP,49,IERR)
        NPER=VAX_TO_LOCAL_INT(ITEMP(1))
C  NSP1,NTC1
        CALL GETSECT(IFORMAT(5)+260,2,ITEMP,49,IERR)
	CALL VAX_TO_LOCAL_INTS(ITEMP, 2)
        NSP1=ITEMP(1)
        NTC1=ITEMP(2)
C--  USER section
        IF (VER1 .EQ. 1)  THEN
		ULEN=0
		IVER(6)=0
        ELSE
                CALL GETSECT(IFORMAT(6),1,IVER(6),49,IERR)
		CALL VAX_TO_LOCAL_INTS(IVER(6), 1)
                CALL GETSECT(IFORMAT(6)+1,1,ITEMP,49,IERR)
                ULEN=VAX_TO_LOCAL_INT(ITEMP(1))
        ENDIF
C--  DATA and NOTES section
        IF (VER1 .EQ. 1) THEN
                CALL GETSECT(IFORMAT(6),1,IVER(7),49,IERR)
		IF (IFORMAT(7) .NE. 0) THEN
			CALL GETSECT(IFORMAT(7),1,IVER(8),49,IERR)
			CALL VAX_TO_LOCAL_INTS(IVER(8), 1)
		ELSE
			IVER(8) = 0
		ENDIF
        ELSE
                CALL GETSECT(IFORMAT(7),1,IVER(7),49,IERR)
                IF (IFORMAT(8) .NE. 0) THEN
			CALL GETSECT(IFORMAT(8),1,IVER(8),49,IERR)
			CALL VAX_TO_LOCAL_INTS(IVER(8), 1)
		ELSE
			IVER(8) = 0
		ENDIF
        ENDIF
	CALL VAX_TO_LOCAL_INTS(IVER(7), 1)

C DATA section header
        do i=1,32
                data_header(i)=0
        enddo   
        if ( (IVER(7) .GE. 2) .AND. (.NOT. DAE_ACCESS(RUNID)) ) then
                CALL GETSECT(IFORMAT(7)+1,32,data_header,49,IERR)
		CALL VAX_TO_LOCAL_INTS(data_header, 32)
        endif
C  finally store the file name and set flags to true
        FILNAM=RUNID
        FOUND=.TRUE.
C       WRITE(6, *)' FILE OPENED'
        RETURN
        END
C
      SUBROUTINE READ_DATA(RUNID,ERRCODE,ISPEC,DELT_WORK,SPEC_WORK,
     +                     TTHE_WORK,
     +                     L2_WORK,NDETMAX,TCB,IDAT,NTCMAX,L1,L2,
     +                     TTHE,DELT,PHI,RUN_TITLE,DURATION,
     +                     COMBINED_TIME,USER_NAME,INST_NAME,RUN_NO,
     +			   USER_TABLES,NUSE,QUICK)
C
C *** THIS SHOULD ONLY BE CALLED AFTER A CALL TO OPEN_DATA_FILE
C
C *** if (QUICK .EQ. 1)  only get spectrun data and not titles etc
C
C *** INPUT PARAMETERS (UNCHANGED BY THIS ROUTINE) ***
C
C RUNID
C ISPEC - the identifier of the spectrum to retrieve 
C NDETMAX - the size of the work arrays DELT_WORK, SPEC_WORK, 
C                                     L2_WORK and TTHE_WORK
C NTCMAX - maximum size of arrays IDAT and TCB
C
C *** RETURNED VALUES ***
C
C L1 - Primary flight path
C L2 - Secondary flight path (averaged over detectors)
C TTHE - two theta (averaged over detectors)
C DELT -
C IDAT - returned spectra file (of length NTC1+1)
C TCB - array of NTC1+1 time channel boundaries
C USER_TABLES - contains the UTn tables (NUSE of them) averaged over detectors
C ERRCODE - returned error code:   0 = All OK
C                                  1 = file not open
C                                  2 = NTCMAX, NDETMAX or NUSE inconsistent
C                                  3 = ISPEC out of range
C                                  4 = other read error
C				   5 = NUSE out of range
C
C
C *** ADDITIONAL PARAMETERS ***
C
C REQUIRES WORK ARRAYS SPEC_WORK(NDETMAX), DELT_WORK(NDETMAX), 
C  TTHE_WORK(NDETMAX) AND L2_WORK(NDETMAX)
C
      IMPLICIT NONE
      CHARACTER*(*) RUNID
      INTEGER*4 QUICK,NDETMAX,NTCMAX,IDAT(NTCMAX),ERRCODE,ISPEC,
     1       NUSE,IERR
      INTEGER*4 I,J,SPEC_WORK(NDETMAX),NDETMATCH
      REAL*4 DELT_WORK(NDETMAX),L2_WORK(NDETMAX),TTHE_WORK(NDETMAX),
     1       USER_TABLES(NUSE)
C
        REAL*4 TCB(NTCMAX)
        REAL*4 TTHE,L2,DELT,PHI,RNDET,DURATION
        REAL*4 L1,RVPBWK(64)
        INTEGER*4 IVPBWK(64),RPB(32),LENGTH_OUT
	INTEGER TRUELEN
	CHARACTER*8 ERROR1, ERROR2
	CHARACTER*2 UTNUM
	CHARACTER*5 RUN_NO
        CHARACTER*80 HEADER
        CHARACTER*20 USER_NAME
        CHARACTER*8 RUN_IDENT,INST_NAME
        CHARACTER*24 EXPT_TITLE
        CHARACTER*12 START_DATE
        CHARACTER*8 START_TIME,RUN_DURATION
	CHARACTER*21 COMBINED_TIME
        CHARACTER*80 RUN_TITLE
C
        CHARACTER*256 FILNAM
	CHARACTER*256 MESS
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,
     +            NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
C *** CONV_ERR controls if an error occurs duing VAXF_TO_LOCAL
	INTEGER CONV_ERR,IDUM
        LOGICAL OPENFILE
        EQUIVALENCE (IVPBWK,RVPBWK)
C
        COMMON/CRPT_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,
     +          OPENFILE,NPER,DATA_HEADER
        SAVE /CRPT_SPECIALS/
C
      EXTERNAL GETDAT,GETPARR,GETPARI,GETPARC,VAXF_TO_LOCAL,
     1         TRUELEN,FERROR_ADD
      CONV_ERR = 0
      ERRCODE = 0
      IF ((RUNID .NE. FILNAM) .OR. (FILNAM .EQ. ' ')) THEN
	 CALL FERROR_ADD('READ_DATA', 
     1                   'Error in file specification', ' ')
         ERRCODE=1
         RETURN
      ENDIF
      IF (NTC1+1 .GT. NTCMAX) THEN
	 WRITE(ERROR1, '(I8)') NTC1
	 WRITE(ERROR2, '(I8)') NTCMAX
         MESS = 'Too many time channels: NTC1 = '//ERROR1//
     1		', NTCMAX = '//ERROR2
	 CALL FERROR_ADD('READ_DATA', MESS, ' ')
         ERRCODE=2
         RETURN
      ENDIF

      IF ((ISPEC .GT. ((NSP1+1)*NPER)-1) .OR. (ISPEC .LT. 0)) THEN
	 WRITE(MESS, 175) ISPEC, ((NSP1+1)*NPER)-1
 175	 FORMAT('Invalid spectrum number = ',I5,
     1   ' (spectra must be in the range 0 - ',I5,')')
	 CALL FERROR_ADD('READ_DATA', MESS, ' ')
         ERRCODE=3
         RETURN
      ENDIF
      IF (NDET .GT. NDETMAX) THEN
	 WRITE(ERROR1, '(I8)') NDET
	 WRITE(ERROR2, '(I8)') NDETMAX
	 MESS = 'Too many detectors: '//ERROR1//' > '//ERROR2
	 CALL FERROR_ADD('READ_DATA', MESS, ' ')
         ERRCODE=2
         RETURN
      ENDIF
      IF (NUSE .NE. NEFF) THEN
	 WRITE(ERROR1, '(I8)') NUSE
	 WRITE(ERROR2, '(I8)') NEFF
         MESS = 'Invalid number of user parameters: '//
     1		ERROR1//' .NE. '//ERROR2
	 CALL FERROR_ADD('READ_DATA', MESS, ' ')
         ERRCODE=2
         RETURN
      ENDIF
      IF (QUICK .EQ. 0) THEN
      CALL GETPARC(RUNID,'HDR',HEADER,1,LENGTH_OUT,IERR)
      IF (IERR .NE. 0) GOTO 999
C Run identifier e.g. LAD12345
      RUN_IDENT=HEADER(1:8)
      RUN_NO=HEADER(4:8)
C User Name
      USER_NAME=HEADER(9:28)
C Experiment short title
      EXPT_TITLE=HEADER(29:52)
C Start date
      START_DATE=HEADER(53:64)
C Start Time
      START_TIME=HEADER(65:72)
      COMBINED_TIME = START_DATE(1:TRUELEN(START_DATE))//' '//
     1                START_TIME(1:TRUELEN(START_TIME))
C
C      CALL GETPARC(RUNID,'TITL',RUN_TITLE,1,LENGTH_OUT,IERR)
C      IF (IERR .NE. 0) GOTO 999
      RUN_TITLE=EXPT_TITLE
C
      CALL GETPARC(RUNID,'NAME',INST_NAME,1,LENGTH_OUT,IERR)
      IF (IERR .NE. 0) GOTO 999
C time channel boundaries
      CALL GETPARR(RUNID,'TCB1',TCB,NTCMAX,IDUM,IERR)
      IF (IERR .NE. 0) GOTO 999
C do not need VAXF_TO_LOCAL as stored as integers
      ENDIF
C *** end of >>> if ( quick=0) <<<    
      CALL GETPARR(RUNID,'RVPB',RVPBWK,64,IDUM,IERR)
      IF (IERR .NE. 0) GOTO 999
C obtained via equivalence statement
      L1=RVPBWK(23)
C spectrum table
      CALL GETPARI(RUNID,'SPEC',SPEC_WORK,NDETMAX,IDUM,IERR)
      IF (IERR .NE. 0) GOTO 999
C run duration (s)
      CALL GETPARI(RUNID, 'RPB', RPB, 32, IDUM,IERR)
      IF (IERR .NE. 0) GOTO 999
      DURATION=FLOAT(RPB(13))
C
      CALL GETPARR(RUNID,'DELT',DELT_WORK,NDETMAX,IDUM,IERR)
      IF (IERR .NE. 0) GOTO 999
C      CALL VAXF_TO_LOCAL(DELT_WORK,NDET,IERR)
C      IF (IERR .NE. 0) CONV_ERR = 1
      CALL GETPARR(RUNID,'LEN2',L2_WORK,NDETMAX,IDUM,IERR)
      IF (IERR .NE. 0) GOTO 999
C      CALL VAXF_TO_LOCAL(L2_WORK,NDET,IERR)
C      IF (IERR .NE. 0) CONV_ERR = 1
      CALL GETPARR(RUNID,'TTHE',TTHE_WORK,NDETMAX,IDUM,IERR)
      IF (IERR .NE. 0) GOTO 999
C      CALL VAXF_TO_LOCAL(TTHE_WORK,NDET,IERR)
C      IF (IERR .NE. 0) CONV_ERR = 1
C
C average TTHE, L2 and DELT over detectors used for given spectrum
C
      TTHE=0.0
      L2=0.0             
      DELT=0.0
      J = 0
      DO I=1,NDET
         IF (SPEC_WORK(I) .EQ. ISPEC) THEN
            TTHE=TTHE+TTHE_WORK(I)
            DELT=DELT+DELT_WORK(I)
            L2=L2+L2_WORK(I)
	    J = J + 1
         ENDIF
      ENDDO   
      NDETMATCH=J
C *** PHI ***
      CALL GETPARR(RUNID, 'PHI', DELT_WORK, NDETMAX, IDUM, IERR)
      PHI=0.0
      IF (IERR .EQ. 0) THEN
          DO I=1,NDET
             IF (SPEC_WORK(I) .EQ. ISPEC) THEN
                PHI=PHI+DELT_WORK(I)
             ENDIF
          ENDDO
      ENDIF
C
      DO I=1,NUSE
         USER_TABLES(I)=0.0
      ENDDO
      DO I=1,NUSE
         WRITE(UTNUM, '(I2.2)') I
         CALL GETPARR(RUNID, 'UT'//UTNUM, DELT_WORK, NDETMAX,IDUM, IERR)
         IF (IERR .NE. 0) GOTO 999
C         CALL VAXF_TO_LOCAL(DELT_WORK, NDET, IERR)
C         IF (IERR .NE. 0) CONV_ERR = 1
         DO J=1,NDET
            IF (SPEC_WORK(J) .EQ. ISPEC) THEN
               USER_TABLES(I)=USER_TABLES(I)+DELT_WORK(J)
            ENDIF
         ENDDO
      ENDDO
      IF (NDETMATCH .NE. 0) THEN
         RNDET=FLOAT(NDETMATCH)
         TTHE=TTHE/RNDET
         DELT=DELT/RNDET
         L2=L2/RNDET
         PHI=PHI/RNDET
         DO I=1,NUSE
            USER_TABLES(I)=USER_TABLES(I)/RNDET
         ENDDO
      ENDIF

      CALL GETDAT(RUNID,ISPEC,1,IDAT,NTCMAX,IERR)
      IF (IERR .NE. 0) GOTO 999
C
      IF (CONV_ERR .NE. 0) THEN
          CALL FERROR_ADD('INFORMATION', 
     1                    'Error during convertion VAX format to IEEE', 
     2                    'May be unimportant - check data')
      ENDIF
      RETURN
999   ERRCODE=4
      CALL FERROR_ADD('READ_DATA', 'Some other error', ' ')
      RETURN
      END
C
	INTEGER FUNCTION TRUELEN(STRING)
C
C *** return string length when trailing blanks are discarded
C
	IMPLICIT NONE
	CHARACTER*(*) STRING
	INTEGER I
	I=LEN(STRING)
 10	IF (STRING(I:I) .EQ. ' ') THEN
	   I=I-1
	   IF (I .GT. 0) GOTO 10
	ENDIF
	TRUELEN=I
	RETURN
	END
C
	SUBROUTINE DEFAULT_FILE_NAME(IUNIT, FILE)
C
C *** this routine will find the default name for unit IUNIT
C
	IMPLICIT NONE
	INTEGER*4 IUNIT
	CHARACTER*(*) FILE
c	OPEN(UNIT=IUNIT, STATUS='NEW')
	OPEN(UNIT=IUNIT, STATUS='UNKNOWN')
	INQUIRE(UNIT=IUNIT, NAME=FILE)
	CLOSE(UNIT=IUNIT, STATUS='DELETE')
	RETURN
	END
C
	SUBROUTINE FORT_FILE(IUNIT, WORK)
	BYTE WORK(256)
	CHARACTER*256 FILE_NAME
	INTEGER I, STAT, IUNIT, TRUELEN
	EXTERNAL TRUELEN
	FILE_NAME = ' '
	INQUIRE(UNIT=IUNIT, NAME=FILE_NAME, IOSTAT=STAT)
	DO I=1,TRUELEN(FILE_NAME)
	    WORK(I) = ICHAR(FILE_NAME(I:I))
	ENDDO
	DO I=TRUELEN(FILE_NAME)+1,256
	    WORK(I) = 0
	ENDDO
	RETURN
	END
C
        SUBROUTINE BYTE_REL_EXPN(INDATA,NIN,NFROM,OUTDATA,NOUT,ISTATUS)
        IMPLICIT NONE
C
C Expansion of byte-relative format into 32bit format
C
C Each integer is stored relative to the previous value in byte form.The first
C is relative to zero. This allows for numbers to be within + or - 127 of
C the previous value. Where a 32bit integer cannot be expressed in this way a 
C special byte code is used (-128) and the full 32bit integer stored elsewhere.
C The final space used is (NIN-1)/4 +1 + NEXTRA longwords, where NEXTRA is the
C number of extra longwords used in giving absolute values.
C
C
C Type definitions
C   Passed parameters
        INTEGER NIN
        BYTE    INDATA(NIN)     
C                               Array of NIN compressed bytes
        INTEGER NFROM           
C                               pass back data from this 32bit word onwards
C
        INTEGER NOUT            
C                               Number of 32bit words expected
        INTEGER OUTDATA(NOUT)   
C                               outgoing expanded data
        INTEGER ISTATUS         
C                               Status return
C                                      =1  no problems!
C                                      =3  NOUT .lt.NIN/5
C                                      =2  NIN .le.0
C                                      =4  NOUT .gt.NIN
C                                      =6  number of channels lt NOUT
C   Local definitions
        INTEGER I,J,ITEMP
        BYTE BTEMP(4)
        EQUIVALENCE (ITEMP,BTEMP)       
C                               For byte UNpacking
c
c
C Assume innocent until proven guilty
        ISTATUS=1
c First check no slip-ups in the input parameters
        IF (NIN.LE.0) THEN
                ISTATUS=2
                GOTO 100
        ENDIF
        IF (NOUT+NFROM-1.gt.NIN) THEN
                ISTATUS=4
                GOTO 100
        ENDIF
C
C Set initial absolute value to zero and channel counter to zero
        ITEMP=0
        J=0
C
C Loop over all expected 32bit integers
        DO I=1,NFROM+NOUT-1
            J=J+1
            IF (J.GT.NIN) THEN  
C check there are enough bytes
                ISTATUS=6
                GOTO 100
            ENDIF
C
C if number is contained in a byte
            IF (INDATA(J).NE.-128) THEN
C add in offset to base
                ITEMP=ITEMP+INDATA(J)   
            ELSE
C Else skip marker and pick up new absolute value
                IF (J+4.GT.NIN) THEN
C check there are enough bytes
                    ISTATUS=6
                    GOTO 100
                ENDIF
C unpack 4 bytes
                BTEMP(1)=INDATA(J+1)    
                BTEMP(2)=INDATA(J+2)
                BTEMP(3)=INDATA(J+3)
                BTEMP(4)=INDATA(J+4)
		CALL VAX_TO_LOCAL_INTS(ITEMP, 1)
                J=J+4
            ENDIF
C update current value
            IF (I.GE.NFROM) THEN
                OUTDATA(I-NFROM+1)=ITEMP
            ENDIF
        ENDDO
c
c
C expansion OK, but excessive number of bytes given to the routine
 100    IF (NOUT.lt.NIN/5) ISTATUS=3
C       TYPE *,'J=',J,' I=',I
        RETURN
        END
C---------------------------------------------------------------------
        SUBROUTINE GETPARR(RUNID,NAME,RVALUE,LENGTH_IN,
     1			   LENGTH_OUT,ERRCODE)
        IMPLICIT NONE
C
C *** ERROR CODES ***
C
C    0 = all OK
C    1 = file RUNID not currently open
C    2 = file number is '00000' and so no access available
C    3 = asked for non-existent parameter
C    4 = other error
C

C  Gets REAL parameter(s) from RAW data file
C
        CHARACTER*(*) RUNID
        CHARACTER*(*) NAME
        INTEGER*4 ERRCODE,LENGTH_IN,IERR
	INTEGER*4 LENGTH_OUT,LENGTH
	INTEGER*4 RVALUE(LENGTH_IN)
C
C
        INTEGER*4 IWORK(256),ISTORE(64)
        REAL*4 WORK(256)
        COMMON/crpt_WORK/IWORK
        SAVE /crpt_WORK/
        EQUIVALENCE (IWORK,WORK)
C
        EQUIVALENCE (TEMP,ITEMP)
C
        INTEGER*4 ITABLE,IPRE1,I,ITEMP,VAX_TO_LOCAL_INT
        REAL*4 TEMP,EXTRA
C
        CHARACTER*256 FILNAM
	CHARACTER*256 MESS
        INTEGER*4 VER1,IFORMAT(10),IVER(10),NDET,NMON,NEFF,
     +  NSEP,NTRG,NSP1,NTC1,ULEN,NPER,DATA_HEADER(32)
        LOGICAL OPENFILE
        COMMON/crpt_SPECIALS/ FILNAM,VER1,IFORMAT,IVER,
     +          NDET,NMON,NEFF,NSEP,NTRG,NSP1,NTC1,ULEN,OPENFILE,NPER,
     +          DATA_HEADER
        SAVE /crpt_SPECIALS/
	EXTERNAL VAX_TO_LOCAL_INT,VAX_TO_LOCAL_INTS
        ERRCODE=0
        IERR=0
	LENGTH_OUT=0
	LENGTH=LENGTH_IN
C
C
C  decide whether it's CRPT or just a f
        IF (RUNID(4:8).EQ.'00000') THEN
C               WRITE(6,6000)
C 6000          FORMAT(' Accessing of current run not available')
C               GOTO 999
	   	CALL FERROR_ADD('GETPARR',
     1     	'Runid is 00000',
     2     	' ')
           ERRCODE=2
           RETURN
        ENDIF
C
C  check name is valid & open file according to RUNID
        IF ((FILNAM .NE. RUNID) .OR. (FILNAM .EQ. ' ')) THEN
           MESS = 'Access requested to '//RUNID//
     1		  ' when '//FILNAM//' open'
	   CALL FERROR_ADD('GETPARR',
     1     MESS,
     2     ' ')
           ERRCODE=1
           RETURN
        ENDIF
C
C  read variables into RVALUE
C ----------------------------
C
C  Instrument section
        ITABLE=IFORMAT(2)+70+2*NMON
        IF (NAME.EQ.'LEN2') THEN
                CALL GETSECT(ITABLE+2*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'OMEG') THEN
          IF (IVER(2).EQ.1) THEN
                CALL GETSECT(ITABLE+3*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
	  ELSE
          	MESS = 'Cannot access '//NAME
	   	CALL FERROR_ADD('GETPARR',
     1     	MESS,
     2     	' ')
	        ERRCODE=4
          ENDIF
        ELSEIF (NAME.EQ.'DELT') THEN
          IF (IVER(2).NE.1) THEN
                CALL GETSECT(ITABLE+NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
	  ELSE
          	MESS = 'Cannot access '//NAME
	   	CALL FERROR_ADD('GETPARR',
     1     	MESS,
     2     	' ')
	        ERRCODE=4
          ENDIF
        ELSEIF (NAME.EQ.'TTHE') THEN
                CALL GETSECT(ITABLE+4*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'PHI') THEN
          IF (IVER(2).EQ.1) THEN
                CALL GETSECT(ITABLE+5*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
	  ELSE
          	MESS = 'GETPARR: No item '//NAME//' in RAW file'
C *** We do not flag an error here as then 'GET:SPECTRUM' would return an error vie READ_DATA
C *** comment out message until i find a RAW file that does have it set!
C	   	CALL FERROR_ADD('INFORMATION', MESS, ' ')
	        ERRCODE=4
          ENDIF
        ELSEIF (NAME.EQ.'EF01') THEN
                CALL GETSECT(ITABLE+6*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'EF02') THEN
                CALL GETSECT(ITABLE+7*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'EF03') THEN
                CALL GETSECT(ITABLE+8*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'EF04') THEN
                CALL GETSECT(ITABLE+9*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'EF05') THEN
                CALL GETSECT(ITABLE+10*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'EF06') THEN
                CALL GETSECT(ITABLE+11*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'RVPB') THEN
                CALL GETSECT(IFORMAT(2)+3,64,RVALUE,49,IERR)
		LENGTH = 64
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
	ELSEIF (NAME(1:2) .EQ. 'UT') THEN
		READ(NAME(3:4),'(I2)') I
                CALL GETSECT(ITABLE+(4+I)*NDET,NDET,RVALUE,49,IERR)
		LENGTH=NDET
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
C  Time channel boundaries section
C     time channel area definition
        ELSEIF (NAME.EQ.'TCP1') THEN
                CALL GETSECT(IFORMAT(5)+267,20,RVALUE,49,IERR)
		LENGTH=20
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
        ELSEIF (NAME.EQ.'TCB1') THEN
		IF (VER1 .NE. 1) THEN
                    CALL GETSECT(IFORMAT(4)+1,64,ISTORE,49,IERR)
		    CALL VAX_TO_LOCAL_INTS(ISTORE, 64)
		    EXTRA=FLOAT(ISTORE(24)*4)
		ELSE
		    EXTRA=0.0
		ENDIF
C  if tcb's requested in real form then return as microsecs
                CALL GETSECT(IFORMAT(5)+287,1,ISTORE,49,IERR)
                IPRE1=VAX_TO_LOCAL_INT(ISTORE(1))
                CALL GETSECT(IFORMAT(5)+288,NTC1+1,RVALUE,49,IERR)
		CALL VAX_TO_LOCAL_INTS(RVALUE, NTC1+1)
C  conversion loop - from clock pulses to microsecs
                DO 100 I=1,NTC1+1
                   TEMP=FLOAT(RVALUE(I))/32.0*FLOAT(IPRE1) + EXTRA
                   RVALUE(I)=ITEMP
 100            CONTINUE
		LENGTH=NTC1+1
	ELSEIF (NAME .EQ. 'DAT1') THEN
		IF (VER1 .NE. 1) THEN
		    CALL GETSECT(IFORMAT(6)+2, ULEN, RVALUE, 49, IERR)
		    LENGTH = ULEN
		    CALL VAXF_TO_LOCAL(RVALUE, LENGTH, IERR)
		ELSE
	    	    CALL FERROR_ADD('GETPARR', 
     1      		    'No user section in this file', 
     2			    ' ')
		ENDIF
C  non existant requests
        ELSEIF (NAME .EQ. 'RRPB') THEN
                CALL GETSECT(IFORMAT(1)+62,32,RVALUE,49,IERR)
		LENGTH = 32
		CALL VAXF_TO_LOCAL(RVALUE, LENGTH,IERR)
        ELSE
C     	        MESS = 'No such REAL parameter as '//NAME
C	        CALL FERROR_ADD('GETPARR', 
C     1      		    MESS,
C     2			    ' ')
                ERRCODE=3
		LENGTH = 0
                RETURN
        ENDIF
C
	IF ((IERR .NE. 0) .AND. (ERRCODE .EQ. 0)) THEN
       	    MESS = 'Error in reading data from file '//RUNID
	    CALL FERROR_ADD('GETPARR', 
     1			    MESS,
     2			    ' ')
	    ERRCODE = 4
	ENDIF
C
	LENGTH_OUT=LENGTH
 999    RETURN
        END
C
	LOGICAL FUNCTION DAE_ACCESS(NAME)
	IMPLICIT NONE
	CHARACTER*(*) NAME
	LOGICAL DCE_NAME
	INTEGER I, TRUELEN
	EXTERNAL TRUELEN
	I = TRUELEN(NAME)
	DCE_NAME = (NAME(1:4) .EQ. '/.:/') .OR. (NAME(1:5) .EQ. '/.../')
	IF ( DCE_NAME .AND. (NAME(I-3:I) .EQ. '_dae') ) THEN
	    DAE_ACCESS = .TRUE.
	ELSE
	    DAE_ACCESS = .FALSE.
	ENDIF
	RETURN
	END
C	
	LOGICAL FUNCTION CRPT_ACCESS(NAME)
	IMPLICIT NONE
	CHARACTER*(*) NAME
	LOGICAL DCE_NAME
	INTEGER I, TRUELEN
	EXTERNAL TRUELEN
	I = TRUELEN(NAME)
	DCE_NAME = (NAME(1:4) .EQ. '/.:/') .OR. (NAME(1:5) .EQ. '/.../')
	IF ( DCE_NAME .AND. (NAME(I-4:I) .EQ. '_crpt') ) THEN
	    CRPT_ACCESS = .TRUE.
	ELSE IF ( (INDEX(NAME, 'CURRENT.RUN') .NE. 0) .OR. 
     +            (INDEX(NAME, 'current.run') .NE. 0) ) THEN
	    CRPT_ACCESS = .TRUE.
	ELSE
	    CRPT_ACCESS = .FALSE.
	ENDIF
	RETURN
	END
C
C---------------------------------------------------------------------
* DEC/CMS REPLACEMENT HISTORY, Element BYTE_REL_EXPN.FOR
* *1    11-AUG-1989 09:28:08 KJK "New routine coming in at 2.4.5"
* DEC/CMS REPLACEMENT HISTORY, Element BYTE_REL_EXPN.FOR
