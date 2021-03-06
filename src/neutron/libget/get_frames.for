	PROGRAM get_frames

      CHARACTER*80 runid
	INTEGER good_frames, total_frames
	INTEGER rpb1(32), length_out, errcode
	REAL rpb2(32), total_uamps, good_uamps, perc_bad
	LOGICAL found

	WRITE(6,*) 'Input data file:'
	READ(5,'(A)') runid

C Open the RAW file
	CALL OPEN_FILE(runid,found)
	
	IF (.NOT.found) STOP 'File not found!'
	
C Get the array of integers from the RPB block

	CALL GETPARI(runid,'RPB',rpb1,32,length_out,errcode)
	good_frames = rpb1(10)
	total_frames = rpb1(11)

C Now get RPB Block as reals in order to extract the uAmps
C (for this case the block to get is RRPB)
	
	CALL GETPARR(runid,'RRPB',rpb2,32,length_out,errcode)
      good_uamps = rpb2(8)
      total_uamps = rpb2(9)

	perc_bad = (real(total_frames)-real(good_frames))/real(total_frames)
 	perc_bad = perc_bad * 100.

	WRITE(6,*) ' Total No. of frames : ',total_frames
	WRITE(6,*) '         Good frames : ',good_frames
	WRITE(6,*) ' % of bad frames     : ',perc_bad
	WRITE(6,*)
	WRITE(6,*) ' Total Beam current  :',total_uamps
	WRITE(6,*) '  Good Beam current  :',good_uamps

	END
              