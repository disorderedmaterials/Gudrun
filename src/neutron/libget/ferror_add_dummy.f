C
C       Module:         FERROR_ADD_DUMMY
C       Author:         Freddie Akeroyd, ISIS
C       Purpose:        FORTRAN interface to dummy Genie error handling interface
C
C	$Id: ferror_add_dummy.f,v 1.1 1997/01/01 14:47:00 faa Exp $
C
	SUBROUTINE FERROR_ADD(OBJECT, FAIL, SOLUTION)
	IMPLICIT NONE
	CHARACTER*(*) OBJECT, FAIL, SOLUTION
	WRITE(6,*) 'Error returned from ',OBJECT
	WRITE(6,*) FAIL
	IF (SOLUTION .NE. ' ') THEN
    	    WRITE(6,*) 'Solution: ',SOLUTION
	ENDIF
	RETURN
	END
C
