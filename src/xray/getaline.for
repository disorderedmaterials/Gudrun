	subroutine getaline(unit)

	include 'inputfilestrings.inc'
	integer*4 unit

c Procedure to read and parse character strings, using the spc2 string for between
c words and spc5 to indicate end of line.

c Find the next valid line

	read(unit,'(a)',iostat=ierr) line
	nwords=0
	do while (nwords.eq.0.and.ierr.eq.0)

c Parse line into words

	   call parse(line,nwords,ncf,ncl,nchartext,nwordsdim
     1,spc2,spc5,lenspc2,lenspc5,ierrline)
	   if(nwords.eq.0) read(unit,'(a)',iostat=ierr) line
	end do

	return
	end
