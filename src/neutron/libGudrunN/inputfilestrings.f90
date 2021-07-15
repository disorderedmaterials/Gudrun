!     
! File:   inputfilestrings.mod
! Author: aks45
!
! Created on 08 March 2012, 13:13
!

MODULE inputfilestrings

! variables to define the strings used to interpret the Gudrun input file

    implicit none

    integer, parameter :: spcdim=20,nwordsdim=250
    character(len=spcdim) :: spc2=' ',spc5='          '   !spacing strings between values and between last value and comment
    integer :: nchar,index1,index2
    integer :: lenspc2=1,lenspc5=10              !no. of characters in respective strings
    integer :: nwords
    integer :: nchartext,ierrline
    integer, dimension(nwordsdim) :: ncf,ncl
    character(len=256) :: line
    public :: parse,nwords,ncf,ncl


 CONTAINS

    subroutine parse(linestring)

!c parses the text in line to separate words with the assumed delimiter spc2
!c and spc5 is the string between the last value and any comments in the line.
!c ncf and ncl contain the first and last character numbers of each word
!c ncharword is the maximum number of characters in a word, nworddim is dimension of the
!c word array, nspcword is the number of spaces between words needed to decide that the
!c end of text has been reached.
!Adapted 6/3/2014 to accommodate different line endings

        character(len=256), intent(in), optional :: linestring
        integer :: i

        nwords=0

        if(present(linestring)) then
            line=linestring
        end if

!c step through the characters and identify words (separated by spc2 strings). Extra spaces
!c are ignored unless the number of spaces exceeds nspcword, in which case the search
!c for words terminates
!c Get number of characters up to end of last valid variable

        nchartext=index(line,spc5(1:lenspc5))-1
        if(nchartext.lt.1) then !If spc5 is not found, then search for a !: ignore anything beyond this
            nchartext=index(line,'!')-1 
            if(nchartext.lt.1) then !If '!' is not found, assume this line contains only data - no comment is appended
                nchartext=len_trim(line)
            end if
        end if
        if(nchartext.lt.1) then ! If there is nothing on the line will report an error
            ierrline=1
            return
        endif
!        write(6,*) nchartext

!c Search for words in the current string

        i=1
        do while(nwords.lt.nwordsdim.and.i.le.nchartext)

!c Ignore empty space before a word

            do while(line(i:i).eq.' '.and.i.le.nchartext)
                i=i+1
            end do
            if(i.le.nchartext) then
                index1=index(line(i:nchartext),spc2(1:lenspc2))-1
                if(index1.gt.0) then
                    nwords=nwords+1
                    ncf(nwords)=i
                    ncl(nwords)=i+index1-1
                    i=ncl(nwords)+1
                else
!c If index1.lt.1 it can only mean this is the last word of the set
                    nwords=nwords+1
                    ncf(nwords)=i
                    ncl(nwords)=nchartext
                    i=ncl(nwords)+1
                endif
            endif
        end do
        ierrline=0
!        write(6,*) nwords,ncf(1),ncl(nwords)

    end subroutine parse

    function get_ext(filename)
        !Searches for the last '.' of filename, and uses this to create the extension
        character(len=*), intent(in)    :: filename
        character(len=256)              :: get_ext
        integer :: lenfilename,i
        logical :: found
        lenfilename=len_trim(filename)
        i=lenfilename
        found = filename(i:i)=='.'
        do while (i.gt.0.and..not.found)
            i=i-1
            found = filename(i:i)=='.'
        end do
        if(i.gt.1) then
            get_ext=filename(i+1:lenfilename)
        else
            get_ext=' '
        end if
    end function get_ext

    subroutine change_ext(filename,newext)
        !Searches for the last '.' of filename, and replaces with newext from that point on
        character(len=*), intent(inout)     :: filename
        character(len=*), intent(in)        :: newext
        integer                             :: lentotal,lenlastdot
        integer                             :: lennewext
        logical                             :: found

        lentotal=len_trim(filename)
        lennewext=len_trim(newext)
        if(lentotal.gt.0) then
            lenlastdot=lentotal
            !Does this filename have a dot?
            if(index(filename,'.').gt.0) then
                !If so remove everything after and including the last dot
                lenlastdot=lentotal
                do while (index(filename(lenlastdot:lenlastdot),'.').eq.0)
                    lenlastdot=lenlastdot-1
                end do
                !Don't include the dot!
                lenlastdot=lenlastdot-1
            end if
            filename=filename(1:lenlastdot)//'.'//newext(1:lennewext)
            return
        else
            filename='junk.'//newext(1:lennewext)
            return
        end if

        return

    end subroutine change_ext

    subroutine blocklines(linestrings,linecountin,nlinesblock,nwordsmax)
        !Determine the number of lines until the specified end block marker is found
        character(len=*), dimension(:), intent(in) :: linestrings
        integer, intent(in) :: linecountin
        integer, intent(out) :: nlinesblock,nwordsmax
        integer :: i,j
        logical :: found
        i=0
        j=linecountin
        nwordsmax=0
        found=.false.
        !Make sure we don't step off the end of the array, nor go beyond the current group marker('}').
        do while (j < size(linestrings).and..not.found )
            j=j+1
            found=index(linestrings(j),'}') > 0
            if(.not.found) then
                i=i+1
                call parse(linestrings(j))
                if (nwords.gt.nwordsmax) nwordsmax=nwords
            end if
        end do
        nlinesblock=i
        write(6,'(a,1x,i4,1x,a)') 'blocklines> found',nlinesblock,'lines.'
    end subroutine blocklines
    
    subroutine getaline(unit)

        integer, intent(in) :: unit
        integer :: ierr

! Procedure to read and parse character strings, using the spc2 string for between
! words and spc5 to indicate end of line.

! Find the next valid line

        read(unit,'(a)',iostat=ierr) line
        nwords=0
        do while (nwords.eq.0.and.ierr.eq.0)
! Parse line into words
            call parse()
            if(nwords.eq.0) read(unit,'(a)',iostat=ierr) line
        end do
    return
    end subroutine getaline
    
END MODULE inputfilestrings
