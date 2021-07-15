c
c variables to define the strings used to interpret the Gudrun input file
c
	integer*4 nworddim,nwordsdim
	parameter (nworddim=256,nwordsdim=20)
	byte test1(nworddim)
	character*256 line		!title of block of data
	equivalence(line,test1)
	character*20 spc2,spc5,pathseparator   !spacing strings between values and between last value and comment
	integer*4 nchar,index1,index2
	integer*4 lenspc2,lenspc5	!no. of characters in respective strings
	integer*4 nwords,ncf(nwordsdim),ncl(nwordsdim),nchartext,ierrline
	common/InputFileStrings_c/line,spc2,spc5,pathseparator
	common/InputFileStrings_i/nchar,index1,index2,lenspc2,lenspc5
     1,nwords,ncf,ncl,nchartext,ierrline
