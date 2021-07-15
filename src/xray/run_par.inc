c
c global run variables common to all routines
c
	character*80 info		!full title of run
	character*80 headinfo		!extra header info
	character*256 directinp		!directory where Gudrun is being run
	character*256 direct		!name of directory containing raw files
	character*80 inst		!instrument name
	character*3 ext			!raw file extension
	character*256 deadtimefile	!name of file containing detector deadtimes
	character*3 pathseparator	!string used to described the character(s) between directory
c names

	integer*4 nspec			!no. of spectra in 1 period (not including spectrum 0)		
	integer*4 nspecb		!no. of spectar in 1 period including spectrum 0 if present
	integer*4 nproc			!no. of spectra to process at one time
	integer*4 NCHAN,nchanv,nchans	!no. of time channels (TC's)
	integer*4 nchanb,nchanbv,nchanbs !no. of TC boundaries ( = nchan+1)
	integer*4 nspch,nspchv,nspchs	 !product of no. spec. times no. of tcb's
	integer*4 nper			!no. of periods in this file
	integer*4 ndet			!no. of detectors in this file
	integer*4 nmod			!number of detector modules (DAE inputs)
	integer*4 specno(mdet)		!spectrum numbers of detectors
	integer*4 detno(mspec)		!detector numbers for spectra
	integer*4 udet(mdet)		!detector labels
	integer*4 crateno(mmod)		!crate numbers of modules
	integer*4 slotno(mmod)		!slot numbers of modules
	integer*4 ndetpermodule(mmod)	!no. of detectors in a module
	integer*4 crate(mdet)		!crate number for each detector
	integer*4 slot(mdet)		!slot number of each detector
	integer*4 module(mdet)		!module number of each detector
	integer*4 lenpathseparator	!number of characters in path separator
	integer*4 lendirect,lendirectinp	!length in characters of direct and directinp

	real*4 tcb(mchan),tcbv(mchan),tcbs(mchan) !time channel boundaries (musec)
	real*4 tcw(mchan),tcwv(mchan),tcws(mchan) !time channel widths (musec)
	real*4 ddeadtime,mdeadtime,adeadtime	! default detector, module and acquisition deadtimes
	real*4 ddeadtimem(mmod),mdeadtimem(mmod)	!detector, module, deadtimes for specific modules
	real*4 modulecnts(mchan,mmod,mcont+3) 	!sum of counts in each channel and module, for each sample type.
	real*4 moduleangle(mmod)	!module scattering angles 

	common/run_par_c80/info,headinfo,inst
	common/run_par_c3/ext,pathseparator
	common/run_par_c256/directinp,direct,deadtimefile
	common/run_par_i/nspec,nspecb,nproc,nchan,nchanv,nchans
     *,nchanb,nchanbv,nchanbs,nspch,nspchv,nspchs,nper,ndet
     *,nmod,specno,detno,udet,crateno,slotno,ndetpermodule,crate
     *,slot,module,lenpathseparator,lendirect,lendirectinp
	common/run_par_r/tcb,tcbv,tcbs,tcw,tcwv,tcws
     *,ddeadtimem,mdeadtimem,adeadtime,ddeadtime,mdeadtime,modulecnts
     *,moduleangle
