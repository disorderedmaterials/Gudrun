c
c arrays to define the spectra grouping
c
	character*80 fname_grps		!name of groups file
	integer*4 ngroup		!number of groups of spectra 
	integer*4 nbread		!number of spectra read in spec.bad
	integer*4 ndetgr(mgroup)	!no. of spectra in each group
	integer*4 indgr(mspcgr,mgroup)	!spectrum number for each group
	integer*4 igrp(mspec)		!group number for each spectrum
	integer*4 ibad(mspec)		!bad spectrum index
	integer*4 spike(mspec)		!number of spikes in a spectrum
	integer*4 nchfir,nchlas		!range of channels for spike search
	integer*4 flagsp		!non-zero to treat spikes and steps
	real*4 lengrp(mgroup)		!average secondary flight path for group
	real*4 tthgrp(mgroup)		!average scattering angle for group
	real*4 phigrp(mgroup)		!average azimuthal angle for group
	real*4 rmsfac,stdfac	!number of rms deviations to allow and range of standard devs.
	real*4 qmingrp(mgroup),qmaxgrp(mgroup)	!Min and max Q value for each group
	real*4 background_factor(mgroup)
	real*4 spikedev			!square of deviation for spike checking
	common/groups_c/fname_grps
	common/groups_i/ngroup,nbread,ndetgr,indgr,igrp,ibad,spike
     *,nchfir,nchlas,flagsp
	common/groups_r/lengrp,tthgrp,phigrp,rmsfac,stdfac,qmingrp
     $,qmaxgrp,background_factor,spikedev
