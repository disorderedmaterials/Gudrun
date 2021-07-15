c
c global calibration variables common to all routines
c
	integer*4 phiutno		!User table column containing phi values
	real*4 lenin			!incident flight path
	real*4 incidentfp		!Incident flight path from input file
	real*4 incidentwl		!Incident wavelength (only used for reactor instruments)
	real*4 deltat(mdet)		!time offset in musec for each detector
	real*4 lendet(mdet)		!scattered flight path
	real*4 ttheta(mdet)		!detector scattering angle
	real*4 phi(mdet)		!detector azimuthal angles
	common/calibration_i/phiutno
	common/calibration_r/lenin,incidentfp,incidentwl,deltat,lendet
     1,ttheta,phi
