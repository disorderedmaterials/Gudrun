	subroutine fltabstof(ncan,abscort,rad1,rad2,rho,sigtl,abscs)
c
c to calculate flat plate absorption factors
c
c	original   by aks
c	modified   by wsh    10-3-88   for use with coral
c	modified by aks 16-7-90 to put out a unity absorption factor
c if secondary angle is close to 90 degrees.
c
c	completely revamped by aks 17-05-01 to become a subroutine of the gudrun
c suite
c
c	modified 5/02/2003 to include situation with azimuthal detectors
c
	implicit none
	include 'dimension.inc'
	include 'beam.inc'
	include 'abs_corr_azi.inc'
	integer ncan		!no. of containers + sample
	integer ic,ic1,ib,il
	real rad1(mcont)	!inner dimension of each slab
	real rad2(mcont)	!outer dimension of each slab
	real width		!lateral half width of sample and containers
	real rho(mcont)		!atomic number density of sample and containers
	real sigtl(mcorrwav,mcont)	!total cross-sections
	real abscs(mcont)		!capture cross section (not used)
	real mut(mcont),amut	!temporary atten. coeff. 
	real theta,theta1	!temporary angle values
	real fs,exp_dn,exp_up	!temporary exponential values
	real pi,piconv	
	real tsec,sec1,sec2	!secants of angles
	real ts(mcont)		!wall thickness of each container
	real ts_in,ts_out	!in and out flight path through each container
	real pl_in_dn(mcont)	!total flight path into downstream slab
	real pl_out_dn(mcont)	!total flight path out of downstream slab
	real pl_in_up(mcont)	!total flight path into uptream slab
	real pl_out_up(mcont)	!total flight path out of upstream slab
	real f		!function to calculate integral of exponential
	real abscort(mcorrwav,mcorrang,mcorrang,mcont)
	pi=3.141592653
	piconv=pi/180.
c
c step through banks and calculate corrections at each wavelength and
c scattering angle
c
	do ib=1,nangabs
		theta1=angabs(ib)
            rot=rotabsang(ib)*piconv
		sec1=1./cos(rot)
		tsec=theta1-rot1
c
c tsec is the angle the scattered beam makes with the normal to the sample
c surface.  if abs(tsec) is close to 90 deg. calculation of absorption
c coefficients is unreliable
c
		if(abs(abs(tsec)-90.).gt.1.) then
		tsec=tsec*piconv
		sec2=1./cos(tsec)
c
c first form the sum of path lengths (pl) through the other slabs into and out
c of each slab. the sample and containers are divided into two slabs, dn and up.
c up refers to the slab on the upstream side of the centre of the sample (i.e.
c facing towards the neutron source). dn refers to slabs on the downstream side 
c of the sample centre. in refers to the flight path into the slab. out refers
c to the flight out of the slab towards the detector.
c
		do il=1,nwavabs
		do ic=1,ncan
			pl_in_up(ic)=0.0
			pl_out_up(ic)=0.0
			pl_in_dn(ic)=0.0
			pl_out_dn(ic)=0.0
		end do
		do ic=1,ncan
			mut(ic)=rho(ic)*sigtl(il,ic)
			ts(ic)=(rad2(ic)-rad1(ic))
			ts_in=mut(ic)*ts(ic)*sec1
			ts_out=mut(ic)*ts(ic)*abs(sec2)
			do ic1=1,ncan
				if(sec2.gt.0.0) then
					if(ic1.lt.ic) then
						pl_in_up(ic1)=pl_in_up(ic1)+ts_in
						pl_out_up(ic1)=pl_out_up(ic1)+ts_out
						pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in
						pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out
					else
						pl_out_up(ic1)=pl_out_up(ic1)+ts_out
						pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in
					endif
				else if(sec2.lt.0.0) then
					if(ic1.gt.ic) then
						pl_in_dn(ic1)=pl_in_dn(ic1)+2.0*ts_in
						pl_out_dn(ic1)=pl_out_dn(ic1)+2.0*ts_out
					else if(ic1.eq.ic) then
						pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in
						pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out
					else
						pl_in_up(ic1)=pl_in_up(ic1)+ts_in
						pl_out_up(ic1)=pl_out_up(ic1)+ts_out
						pl_in_dn(ic1)=pl_in_dn(ic1)+ts_in
						pl_out_dn(ic1)=pl_out_dn(ic1)+ts_out
					endif
				endif
			end do
		end do
c
c calculate absorption integrals for each slab and form attenuation correction.
c
		do ic=1,ncan
			amut=mut(ic)
			fs=f(amut,ts(ic),sec1,sec2)
			exp_dn=exp(-(pl_in_dn(ic)+pl_out_dn(ic)))
			exp_up=exp(-(pl_in_up(ic)+pl_out_up(ic)))
			abscort(il,ib,1,ic)=fs*(ts(ic)*exp_up+ts(ic)*exp_dn)
     *	/(ts(ic)+ts(ic))
		end do
		end do
c
c case where tsec is close to 90 - no attenuation correction can be calculated
c
		else
			do ic=1,ncan
				do il=1,nwavabs
					abscort(il,ib,1,ic)=1.0
				end do
			end do
		endif
   	end do
	return
	end
