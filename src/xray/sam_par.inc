c
c 	sam_par.inc
c
c arrays to define the sample parameters
c
	integer*4 sgeom		!1 for cylindrical, = 2 for flat plate sample
	integer*4 nlambs	!no. of wavelengths for Sample corrections
	integer*4 ncont		!no. of cylinders plus sample
	integer*4 nselement(mcont) !no. of elements in each sample
	integer*4 smass(melement,mcont)	!mass numbers - 0 means natural element
	integer*4 stmon(mcont)	!1 = use transmission monitor, else 0
	integer*4 nreson	!no. of resonances for this sample
	real*4 sdimen1(mcont)	!Inner dimension of sample or container
	real*4 sdimen2(mcont)	!Outer dimension of sample or container
	real*4 srho(mcont)	!Number densities
	real*4 stemp(mcont)	!Sample temperature (for Placzek correction)
	real*4 sfrac(melement,mcont)	!atomic fraction of components
	real*4 satwt(melement,mcont)	!At wt.s of components
	real*4 scatlen(melement,mcont)	!bound scattering lengths of components
	real*4 sscatcs(melement,mcont)	!Bound scattering C/s of components
	real*4 scaptcs(melement,mcont)	!Capture C/S of components
	real*4 satwtav(mcont)		!average sample atomic weight
	real*4 scatlenav(mcont)		!average scattering length (sum c_i*b_i)
	real*4 scatlensqav(mcont)		!average square of the scattering length (sum c_i*b_i**2)
	real*4 sscatav(mcont)		!average bound scattering cross section 
	real*4 scaptav(mcont)		!average capture cross section 
	real*4 slamb(mcorrwav)	!sample wavelengths
	real*4 stscat(mcorrwav,mcont) !Tot. cross section for Sample + can(s)
	real*4 sheight		!overall height of sample
	real*4 tweak(mcont)	!Sample or container tweak factor
	real*4 snorm(mcont)	!sample normalisation factor
	real*4 srot(mcont)	!angle of rotation of sample
	real*4 swidth(mcont)	!full width of sample or can
	real*4 wavresf(mres),wavresl(mres) !min. and max. wavelengths of resonances
	real*4 scaptav_out(mcont)		!average capture cross section for outgoing radiation
	real*4 stscat_out(mcorrwav,mcont) !Tot. cross section for Sample + can(s) for outgoing radiation
	character*2 ssymbol(melement,mcont)	!Symbols of elements
	common/sam_par_i/sgeom,nlambs,ncont,nselement,smass,stmon,nreson
	common/sam_par_r/sdimen1,sdimen2,srho,stemp,sfrac,satwt,scatlen,sscatcs,scaptcs
     *,satwtav,scatlenav,scatlensqav,sscatav,scaptav
     *,slamb,stscat,sheight,tweak,snorm,srot,swidth,wavresf,wavresl
     *,scaptav_out,stscat_out
	common/sam_par_c/ssymbol
