c
c dimension sizes
c
      integer mspec            !dimensions of spectrum arrays
      integer mproc            !max. no. of spectra to process at one time
      integer mdet            !dimensions of detector arrays
      integer mchan            !dimensions of time channel arrays
      integer mmod            !dimensions of module arrays
      integer mcount            !dimensions of counts arrays
      integer mruns            !dimensions of file arrays
      integer mgroup            !dimensions of group and module arrays
      integer mperiod            !max. no. of DAE periods
      integer mspcgr            !max. no of spectra in a group
      integer mcorrwav      !max. no. of wavelengths in corrections.
      integer mcorrang      !max. no. of angles in corrections
      integer mq            !max. no. of Q values in final output
      integer mprof            !max. no. of beam profile values
      integer mcont            !max. no. of containers (incl. sample)
      integer mcont1            !no. of attenuation correction arrays needed
      integer mslice            !no. of slices for flat plate m.s. calc.
      integer melement      !max. no. of elements in a sample
      integer mcoeff            !max. no. of coefficients in polynomial
      integer mres            !no. of resonance regions to specify
      parameter (mspec=7500,mproc=600,mdet=10000,mchan=7000
     *,mmod=140,mcount=mchan*mproc
     *,mruns=100,mgroup=100,mperiod=20,mspcgr=2000,mcorrwav=1
     *,mcorrang=1001,mq=10000
     *,mprof=50,mcont=5,mcont1=(mcont*(mcont+1))/2,mslice=1000
     *,melement=20,mcoeff=20,mres=10)
