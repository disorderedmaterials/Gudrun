      integer nwavmul            !number of wavelengths for m.s.
      integer nangmul            !number of angles for m.s.
      integer nazimul            !number of azimuthal angles for m.s.
      integer nanglesmul
      real wavmul(mcorrwav)            !wavelengths for m.s.
      real angmul(mcorrang)                  !angles for m.s.
      real azimul(mcorrang)                  !azimuthal angles for m.s.
      real rotmulang(mcorrang)       !Used to store rotation angle of flat plate samples
      real tthetamul(mcorrang)      !Actual scattering angles to be used
      real onescat(mcorrwav,mcorrang,mcorrang)      !single scattering d.s.c.s.
      real mulscat(mcorrwav,mcorrang,mcorrang)      !multiple scattering d.s.c.s.
      common/mul_corr_i/nwavmul,nangmul,nazimul,nanglesmul
      common/mul_corr_r/wavmul,angmul,azimul,rotmulang,tthetamul
     1,onescat,mulscat
     
