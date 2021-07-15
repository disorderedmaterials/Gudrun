      integer nwavabs            !number of wavelengths for a.c.
      integer nangabs            !number of angles for a.c.
      integer naziabs            !number of azimuthal angles for a.c.
      integer nanglesabs
      real wavabs(mcorrwav)            !wavelengths for a.c.
      real angabs(mcorrang)                  !angles for a.c.
      real aziabs(mcorrang)
      real rotabsang(mcorrang)       !Used to store rotation angle of flat plate samples
      real tthetaabs(mcorrang)      !Actual scattering angles to be used
      real abscor(mcorrwav,mcorrang,mcorrang,mcont1) !full corrections array
      common/abs_corr_i/nwavabs,nangabs,naziabs,nanglesabs
      common/abs_corr_r/wavabs,angabs,aziabs,tthetaabs,rotabsang,abscor
