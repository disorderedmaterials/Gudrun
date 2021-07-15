!     
! File:   multiple_scattering_modules.f90
! Author: aks45
!
! Created on 11 November 2013, 14:26
!

MODULE mul_corr_azi
    
    implicit none
    
    integer                                     :: nwavmul		!number of wavelengths for m.s.
    integer                                     :: nangmul		!number of angles for m.s.
    integer                                     :: nazimul		!number of azimuthal angles for m.s.
    integer                                     :: nan			!temporary array of number of annuli
    integer                                     :: ncval                !Number of cos(theta) values for arrays pc and czval
    real, dimension(:), allocatable             :: wavmul!(mcorrwav)		!wavelengths for m.s.
    real, dimension(:), allocatable             :: angmul!(mcorrang)			!angles for m.s.
    real, dimension(:), allocatable             :: azimul!(mcorrang)			!azimuthal angles for m.s.
    real, dimension(:,:,:), allocatable         :: onescat!(mcorrwav,mcorrang,mcorrang)	!single scattering d.s.c.s.
    real, dimension(:,:,:), allocatable         :: mulscat!(mcorrwav,mcorrang,mcorrang)	!multiple scattering d.s.c.s.
    real, dimension(:,:,:), allocatable         :: totscat!(mcorrwav,mcorrang,mcorrang)	!total scattering d.s.c.s.
    real, dimension(:,:,:), allocatable         :: pc!(mcval,mcorrwav,mcont) !Stores the current set of DCS datasets for each sample or container
    real, dimension(:), allocatable             :: czval!(mcval)  !Cos(theta) values for DCS for sample for m.s. calculation.
    real, dimension(:,:), allocatable           :: SIGSL,SIGTL!(mcorrwav,mcont)
    real, dimension(:), allocatable             :: rad1,rad2!(mcont)	!temporary dimension values
    real, dimension(:), allocatable             :: den!(mcont)		!density of each sample and cont.
    real, dimension(:), allocatable             :: captcs!(mcont)			!capture cross section of sample or cont
    real, dimension(:), allocatable             :: SIGS,SIGT,MUS!(mcont)
    real                                        :: theight    
    
END MODULE mul_corr_azi
    
MODULE van_mul_azi
       
    implicit none
    
    integer                                     :: nvwavmul		!number of wavelengths for vanadium m.s.
    integer                                     :: nvangmul		!number of angles for vanadium m.s.
    integer                                     :: nvazimul		!number of azimuthal angles for vanadium m.s.
    real, dimension(:), allocatable             :: vwavmul!(mcorrwave)	!wavelengths for vanadium m.s.
    real, dimension(:), allocatable             :: vangmul!(mcorrang)	!angles for vanadium m.s.
    real, dimension(:), allocatable             :: vazimul!(mcorrang)	!azimuthal angles for vanadium m.s.
    real, dimension(:,:,:), allocatable         :: vonescat!(mcorrwave,mcorrang,mcorrang) !single scattering d.s.c.s.
    real, dimension(:,:,:), allocatable         :: vmulscat!(mcorrwave,mcorrang,mcorrang) !multiple scattering d.s.c.s.

END MODULE van_mul_azi
    
MODULE sam_mul_azi
    
    implicit none
    
    integer                                     :: nswavmul		!number of wavelengths for vanadium m.s.
    integer                                     :: nsangmul		!number of angles for vanadium m.s.
    integer                                     :: nsazimul		!number of azimuthal angles for vanadium m.s.
    real, dimension(:), allocatable             :: swavmul!(mcorrwave)	!wavelengths for vanadium m.s.
    real, dimension(:), allocatable             :: sangmul!(mcorrang)	!angles for vanadium m.s.
    real, dimension(:), allocatable             :: sazimul!(mcorrang)	!azimuthal angles for vanadium m.s.
    real, dimension(:,:,:,:), allocatable         :: sonescat!(mcorrwave,mcorrang,mcorrang) !single scattering d.s.c.s.
    real, dimension(:,:,:,:), allocatable         :: smulscat!(mcorrwave,mcorrang,mcorrang) !multiple scattering d.s.c.s.

END MODULE sam_mul_azi
    

