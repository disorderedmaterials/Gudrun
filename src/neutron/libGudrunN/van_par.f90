!     
! File:   van_par.f90
! Author: aks45
!
! Created on 01 November 2013, 13:00
!

MODULE van_par
    
    implicit none

!c
!c       van_par.inc
!c
!c arrays to define the vanadium parameters
!c
    integer vgeom            !1 for cylindrical, = 2 for flat plate V
    integer nlambv,nlambtransv      !no. of wavelengths for V corrections
    integer nvelement,mvelement      !no. of elements in vanadium
    integer, dimension(:), allocatable :: vmass!(melement)      !mass numbers - 0 means natural element
    integer vtmon            !1 = use transmission monitor, else 0
    integer nqv            !no. of q-values in V dcs file
    integer nbraggv      !no. of vanadium Bragg values
    integer ndettested,ndetrejected     !Used to store detector information about vanadium smooth.


    real vdimen1,vdimen2      !Inner and outer dimension of vanadium
    real vrho            !vanadium density
    real vtemp            !vanadium temperature for Placzek correction
    real, dimension(:), allocatable :: vfrac!(melement)      !atomic fraction of V components
    real, dimension(:), allocatable :: vatwt!(melement)      !At wt.s of V components
    real, dimension(:), allocatable :: vscatlen!(melement) !Scattering length of vanadium components
    real, dimension(:), allocatable :: vscatcs!(melement)      !Bound scattering C/s of V components
    real, dimension(:), allocatable :: vcaptcs!(melement)      !Capture C/S of V components
    real vatwtav            !average vanadium atomic weight
    real vscatlenav      !average scattering length of vanadium
    real vscatav            !average bound scattering cross section for V
    real vcaptav            !average capture cross section for V
    real, dimension(:), allocatable :: vlamb,vlambtrans!(mcorrwav)      !vanadium wavelengths
    real, dimension(:), allocatable :: vtscat,vsscat!(mcorrwav)      !Tot. and scattering cross-section for Van.
    real vheight            !Overall height of vanadium
    real vrot            !for flat plates angle of rotation
    real vwidth            !for flat plates full width of vanadium
    real, dimension(:), allocatable :: vqdcs!(mq)      !qvalues for dcs file
    real, dimension(:,:), allocatable :: vdcs!(mq,mgroup)      !differential c/s for vanadium (if present)
    real vsmear,vrmin      !values to 'unsmear' the vanadium dcs if appropriate
    real, dimension(:), allocatable :: braggq,braggh,braggw!(mq)
    real braggf !Used to generate V dcs
    real, dimension(:), allocatable :: vwavratio!(mchan)    !Wavelength values for vanadium smooth ratio
    real, dimension(:), allocatable :: vwavratioread!mchan)    !Wavelength values for vanadium smooth ratio read from previous run
    real, dimension(:), allocatable :: vsmoratio!(mchan)    !Summed ratio of raw counts to smoothed counts for vanadium
    real, dimension(:), allocatable :: vsmoratioread!(mchan) !Mean ratio of raw counts to smoothed counts for vanadium read from a previous run 
    real vanbackfraction    !Sets the lower limit of (van-background)/background

    character(len=2), dimension(:), allocatable :: vsymbol!(melement)      !list of chemical symbols
    character(len=256) vdcsname      !filename for calibration .DCS file

END MODULE
