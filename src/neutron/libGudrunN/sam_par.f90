!     
! File:   sam_par.f90
! Author: aks45
!
! Created on 01 November 2013, 13:05
!

MODULE sam_par

!c
!c       sam_par.inc
!c
!c arrays to define the sample parameters
!c
      integer sgeom            !1 for cylindrical, = 2 for flat plate sample
      integer nlambs,nlambtranss      !no. of wavelengths for Sample corrections
      integer ncont            !no. of cylinders plus sample
      integer, dimension(:), allocatable :: nselement !(mcont) !no. of elements in each sample
      integer, dimension(:,:), allocatable :: smass!(melement,mcont)      !mass numbers - 0 means natural element
      integer, dimension(:), allocatable :: stmon!(mcont)      !1 = use transmission monitor, else 0
      integer nreson      !no. of resonances for this sample
      integer nxdcs,ngrpdcs
      integer, dimension(:), allocatable    :: ncvals!Number of cos(theta) values for using mint file to calculate multiple scattering
      integer mscont,mselement !Parameters for no of containers and rank of elements in arrays
      integer nexpon              !No. of exponential decays to use in subtracting background

      real, dimension(:), allocatable :: sdimen1!(mcont)      !Inner dimension of sample or container
      real, dimension(:), allocatable :: sdimen2!(mcont)      !Outer dimension of sample or container
      real, dimension(:), allocatable :: srho!(mcont)      !Number densities
      real, dimension(:), allocatable :: stemp!(mcont)      !Sample temperature (for Placzek correction)
      real, dimension(:,:), allocatable :: sfrac!(melement,mcont)      !atomic fraction of components
      real, dimension(:,:), allocatable :: satwt!(melement,mcont)      !At wt.s of components
      real, dimension(:,:), allocatable :: sscatlen!(melement,mcont)      !bound scattering lengths of components
      real, dimension(:,:), allocatable :: sscatcs!(melement,mcont)      !Bound scattering C/s of components
      real, dimension(:,:), allocatable :: scaptcs!(melement,mcont)      !Capture C/S of components
      real, dimension(:), allocatable :: satwtav!(mcont)            !average sample atomic weight
      real, dimension(:), allocatable :: sscatlenav!(mcont)            !average scattering length (sum c_i*b_i)
      real, dimension(:), allocatable :: sscatlensqav!(mcont)            !average square of the scattering length (sum c_i*b_i**2)
      real, dimension(:), allocatable :: sscatav!(mcont)            !average bound scattering cross section 
      real, dimension(:), allocatable :: scaptav!(mcont)            !average capture cross section 
      real, dimension(:), allocatable :: slamb,slambtrans!(mcorrwav)      !sample wavelengths
      real, dimension(:,:), allocatable :: stscat,ssscat!(mcorrwav,mcont) !Tot. and scattering cross section for Sample + can(s)
      real sheight            !overall height of sample
      real, dimension(:), allocatable :: tweak!(mcont)      !Sample or container tweak factor
      real, dimension(:), allocatable :: snorm!(mcont)      !sample normalisation factor
      real, dimension(:), allocatable :: srot!(mcont)      !angle of rotation of sample
      real, dimension(:), allocatable :: swidth!(mcont)      !full width of sample or can
      real, dimension(:), allocatable :: wavresf,wavresl!(mres) !min. and max. wavelengths of resonances
      real, dimension(:), allocatable :: sxdcs!(mq)      !qvalues for dcs file
      real, dimension(:,:), allocatable :: sdcs,esdcs!(mgroup,mq)      !differential c/s for vanadium (if present)
      real, dimension(:), allocatable :: czvals!(mcval)  !Cos(theta) values for DCS for sample for m.s. calculation.
      real, dimension(:,:,:), allocatable :: pcs!(mcval,mcorrwav,mcont) !Stores the DCS for this sample as derived from a mint file if input.
      real, dimension(:), allocatable :: sampleenvscatfrac,sampleenvamut!(mcont)    !Used to calculate effect of altered geometry between sample and containers
      real, dimension(:), allocatable :: exponamp,expondecay      !Stretched exponential parameters


      character(len=2), dimension(:,:), allocatable :: ssymbol!(melement,mcont)      !Symbols of elements
      character(len=256), dimension(:), allocatable :: smintname!(mcont) !Merged interference DCS filename to be used to estimate mulitple scattering. Default is *

END MODULE
