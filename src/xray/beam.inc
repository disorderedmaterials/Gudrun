c
c       beam.inc
c
c
c arrays to define the beam parameters to be used for
c calculating corrections
c
      integer nprof            !no. of beam profile values
      integer nslice      !no. of slices for flat plate m.s. calculation
      integer ndeg            !no. of degrees between corrections
      real profil(mprof)      !Beam profile values
      real prstep            !step between profile values
      real stepa,stepm      !step size for attenuation and m.s. calculations
      real a,b            !beam width parameters (a>b)
      real hdown,hup      !bottom and top of beam from sample bottom
      real a1,b1            !scattered beam width parameters (a1 > b1)
      real hsdown,hsup      !bottom and top of scattered beam from sample b.
      real height            !overall height of sample
      real hbdown,hbup      !bottom and top of beam from sample centre
      real hsbdown,hsbup      !bottom and top of scattered beam from s. c.
      real bakfac,depression            !background factor and sample depression
      logical theta2theta,beamcompensation,detectorcompensation
      COMMON/ANGBLO_i/nslice,ndeg
      common/angblo_r/stepa,stepm
      COMMON/PROFBL_i/NPROF
      common/profbl_r/PROFIL,PRSTEP
      common/beam_r/a,b,hdown,hup,a1,b1,hsdown,hsup,height,hbdown,hbup
     *,hsbdown,hsbup,bakfac,depression
      common/beam_l/theta2theta,beamcompensation,detectorcompensation
