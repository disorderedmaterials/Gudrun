!     
! File:   beam.f90
! Author: aks45
!
! Created on 27 April 2012, 12:52
!

MODULE beam_routines
      
    implicit none

!c arrays to define the beam parameters to be used for
!c calculating corrections
!c
    integer                         :: nprof        !no. of beam profile values
    integer                         :: nslice    !no. of slices for flat plate m.s. calculation
    integer                         :: ndeg        !no. of degrees between corrections
    real, dimension(:), allocatable :: profil    !Beam profile values
    real                            :: prstep        !step between profile values
    real                            :: stepa,stepm    !step size for attenuation and m.s. calculations
    real                            :: a,b        !beam width parameters (a>b)
    real                            :: hdown,hup    !bottom and top of beam from sample bottom
    real                            :: a1,b1        !scattered beam width parameters (a1 > b1)
    real                            :: hsdown,hsup    !bottom and top of scattered beam from sample b.
    real                            :: hbdown,hbup    !bottom and top of beam from sample centre
    real                            :: hsbdown,hsbup    !bottom and top of scattered beam from s. c.
    real                            :: bakfac        !background factor

CONTAINS
!***********************************************************************************
!*
!*	init_beam.FOR
!*
!*	A K Soper, March 2001
!*	
!*	Reads the beam profile and other parameters from unit nin 
!*       to be used for the calculation of attenuation and multiple 
!*	scattering corrections.
!*
!***********************************************************************************
    subroutine init_beam(nin)

        use inputfilestrings
        use reallocation_routines
        use run_par
        use spectrum_par

        implicit none
      	integer nin		!unit number to read data from 
        integer i		!do loop counter
        real am		!maximum value of beam profile
        real lowang,highang	!lowest and highest scattering angles for corrs.
        real range		!range of scattering angles to be used.
        real pi		!pi
!C
!C READ NO. OF PROFILE VALUES AND VALUES
!C
        READ(NIN,*) NPROF
        call reallocate1d_r(profil,nprof)
        READ(NIN,*) (PROFIL(I),I=1,NPROF)
!C
!C NORMALIZE PROFILE VALUES TO MAXIMUM VALUE
!C
        AM=0.
        DO I=1,NPROF
            AM=MAX(AM,PROFIL(I))
        end do
        DO I=1,NPROF
            PROFIL(I)=PROFIL(I)/AM
        end do
!c
!C INPUT INTEGRATION PARAMETERS
!C
        READ(NIN,*) stepa,stepm,nslice
!C
!C INPUT NO. OF degrees between corrections
!C
        READ(NIN,*) NDEG
!C
!C READ position of edges of incident beam, and scattered beam relative to 
!c centre of sample.
!C
        READ(NIN,*) B,A,HBDOWN,HBUP
        read(nin,*) B1,A1,HSBDOWN,HSBUP
!c
!c define profile step
!c
        PRSTEP=(A-B)/(NPROF-1)
!c
!c Name of file containing the incident spectrum parameters
!c
        call getaline(nin)
        specname=line(ncf(1):ncl(nwords))
!c           write(6,*) specname
!c
!c Read the parameters
!c
        call get_spectrum_par
        write(6,*) ef0,phiepi,et,phimax,w1,w2,alpha,sqrtk,sqrtl
        write(6,*) specname
!c
!c Read background factor
!c
        read(nin,*) bakfac
!c
!c Set the individual background factors for those groups which have not already been set
!c
!            do i=1,mgroup
!                if(setbackground_factor(i).eq.0) background_factor(i)=bakfac
!            end do
        return
    end subroutine init_beam

    FUNCTION PROBE(X,PROFIL,NPROF,PRSTEP,A1,B1)                       !MUL06180

        integer nprof
        real PROFIL(*)
        real probe,x,prstep,a1,b1,diff,apos
        integer npos,npos1,npos2
        
        IF(X.GT.A1.OR.X.LT.B1) then
            PROBE=0.
            RETURN
        else
            DIFF=A1-X
            APOS=DIFF/PRSTEP
            NPOS=INT(APOS)
            DIFF=APOS-FLOAT(NPOS)
            NPOS1=NPOS+1
            NPOS2=NPOS+2
            IF(DIFF.EQ.0.) then
                PROBE=PROFIL(NPOS1)
                RETURN
            else
                PROBE=PROFIL(NPOS1)+(PROFIL(NPOS2)-PROFIL(NPOS1))*DIFF
                RETURN
            end if
        end if
    END function probe

END MODULE beam_routines
