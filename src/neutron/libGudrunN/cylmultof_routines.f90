!     
! File:   cylmultop_routines.f90
! Author: aks45
!
! Created on 11 November 2013, 12:16
!

MODULE cylmultof_routines
    
    USE reallocation_routines
    use corrections_routines
    
    implicit none
    
!C                                                                     
!C ARRAYS - PRIMARY SCATTERING:-                                       
!C        SUMC = C-C                                                     
!C        SUMS = C-P                                                    
!C        SUMT = P-C                                                     
!C        SUMP = P-P                                                    
!C                                                                      
    real, dimension(:,:,:), allocatable                 :: SUMC!(mcorrang,mcorrang,mcont)
    real, dimension(:,:,:), allocatable                 :: SUMS!(mcorrang,mcorrang,mcont)
    real, dimension(:,:,:), allocatable                 :: SUMT!(mcorrang,mcorrang,mcont)
    real, dimension(:,:,:), allocatable                 :: SUMP!(mcorrang,mcorrang,mcont)
!C  ARRAYS - SECONDARY SCATTERING:-                                      
!C                                                                      
!C        SUMA = C-C                                                     
!C        SUMD = C-P                                                    
!C        SUME = P-C                                                    
!C        SUMB = P-P                                                    
    real, dimension(:,:,:,:), allocatable               :: SUMA!(mcorrang,mcorrang,mcont,mcont)
    real, dimension(:,:,:,:), allocatable               :: SUMD!(mcorrang,mcorrang,mcont,mcont)
    real, dimension(:,:,:,:), allocatable               :: SUME!(mcorrang,mcorrang,mcont,mcont)
    real, dimension(:,:,:,:), allocatable               :: SUMB!(mcorrang,mcorrang,mcont,mcont)

    integer                                             :: NS,nkrad,nomegamax
    integer                                             :: NBLEN1!=2000                                                      
    integer                                             :: NBLEN2!=10000
    integer                                             :: ms,nz,minb,maxb,mins,maxs
    integer, dimension(:), allocatable                  :: NRAD!(20)
    integer, dimension(:), allocatable                  :: NOMEGA!(20)
    integer, dimension(:), allocatable                  :: NT2!(10000)

    real                                                :: HIGHT,HIGHB,HIGHS,HIGHBS,ZSTEP
    real, dimension(:), allocatable                     :: OM!(2000)
    real, dimension(:), allocatable                     :: ALEN1!(2000)
    real, dimension(:,:), allocatable                   :: ALEN2!(mcorrang,10000)
    real, dimension(:), allocatable                     :: RST!(20)
    real, dimension(:), allocatable                     :: R1!(20)
    real, dimension(:), allocatable                     :: RSQ!(20)
    real, dimension(:), allocatable                     :: AREAEL!(20)
    real, dimension(:), allocatable                     :: NTEST!(2000)

    CONTAINS
    
!C                                                                       MUL02620
!C        **********************                                         MUL02630
!C        *                    *                                         MUL02640
!C        *   SUBROUTINE LEN3  *                                         MUL02650
!C        *                    *                                         MUL02660
!C        **********************                                         MUL02670
!C                                                                       MUL02680
    SUBROUTINE LEN3(nw,ASTEP)
        
        use mul_corr_azi

        REAL                :: astep,pi,X,Y,Z,rstep,radd,r,omegst,fac
        integer             :: nw,nom,krad,nomsta,i,j,k,ns,nomeg,nfact
        !C                                                                       MUL02740
!C THIS SUB. SETS UP ARRAY TO BE USED IN LEN4,LEN5 AND SUMMU1            MUL02750
!C NAN IS THE THE NO. OF ANNULI.                                         
!C MS IS THE NUMBER OF STEPS FROM rad1 TO RAD2 - IT IS USED TO           
!C DETERMINE THE SIZE OF THE STEP IN EACH ANNULUS                        MUL02780
!C THE ITH ANNULUS IS DEFINED TO LIE BETWEEEN RAD1(I) AND RAD2(I)        MUL02790
!C WITHE A DENSITY OF DEN(MOLS PER CC), SCATTERING C/S SIGS(I)           MUL02800
!C BARNS AND TOTAL ATTENUATION C/S SIGT(I) BARNS,AND TOTAL               MUL02810
!C ABSORPTION COEFFICIENT MUS(I) PER CM.  ON EXIT, NRAD(I) IS THE        MUL02820
!C NUMBER OF ELEMENTAL RINGS TO BE INTEGRATED OVER,                      MUL02830
!C NOMEGA(I) IS THE NUMBER OMEGA VALUES IN THE ITH RING, AND OM(K)       MUL02840
!C IS THE KTH OMEGA VALUE. RSQ(I) IS R1(I)**2 WHERE R1(I)                MUL02850
!C IS THE RADIUS OF THE ITH RING                                         MUL02860
!C                                                                       MUL02870
!C                                                                       MUL02880
!C SET A NUMBER TO CHECK THE NUMBER OF ELEMENTS IN R1 DOES NOT EXCEED 20 MUL02890
!C                                                                       MUL02900
!        write(6,*) 'len3> HIGHB,HIGHS,HIGHBS,HIGHT: ',HIGHB,HIGHS,HIGHBS,HIGHT
!	write(6,*) 'len3> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
        PI=4.0*atan(1.0)
        call reallocate1d_i(nrad,nan)
        call reallocate1d_r(rst,nan)
        NKRAD=0
        nblen1=0
        NOM=0
        KRAD=0
        NOMSTA=1
        nomegamax=0
        DO I=1,NAN
!C                                                                       MUL02980
!C DO NOT WISH TO INTEGRATE OVER NON-SCATTERING ANNULI                   MUL02990
!C                                                                       MUL03000
            IF(den(I).gt.0.0) then
!C                                                                       MUL03020
!C DEFINE SCATTERING C/S PER UNIT VOLUME OF ITH ANNULUS                  MUL03030
!C                                                                       MUL03040
                FAC=DEN(I)*SIGS(I)/(4*PI)
!C                                                                       MUL03060
!C DEFINE NO. OF RADIUS STEPS IN THE ITH ANNULUS                         MUL03070
!C                                                                       MUL03080
                Y=RAD1(I)                                                         
                Z=RAD2(i)-Y                                                       
                X=Z/ASTEP+0.5
                NS=INT(X)
                IF(NS.LT.1) NS=1
                NRAD(I)=NS
                RSTEP=Z/NS
                RST(I)=RSTEP
                RADD=-0.5*RSTEP+RAD1(I)                                           
!C                                                                       MUL03190
!C SET UP NS RADIUS VALUES                                               MUL03200
!C                                                                       MUL03210
                DO J=1,NS
                    R=J*RSTEP+RADD
                    KRAD=KRAD+1
!C                                                                       MUL03250
!C CHECK THAT KRAD IS NOT GREATER THAN NKRAD                             MUL03260
!C                                                                       MUL03270
                    IF(KRAD.GT.NKRAD) then
                        nkrad=nkrad+20
                        call reallocate1d_r(r1,nkrad)
                        call reallocate1d_r(rsq,nkrad)
                        call reallocate1d_r(areael,nkrad)
                        call reallocate1d_i(nomega,nkrad)
                    end if
                    R1(KRAD)=R
                    RSQ(KRAD)=R*R
!C                                                                       MUL03310
!C DEFINE NO. OF STEPS IN OMEGA FOR EACH RADIUS VALUE                    MUL03320
!C                                                                       MUL03330
!C THE NUMBER OF ELEMENTS IN ANY GIVEN RING IS SET TO BE AN INTEGER       MUL03340
!C MULTIPLE OF THE NUMBER IN THE PREVIOUS RING                           MUL03350
!C                                                                       MUL03360
                    X=PI*R/ASTEP
                    nfact=INT(X/real(NOMSTA)+0.5)
                    NOMEG=nfact*NOMSTA
                    NOMSTA=NOMEG
                    OMEGST=PI/real(NOMEG)
!C                                                                       MUL03460
!C DEFINE THE AREA ELEMENTS AT THIS RADIUS AND MULTIPLY BY               MUL03470
!C THE APPROPRIATE SCATTERING C/S                                        MUL03480
!C                                                                       MUL03490
                    AREAEL(KRAD)=R*OMEGST*RSTEP*FAC
                    NOMEG=2*NOMEG
                    NOMEGA(KRAD)=NOMEG
                    nomegamax=max(nomegamax,nomeg)
                    DO K=1,NOMEG
                        NOM=NOM+1
!C                                                                       MUL03550
!C CHECK THAT NOM IS NOT GREATER THAN NBLEN1                             MUL03560
!C                                                                       MUL03570
                        IF(NOM.GT.NBLEN1) then
                            nblen1=nblen1+200
                            call reallocate1d_r(om,nblen1)
                        end if
                        OM(NOM)=(K-1)*OMEGST
                    end do
                end do
            else
                NRAD(I)=0
            end if
        end do
!        write(6,*) 'len3> array sizes', size(r1),size(rsq),size(areael),size(nomega)
        RETURN
      
    END subroutine len3
    
!C                                                                       MUL03730
!C        *********************                                          MUL03740
!C        *                   *                                          MUL03750
!C        *  SUBROUTINE LEN4  *                                          MUL03760
!C        *                   *                                          MUL03770
!C        *********************                                          MUL03780
!C                                                                       MUL03790
    SUBROUTINE LEN4(nw)

        use beam_routines
        use mul_corr_azi
        
        real r,om1,d,sum,xstart,xnext
        integer i,j,k,l,nom,nra,nr,nw,nomeg
!C                                                                       MUL03860
!C THIS SETS UP THE ARRAY ALEN CORRESPONDING TO THE EXPONENTIAL          MUL03870
!C SUM EXP(-SUMMATION(MUS*L1(R,OMEGA))).  THE SUMMATION IS               MUL03880
!C OVER THE ANNULI AND L1 IS THE DISTANCE OF THE                         MUL03890
!C CURRENT ANNULUS THROUGH THAT ANNULUS IN A DIRECTION PARALLEL          MUL03900
!C TO INCIDENT BEAM.                                                     MUL03910
!C ALSO AN ARRAY NTEST IS SETUP WHICH IS 1 IF A GIVEN ELEMNT IS          MUL03920
!C IN INCIDENT BEAM AND ZERO IF NOT                                      MUL03930
!C                                                                       MUL03940
!C CALL DIST3                                                            MUL03950
!C                                                                       MUL03960
!        write(6,*) 'len4> HIGHB,HIGHS,HIGHBS,HIGHT: ',HIGHB,HIGHS,HIGHBS,HIGHT
!	write(6,*) 'len4> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
        NOM=0
        NRA=0
        call reallocate1d_r(ntest,nblen1)
        call reallocate1d_r(alen1,nblen1)
!C                                                                       MUL03990
!C STEP THROUGH ANNULUS VALUES                                           MUL04000
!C                                                                       MUL04010
        DO I=1,NAN
            NR=NRAD(I)
            IF(NR.gt.0) then
                DO J=1,NR
                    NRA=NRA+1
                    R=R1(NRA)
                    NOMEG=NOMEGA(NRA)
!
!C STEP THROUGH OMEGA VALUES
!C                                                                       MUL04110
                    DO K=1,NOMEG
                        NOM=NOM+1
                        OM1=OM(NOM)
!                                                                         MUL04150
!C TEST WHETHER CURRENT ELEMENT IS IN BEAM                               MUL04160
!C                                                                       MUL04170
                        D=R*SIN(OM1)
                        NTEST(NOM)=PROBE(D,PROFIL,NPROF,PRSTEP,A,B)
!                                                                         MUL04200
!C AT EACH R1,OM1 VALUE L1 IS CALCULATED THROUGH EACH ANNULUS ALONG      MUL04210
!C A LINE PARALLEL TO OM0=0, IN THE DIREDTION LEFT TO RIGHT              MUL04220
!C                                                                       MUL04230
                        SUM=0
                        DO L=1,NAN
!C                                                                       MUL04280
!C CALCULATE THE DISTANCE THROUGH A CIRCLE OF RADIUS RAD(L1)             MUL04290
!C                                                                       MUL04300
                            XSTART=DIST3(R,RAD1(L),OM1)                                       
                            XNEXT=DIST3(R,RAD2(L),OM1)
!C                                                                       MUL04320
!C XSTART IS THE VALUE THAT DISTANCE THROUGH THE NEXT SMALLEST           MUL04330
!C ANNULUS- HENCE IF WE SUBTRACT XSTART FORM XNEXT, THIS WILL LEAVE      MUL04340
!C THE DISTANCE THROUGH THE CURRENT ANNULUS                              MUL04350
!C                                                                       MUL04360
                            SUM=SUM+MUS(L)*(XNEXT-XSTART)
                        end do
!C                                                                       MUL04400
!C BEFORE OUTPUT WE MUST RECONSTRUCT SCATTERING FUNCTIONS FOR EACH       MUL04410
!C ANNULUS USING THE ACTUAL AREAS FOR THE GIVEN GEOMETRY                 MUL04420
!C                                                                       MUL04430
                        ALEN1(NOM)=EXP(-SUM)
                        end do
                end do
            end if
        end do
        RETURN 
    END subroutine len4    
      
!C                                                                       MUL04500
!C         *********************                                         MUL04510
!C         *                   *                                         MUL04520
!C         *  SUBROUTINE LEN5  *                                         MUL04530
!C         *                   *                                         MUL04540
!C         *********************                                         MUL04550
!C                                                                       MUL04560
    SUBROUTINE LEN5(NW)

        use beam_routines
        use mul_corr_azi

        integer nw,nalen,nom,nra,nsave,nomeg,i,j,k,l,iang,iazi
        real pi,piconv,r,angle,om1,d,sum,xstart,xnext,detfac
        !C                                                                       MUL04630
!C THIS SETS UP THE ARRAY ALEN2 CORRESPONDING TO THE VALUES              MUL04640
!C OF EXP(-SUMMATION(MUS*L2(R,OMEGA))), WHERE L2(R,OMEGA) IS THE         MUL04650
!C DISTANCE THROUGH A GIVEN ANNULUS FROM THE SCATTERING POINT            MUL04660
!!C (R,OMEGA) PARALLEL TO THE THE DIRECTION OF THE                        MUL04670
!C COUNTER.  THE SUMMATION IS OVER ALL THE ANNULI, AND THERE IS          MUL04680
!C AN ARRAY ELEMENT FOR EACH SCATTERING POINT WITHIN THE SAMPLE          MUL04690
!C                                                                       MUL04700
!C CALLS DIST2                                                           MUL04710
!        write(6,*) 'len5> HIGHB,HIGHS,HIGHBS,HIGHT: ',HIGHB,HIGHS,HIGHBS,HIGHT
!	write(6,*) 'len5> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
        PI=4.0*atan(1.0)
        piconv=pi/180.0
        nblen2=0
        NALEN=0
        NOM=0
        NRA=0
!C                                                                       MUL04760
!C STEP THROUGH THE ANNULI AND RADIUS VALUES                             MUL04770
!C                                                                       MUL04780
        DO I=1,NAN
            IF(NRAD(I).gt.0) then
                DO J=1,NRAD(I)
                    NRA=NRA+1
                    R=R1(NRA)
                    NOMEG=NOMEGA(NRA)
!C                                                                       MUL04860
!C STEP THROUGH THE OMEGA VALUES                                         MUL04870
!C                                                                       MUL04880
                    NSAVE=NOM
                    DO iazi=1,nazimul
                        NOM=NSAVE
                        ANGLE=azimul(iazi)*piconv
                        DO K=1,NOMEG
                            NOM=NOM+1
                            OM1=OM(NOM)-ANGLE
!C                                                                       MUL04960
!C DECIDE WHETHER SCATTERING POINT IS WITHIN THE SCATTERED BEAM -        MUL04970
!C IF SO WE SET NTEST TO UNITY, OTHERWISE IT IS ZERO.                    MUL04980
!C                                                                       MUL04990
                            NALEN=NALEN+1
!C                                                                       MUL05010
!C CHECK NALEN IS NOT GREATER THAN NBLEN2                                MUL05020
!C                                                                       MUL05030
                            IF(NALEN.GT.NBLEN2) then
                                nblen2=nblen2+200
                                call reallocate1d_i(nt2,nblen2)
                                call reallocate2d_r(alen2,nangmul,nblen2)
                            end if
                            NT2(NALEN)=1
                            D=R*SIN(OM1) 
                            IF(D.GT.A1.OR.D.LT.B1) NT2(NALEN)=0
!C                                                                       MUL05080
!C AT EACH R,OM1 VALUE THE DISTANCE L2 IS CALCULATED THROUGH EACH        MUL05090
!C ANNULUS, ALONG A LINE PARALLEL TO OMEGA=0, IN THE DIRECTION           MUL05100
!C LEFT TO RIGHT                                                         MUL05110
!C                                                                       MUL05120
                            SUM=0.
                            DO L=1,NAN
!C                                                                       MUL05170
!C CALCULATE THE DISTANCE THROUGH A CIRCLE OF RADIUS RAD(L)              MUL05180
!C                                                                       MUL05190
                                Xstart=DIST2(R,RAD1(L),OM1)
                                XNEXT=DIST2(R,RAD2(L),OM1)
!C                                                                       MUL05210
!C XSTART IS THE VALUE OF THAT DISTANCE THROUGH THE NEXT                 MUL05220
!C SMALLEST ANNULUS. HENCE IF WE SUBTRACT XSTART FROM XNEXT,             MUL05230
!C THIS WILL LEAVE THE DISTANCE THROUGH THE CURRENT ANNULUS              MUL05240
!C                                                                       MUL05250
                                SUM=SUM+MUS(L)*(XNEXT-XSTART)
                            end do
!c
!c step through the detector azimuthal angles (corresponds to angmul in 
!c cylindrical sample coordinates)
!c
                            do iang=1,nangmul
                                detfac=abs(sin(angmul(iang)*piconv))
!c
!c perform azimuthal correction
!c
                                if(detfac.gt.1.0e-10) then
                                    ALEN2(iang,NALEN)=EXP(-SUM/detfac)
                                else
                                    alen2(iang,nalen)=0.0
                                endif
                            end do
                        end do 
                    end do
                end do
            end if
        end do
        RETURN
    END subroutine len5
      
      
!C                                                                       MUL07450
!C       ***********************                                         MUL07460
!C       *                     *                                         MUL07470
!C       *  SUBROUTINE SUMMU1  *                                         MUL07480
!C       *                     *                                         MUL07490
!C       ***********************                                         MUL07500
!C                                                                       MUL07510
    SUBROUTINE SUMMU1()

        use beam_routines
        use mul_corr_azi
                
        !Internal variables
        
        real, dimension(:), allocatable                     :: AREAP,AREAS,AREAT,AREAC!mcorrang
 	real, dimension(:), allocatable                     :: tema1,tema2,temb1,temb2,temd1,temd2,teme1,teme2!(mcorrang)
        real, dimension(:), allocatable                     :: STORSQ,STORSM
        real, dimension(:), allocatable                     :: AFACT,DFACT,EFACT,BFACT,ZMULL,ZVALSQ
        real, dimension(:), allocatable                     :: ZERO,EQRAD
        real, dimension(:,:,:), allocatable                 :: SUMHA,SUMHB,SUMHD,SUMHE!(mcorrang,mcorrang,2)
        real, dimension(:,:), allocatable                   :: AREAA,AREAB,AREAD,AREAE!(mcorrang,2)
        real, dimension(:,:), allocatable                   :: PRODC,PRODS,PRODT,PRODP!(mcorrang,mcorrang)
        
        real :: pi,alimit,zstep,zadd,zsteps,highcc,highcp,highpc,highpp,z,arad,r1t,r2t,d,e &
        ,area,area1,area2,arprod,amu,add,adda,prod1a,prod1d,prod1e,prod1b,prod2a,prod2e,prod2d,prod2b &
        ,omcos,omdif,om2,om1,r12,r1sq,r2sq,radd,radnx1,radnx2,radst1,radst2,rdiff,rstep,rt1t,rt1t1,rt1t2 &
        ,sigcos,store,storm,suma1,suma2,sume1,sume2,sumd1,sumd2,sumb1,sumb2,suml,sumr1,sumr2,tstep,xadd,xnext &
        ,xstart,xprod,zmul,zmula,zmule,zmuld,zmulb,zsq,zsqrt,sumar
        integer     :: i,iz1,iz2,i1,i2,j,j1,j2,j2star,k,k1,k2,iazi,iang,nom1,nom2,nom1s,nom2s,nr1,nr1s,nserch,izp,nomeg1,nomeg2 &
        ,nlast,nlast1,ndom1,nalen1,nalen2,na,ian,ntemp1,ntemp2,nt2t,nt2t1,nt2t2,nserc1,nserc2,nsave,nrat,nrad1,nrad2 &
        ,nr2s,nr2
     
!        write(6,*) 'summu1> HIGHB,HIGHS,HIGHBS,HIGHT: ',HIGHB,HIGHS,HIGHBS,HIGHT
!        write(6,*) 'summu1> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
        PI=4.0*atan(1.0)
!C                                                                       MUL07690
!C THIS ROUTINE CALCULATES THE PRIMARY AND SECONDARY SCATTERING          MUL07700
!C FOR A SERIES OF CONCENTRIC INCOHERENT SCATTERERS OF ANNULAR           MUL07710
!C CROSS-SECTION WITH A VARIETY OF THICKNESSES, USING EQ (13) OF         MUL07720
!C BLECH AND AVERBACH 1965 AS A STARTING POINT.  THE ELEMENTAL           MUL07730
!C AREAS ARE DEFINED IN SUB LEN3.  THE GEOMETRY HAS BEEN DEFINED         MUL07740
!C SO AS TO KEEP THE NUMBER OF COMPUTATIONS OF  INTER-                   MUL07750
!C ELEMENTAL DISTANCES (AND ANGLES) TO A MINIMUM --SEE NOTES FOR         MUL07760
!C A FULLER EXPLANATION.                                                 MUL07770
!C   THE INTEGRAL IS CALCULATED FOR A SERIES OF Z' VALUES AND AT         MUL07780
!C THE END THE SUMS ARE WEIGHTED BY THE NUMBER OF TIMES THEY APPEAR      MUL07790
!C IN THE FINAL INTEGRAL AND ARE ADDED TOGETHER                          MUL07800
!C  MOST OF THE VARIABLE NAMES HAVE EITHER BEEN INTRODUCED IN            MUL07810
!C THE CALL ROUTINE OR THEIR USE SHOULD BE DEDUCEABLE FROM THE           MUL07820
!C PROGRAM COMMENTS                                                      MUL07830
!C  THE BEAM LIES BETWEEN X=A AND X=B, WHERE X=0 IS AN AXIS PARALLEL     MUL07840
!C TO THE DIRECTION OF THE BEAM THROUGH CENTRE OF THE SYSTEM.  THE       MUL07850
!C INTEGRAL IS CALCULATED BOTH FOR WHEN THE BEAM IS PARTIALLY            MUL07860
!C EMMERSED AND FOR WHEN COMPLETELY EMMERSED                             MUL07870
!C                                                                       MUL07880
!C   THE VALUES OF ALEN1 AND ALEN2 ARE READ FROM BINARY FILES            MUL07890
!C                                                                       MUL07900
!C INITIALIZE VARIOUS ARRAYS                                             MUL07910
!C                                                                       MUL07920
        ALIMIT=1.0E-30
        ZSTEP=HIGHT/NZ
        ZADD=-ZSTEP
        ZSTEPS=ZSTEP*ZSTEP
!C                                                                       MUL07970
!C CALCULATE HIGHT OF BEAM FOR PARTIAL EMMERSION                        MUL07980
!C                                                                       MUL07990
        HIGHB=(MAXB-MINB+1)*ZSTEP
        HIGHS=(MAXS-MINS+1)*ZSTEP
        HIGHCC=HIGHT*HIGHT
        HIGHCP=HIGHT*HIGHS
        HIGHPC=HIGHB*HIGHT
        HIGHPP=HIGHB*HIGHS
!        write(6,*) 'summu1> HIGHB,HIGHS,HIGHBS,HIGHT: ',HIGHB,HIGHS,HIGHBS,HIGHT
!        write(6,*) 'summu1> minb,maxb,mins,maxs ',minb,maxb,mins,maxs
!c
!c zero accumulators
!c
!        write(6,*) 'summu1> nangmul,nazimul,nan,nz ',nangmul,nazimul,nan,nz
	call reallocate4d_r(suma,nangmul,nazimul,nan,nan)
	call reallocate4d_r(sumb,nangmul,nazimul,nan,nan)
	call reallocate4d_r(sumd,nangmul,nazimul,nan,nan)
	call reallocate4d_r(sume,nangmul,nazimul,nan,nan)
        call reallocate1d_r(zero,size(r1))
        call reallocate1d_r(eqrad,size(r1))
        call reallocate1d_r(zvalsq,nz)
        call reallocate1d_r(zmull,nz)
        call reallocate1d_r(afact,nz)
        call reallocate1d_r(dfact,nz)
        call reallocate1d_r(efact,nz)
        call reallocate1d_r(bfact,nz)
        do k=1,nan
            do j=1,nan
                do iazi=1,nazimul
                    do iang=1,nangmul
                        SUMA(IANG,iazi,J,K)=0.0
                        SUMD(IANG,iazi,J,K)=0.0
                        SUME(IANG,iazi,J,K)=0.0
                        SUMB(IANG,iazi,J,K)=0.0
                    end do
                end do
            end do
        end do
!    2 CONTINUE
        DO I=1,NZ
            Z=I*ZSTEP+ZADD
            ZVALSQ(I)=Z*Z
            AFACT(I)=0.
            DFACT(I)=0.
            EFACT(I)=0.
            BFACT(I)=0.
        end do
!    1 CONTINUE                                                          MUL08230
!C                                                                       MUL08240
!C CALCULATE SUMMATION FACTORS FOR COMPLETE AND PARTIAL BEAMS            MUL08250
!C                                                                       MUL08260
        DO IZ1=1,NZ
            DO IZ2=1,NZ
!C                                                                       MUL08290
!C COMPLETE BEAM, COMPLETE SCATTERING - C-C                              MUL08300
!C                                                                       MUL08310
                IZP=IZ2-IZ1
                IZP=1+IABS(IZP)
                AFACT(IZP)=AFACT(IZP)+1.
!C                                                                       MUL08350
!C COMPLETE BEAM, PARTIAL SCATTERING - C-P                               MUL08360
!C                                                                       MUL08370
                IF(IZ2.ge.MINS.and.IZ2.le.MAXS) then
                    DFACT(IZP)=DFACT(IZP)+1.
                end if
!412   CONTINUE                                                          MUL08400
!C                                                                       MUL08410
!C PARTIAL BEAM, COMPLETE SCATTERING - P-C                               MUL08420
!C                                                                       MUL08430
                IF(IZ1.ge.MINB.and.IZ1.le.MAXB) then
                    EFACT(IZP)=EFACT(IZP)+1.
!C                                                                       MUL08460
!C PARTIAL BEAM, PARTIAL SCATTERING - P-P                                MUL08470
!C                                                                       MUL08480
                    IF(IZ2.ge.MINS.and.IZ2.le.MAXS) then
                        BFACT(IZP)=BFACT(IZP)+1.
                    end if
                end if
            end do
        end do
!        do izp=1,nz
!            write(6,*) 'summu1> ',afact(izp),dfact(izp),efact(izp),bfact(izp)
!        end do
!C                                                                       MUL08530
!C SET UP MEAN VALUES OF 1/(L**2) FOR SINGLE ELEMENTS                    MUL08540
!C                                                                       MUL08550
        NR1=0
        DO I1=1,NAN
            RSTEP=RST(I1)
            NRAD1=NRAD(I1)
            IF(NRAD1.gt.0) then
                DO J1=1,NRAD1
                    NR1=NR1+1
                    TSTEP=2*PI*R1(NR1)/NOMEGA(NR1)
!C                                                                       MUL08630
!C GENERATE AN EQUIVALENT RADIUS BASED ON THE SURFACE AREA               MUL08640
!C OF A CUBE, AND EVALUATE THE ZEROTH CONTRIBUTION                       MUL08650
!C                                                                       MUL08660
                    AREA=2*(RSTEP*TSTEP+RSTEP*ZSTEP+TSTEP*ZSTEP)
                    ZERO(NR1)=33.78/AREA
                    ARAD=2.25/ZERO(NR1)
                    EQRAD(NR1)=SQRT(ARAD)
!		write(6,*) i1,nr1,nomega(nr1),rstep,r1(nr1),tstep,zstep,zero(nr1),eqrad(nr1)
                end do
            end if
        end do
!        stop
!C                                                                       MUL08740
!C***********************************************************************MUL08750
!C                                                                       MUL08760
!C CALCULATE PRIMARY SCATTERING FOR COMPLETE AND PARTIALLY EMMERSED      MUL08770
!C BEAMS                                                                 MUL08780
!C                                                                       MUL08790
        NOM1=1
        NR1=0
        NSERCH=1
        call reallocate1d_r(areac,nazimul)
        call reallocate1d_r(areas,nazimul)
        call reallocate1d_r(areat,nazimul)
        call reallocate1d_r(areap,nazimul)
        call reallocate2d_r(prodc,nangmul,nazimul)
        call reallocate2d_r(prods,nangmul,nazimul)
        call reallocate2d_r(prodt,nangmul,nazimul)
        call reallocate2d_r(prodp,nangmul,nazimul)
        call reallocate3d_r(sumc,nangmul,nazimul,nan)
        call reallocate3d_r(sums,nangmul,nazimul,nan)
        call reallocate3d_r(sumt,nangmul,nazimul,nan)
        call reallocate3d_r(sump,nangmul,nazimul,nan)
        DO I1=1,NAN
!C                                                                       MUL08840
!C INITIALIZE THE AREA SUMS WHICH WILL BE USED TO NORMALIZE SUMC,        MUL08850
!C SUMS, SUMT,  AND SUMB                                                 MUL08860
!C                                                                       MUL08870
            DO IAZI=1,nazimul
                AREAC(IAZI)=0.                                                 
                AREAS(IAZI)=0.                                                 
                AREAT(IAZI)=0.                                                
                AREAP(IAZI)=0.                                                
                do iang=1,nangmul
                    SUMC(IANG,iazi,I1)=0.0
                    SUMS(IANG,iazi,I1)=0.0
                    SUMT(IANG,iazi,I1)=0.0
                    SUMP(IANG,iazi,I1)=0.0
                end do
            end do
            NRAD1=NRAD(I1)
            IF(NRAD1.gt.0) then
                DO J1=1,NRAD1
                    NR1=NR1+1
                    NOMEG1=NOMEGA(NR1)
                    AREA1=AREAEL(NR1)
                    DO IAZI=1,nazimul
                        do iang=1,nangmul
                            PRODC(IANG,iazi)=0.                                            
                            PRODS(IANG,iazi)=0.                                          
                            PRODT(IANG,iazi)=0.                                            
                            PRODP(IANG,iazi)=0. 
                        end do
                    end do
!C                                                                       MUL09100
!C PICK UP THE APPROPRIATE VALUES OF ALEN1 AND ALEN2                     MUL09110
!C                                                                       MUL09120
                    DO IAZI=1,nazimul
                        NDOM1=NOM1
                        DO K1=1,NOMEG1
                            RT1T=NTEST(NDOM1)
                            NT2T=NT2(NSERCH)
                            AREAC(IAZI)=AREAC(IAZI)+AREA1
                            AREAS(IAZI)=AREAS(IAZI)+AREA1*NT2T
                            AREAT(IAZI)=AREAT(IAZI)+AREA1*RT1T
                            AREAP(IAZI)=AREAP(IAZI)+AREA1*RT1T*NT2T
                            do iang=1,nangmul
                                XPROD=ALEN1(NDOM1)*ALEN2(iang,NSERCH)
                                PRODC(IANG,iazi)=PRODC(IANG,iazi)+XPROD
                                PRODS(IANG,iazi)=PRODS(IANG,iazi)+XPROD*NT2T
                                PRODT(IANG,iazi)=PRODT(IANG,iazi)+XPROD*RT1T
                                PRODP(IANG,iazi)=PRODP(IANG,iazi)+XPROD*RT1T*NT2T
                            end do
                            NDOM1=NDOM1+1
                            NSERCH=NSERCH+1
                        end do
                    end do
!c	write(6,*) ((PRODP(iang,kk),kk=1,nazimul),iang=1,nangmul)
                    NOM1=NOM1+NOMEG1
                    DO IAZI=1,nazimul
                        do iang=1,nangmul
                            PRODC(IANG,iazi)=PRODC(IANG,iazi)*AREA1                       
                            PRODS(IANG,iazi)=PRODS(IANG,iazi)*AREA1                     
                            PRODT(IANG,iazi)=PRODT(IANG,iazi)*AREA1                     
                            PRODP(IANG,iazi)=PRODP(IANG,iazi)*AREA1                     
                            SUMC(IANG,iazi,I1)=SUMC(IANG,iazi,I1)+PRODC(IANG,iazi)      
                            SUMS(IANG,iazi,I1)=SUMS(IANG,iazi,I1)+PRODS(IANG,iazi)      
                            SUMT(IANG,iazi,I1)=SUMT(IANG,iazi,I1)+PRODT(IANG,iazi)      
                            SUMP(IANG,iazi,I1)=SUMP(IANG,iazi,I1)+PRODP(IANG,iazi)      
                        end do
                    end do
                end do
!    7 CONTINUE                                                          
                do iang=1,nangmul
                    DO IAZI=1,nazimul
                        SUMC(IANG,iazi,I1)=SUMC(IANG,iazi,I1)/AREAC(IAZI)          
                        IF(AREAS(IAZI).GT.ALIMIT) SUMS(IANG,iazi,I1)=SUMS(IANG,iazi,I1)/AREAS(IAZI)                 
                        IF(AREAT(IAZI).GT.ALIMIT) SUMT(IANG,iazi,I1)=SUMT(IANG,iazi,I1)/AREAT(IAZI)                 
                        IF(AREAP(IAZI).GT.ALIMIT) SUMP(IANG,iazi,I1)=SUMP(IANG,iazi,I1)/AREAP(IAZI)                
                    end do
                end do
            end if
!            write(6,*) 'summu1> areac(1),areas(1),areat(1),areap(1)',areac(1),areas(1),areat(1),areap(1)
!            write(6,*) 'summu1> ian,sumc(1,1,i1),sums(1,1,i1),sumt(1,1,i1),sump(1,1,i1)',i1,sumc(1,1,i1),sums(1,1,i1),sumt(1,1,i1) &
!            ,sump(1,1,i1)
        end do
!C                                                                       MUL09550
!C***********************************************************************MUL09560
!C                                                                       MUL09570
!C  CALCULATE SECONDARY SCATTERING                                       MUL09580
!C                                                                       MUL09590
!C DEFINE R1 VALUES                                                      MUL09600
!C                                                                       MUL09610
        call reallocate2d_r(areaa,nazimul,2)
        call reallocate2d_r(aread,nazimul,2)
        call reallocate2d_r(areae,nazimul,2)
        call reallocate2d_r(areab,nazimul,2)
        call reallocate3d_r(sumha,nangmul,nazimul,2)
        call reallocate3d_r(sumhd,nangmul,nazimul,2)
        call reallocate3d_r(sumhe,nangmul,nazimul,2)
        call reallocate3d_r(sumhb,nangmul,nazimul,2)
        call reallocate1d_r(tema1,nangmul)
        call reallocate1d_r(tema2,nangmul)
        call reallocate1d_r(temd1,nangmul)
        call reallocate1d_r(temd2,nangmul)
        call reallocate1d_r(teme1,nangmul)
        call reallocate1d_r(teme2,nangmul)
        call reallocate1d_r(temb1,nangmul)
        call reallocate1d_r(temb2,nangmul)
        call reallocate4d_r(suma,nangmul,nazimul,nan,nan)
        call reallocate4d_r(sumd,nangmul,nazimul,nan,nan)
        call reallocate4d_r(sume,nangmul,nazimul,nan,nan)
        call reallocate4d_r(sumb,nangmul,nazimul,nan,nan)
        call reallocate1d_r(storsm,nomegamax)
        call reallocate1d_r(storsq,nomegamax)

        
        NR1=0
        NOM1=1
        DO I1=1,NAN
            NRAD1=NRAD(I1)
!C                                                                       MUL09660
!C AVOID INTEGRATING OVER ANNULI WHICH DO NOT SCATTER                    MUL09670
!C                                                                       MUL09680
            IF(NRAD1.gt.0) then
                NR1S=NR1
                NOM1S=NOM1
                NR2=NR1
                NOM2=NOM1
                DO I2=I1,NAN
                    NRAD2=NRAD(I2)
!C                                                                       MUL09760
!C DO NOT BOTHER WITH ELEMENTS WHICH DO NOT SCATTER                      MUL09770
!C                                                                       MUL09780
                    IF(NRAD2.gt.0) then
!C                                                                       MUL09800
!C SAVE VALUE OF NR2,NOM2                                                MUL09810
!C                                                                       MUL09820
                        NR2S=NR2
                        NOM2S=NOM2
!C                                                                       MUL09850
!C RESET NR1                                                             MUL09860
!C                                                                       MUL09870
                        NR1=NR1S
                        
                        NOM1=NOM1S
!C                                                                       MUL09900
!C SET AREA PRODUCTS AND SUMH'S TO ZERO                                  MUL09910
!C                                                                       MUL09920
                        DO J=1,2
                            DO IAZI=1,nazimul 
                                AREAA(IAZI,J)=0.
                                AREAD(IAZI,J)=0.
                                AREAE(IAZI,J)=0.
                                AREAB(IAZI,J)=0.
                                do iang=1,nangmul
                                    SUMHA(IANG,iazi,J)=0.                                         
                                    SUMHD(IANG,iazi,J)=0.                                          
                                    SUMHE(IANG,iazi,J)=0.                                          
                                    SUMHB(IANG,iazi,J)=0.                                          
                                end do
                            end do
                        end do
                        DO J1=1,NRAD1
                            NR1=NR1+1
                            R1T=R1(NR1)
                            R1SQ=RSQ(NR1)
                            R12=2*R1T
                            NOMEG1=NOMEGA(NR1)
                            AREA1=AREAEL(NR1)
!C                                                                       MUL10120
!C SET A TEMPORARY VALUE OF NOM1                                         MUL10130
!C                                                                       MUL10140
                            NTEMP1=NOM1
!C                                                                       MUL10160
!C DEFINE R2 VALUES  THESE ARE ALWAYS THE SAME OR GREATER THEN          MUL10170
!C R1 VALUES                                                             MUL10180
!C                                                                       MUL10190
                            NR2=NR2S
                            NOM2=NOM2S
!C                                                                       MUL10220
!C IF 2ND RING IS IN SAME ANNULUS AS FIRST THEN WE HAVE ENSURED          MUL10230
!C THAT THE INITIAL R2 VAUE IS NOT LESS THAN R1 VALUE.  THIS             MUL10240
!C WILL ALSO REDUCE THE NUMBER OF 2ND RING S TO BE SUMMED OVER           MUL10250
!C                                                                       MUL10260
                            J2STAR=1
                            IF(I1.eq.I2) then
                                NR2=NR1
                                NOM2=NOM1
                                J2STAR=J1
                            end if
                            DO J2=J2STAR,NRAD2
                                R2T=R1(NR2)
                                R2SQ=RSQ(NR2)
                                ADD=R1SQ+R2SQ 
                                NOMEG2=NOMEGA(NR2)
!C                                                                       MUL10370
!C TAKE RATIO OF NOMEG2 TO NOMEG1 - IF THE PROBLEM WAS SET UP CORRECTLY  MUL10380
!C IN SUB LEN3, THIS SHOULD ALWAYS BE AN INTEGER.                        MUL10390
!C                                                                       MUL10400
                                NRAT=NOMEG2/NOMEG1
                                AREA2=AREAEL(NR2)
                                ARPROD=AREA1*AREA2
!C                                                                       MUL10440
!C DEFINE OMEGA VALUES FOR FIRST ELEMENT                                 MUL10450
!C                                                                       MUL10460
                                OM1=OM(NTEMP1)
!C                                                                       MUL10480
!C STEP THROUGH ELEMENTS IN 2ND RING AND CalCULATE HORIZONTAL            MUL10490
!C SEPARATION OF ELEMENTS AND ABSORPTION PATH LENGTH.  SINCE THE         MUL10500
!C THE GEOMETRY IS SYMMETRIC ABOUT OMEGA=0, THIS NEED ONLY BE DONE       MUL10510
!C FOR HALF THE RING                                                     MUL10520
!C                                                                       MUL10530
                                NLAST=NOMEG2/2
                                NLAST1=NLAST+1
!C                                                                       MUL10560
!C FOR THE CASES OMDIF=0 AND OMDIF=PI, THE CALCULATION OF ADDA AND       MUL10570
!C SUML DOES NOT INVOLVE CALCULATION OF SINES AND COSINES                MUL10580
!C                                                                       MUL10590
                                RADD=R1T+R2T
                                RDIFF=R2T-R1T
                                STORSQ(1)=RDIFF*RDIFF
!                                write(6,*) nr1,nr2,ntemp1,nlast1,nomegamax
                                STORSQ(NLAST1)=RADD*RADD
                                SUMR1=0.0
                                RADST1=R1T-RAD1(1)                                                 
                                RADNX1=0.0
!C                                                                       MUL10680
!C THE METHOD IS TO FORM TH SUMMATION OF PRODUCTS                        MUL10690
!C MUS(IAN)*L FROM THE CENTRE TO R1 (SUMR1) AND FROM THE CENTRE          MUL10700
!C TO R2 (SUMR2).  FOR OMDIF=0 SUML=SUMR2-SUMR1, AND FOR OMDIF=PI        MUL10710
!C SUML=SUMR1+SUMR2                                                      MUL10720
!C                                                                       MUL10730
                                ian=0
                                do while (ian.lt.nan.and.radnx1.ge.0.0)
                                    ian=ian+1
                                    AMU=MUS(IAN)
                                    SUMR1=SUMR1+RADST1*AMU
                                    RADNX1=R1T-RAD2(IAN)        
                                    IF(RADNX1.ge.0.0) then
                                        SUMR1=SUMR1-RADNX1*AMU
                                        RADST1=RADNX1
                                    end if
!                                    write(6,*) nr1,nr2,ian,sumr1,radst1,radnx1
                                end do
                                SUMR2=0.0
                                RADST2=R2T-RAD1(1)                                                 
                                RADNX2=0.0 
                                ian=0
                                do while (ian.lt.nan.and.radnx2.ge.0.0)
                                    ian=ian+1
                                    AMU=MUS(IAN)
                                    SUMR2=SUMR2+RADST2*AMU
                                    RADNX2=R2T-RAD2(IAN)                                               
                                    IF(RADNX2.ge.0.0) then
                                        SUMR2=SUMR2-RADNX2*AMU
                                        RADST2=RADNX2
                                    end if
!                                    write(6,*) nr1,nr2,ian,sumr2,radst2,radnx1
                                end do
!C                                                                       MUL10920
!C NORMALIZE THE SUMS TO THE PATH LENGTH IN EACH CASE                    MUL10930
!C                                                                       MUL10940
                                STORSM(1)=MUS(I1)
                                IF(NR1.ne.NR2) then
                                    STORSM(1)=(SUMR2-SUMR1)/RDIFF
                                end if
                                STORSM(NLAST1)=(SUMR1+SUMR2)/RADD
!                                write(6,*) nr1,nr2,nlast1,sumr1,sumr2,rdiff,radd,storsm(nlast1)
!C                                                                       MUL10990
!C SET A TEMPORARY VALUE OF NOM2. THIS IS ONE GREATER THAN THE           MUL11000
!C CURRENT VALUE OF NOM2, WHICH CORRESPONDS TO OMDIF=0 AND WHICH HAS     MUL11010
!C ALREADY BEEN SUMMED FOR.                                              MUL11020
!C (EXCUSE ENGLISH_)                                                     MUL11030
!C                                                                       MUL11040
                                NTEMP2=NOM2+1
                                DO K2=2,NLAST
                                    OM2=OM(NTEMP2)
                                    OMDIF=OM2-OM1
!C
!C CALCULATE HORIZONTAL SEPARATION OF ELEMENTS                            MUL11090
!C                                                                       MUL11100
                                    OMCOS=R2T*COS(OMDIF)
                                    ADDA=ADD-R12*OMCOS
                                    XADD=SQRT(ADDA)
!C                                                                       MUL11140
!C CALCULATE THE PERPENDICULAR DISTANCE OF THE LINE JOINING              MUL11150
!C (R2,OMDIF) AND (R1,0)                                                 MUL11160
!C FROM THE CNETRE                                                       MUL11170
!C                                                                       MUL11180
                                    SUML=0
                                    D=R1T*R2T*SIN(OMDIF)/XADD 
                                    E=R1T-OMCOS
                                    SIGCOS=SIGN(1.,E)
!C                                                                       MUL11230
!C CALCULATE THE SEPARATION OF THE POINTS - IF EITHER OR BOTH             MUL11240
!C LIE OUTSIDE A CIRCLE OF RADIUS  RAD(J), THEN THE LENGHT               MUL11250
!C OF THEIR SEPARAITON WITHIN THAT CIRCLE OS CALCULATED                  MUL11260
!C                                                                       MUL11270
                                    DO IAN=1,NAN
!C                                                                       MUL11310
!C DEFINE LENGTH FO PATH BETWEEN THE TWO POINTS WHICH LIES WITHIN        MUL11320
!C THE CIRCLE FO RAIUS RAD(J)                                            MUL11330
!C                                                                       MUL11340
                                        XSTART=DIST1(R1T,R2T,RAD1(ian),D,SIGCOS)
                                        XNEXT=DIST1(R1T,R2T,RAD2(ian),D,SIGCOS)
!C                                                                       MUL11360
!C THE PATH LENGTH THROUGH ITH ANNULUS IS NOW XNEXT-XSTART  --           MUL11370
!C XSTART WAS THE PATH THROUGH THE NEXT SMALLER CIRCLE                   MUL11380
!C                                                                       MUL11390
                                        SUML=SUML+MUS(IAN)*(XNEXT-XSTART)
                                    end do
!C                                                                       MUL11430
!C CONVERT SUML TO A FRACTIONAL PATH LENGTH OF XADD                      MUL11440
!C                                                                       MUL11450
                                    SUML=SUML/XADD
                                    STORSQ(K2)=ADDA
                                    STORSM(K2)=SUML
                                    NTEMP2=NTEMP2+1
                                end do
!C                                                                       MUL11510
!C GENERATE REMAINING NOMEG2-NLAST ELEMENTS                              MUL11520
!C                                                                       MUL11530
                                NSAVE=NLAST1
                                NLAST=NLAST1+1
                                DO K2=NLAST,NOMEG2
                                    NSAVE=NSAVE-1
                                    STORSQ(K2)=STORSQ(NSAVE)
                                    STORSM(K2)=STORSM(NSAVE)
                                end do
!C                                                                       MUL11610
!C NOW CALCULATE  THE PRODUCTS ALEN1*ALEN2 FOR ALL                       MUL11620
!C COMBINATIONS OF THE ELEMENTS IN 1ST AND 2ND RINGS - THESE ARE         MUL11630
!C STORED IN NOMEG2 TEMPORARY SUMS TEMPSU (PARTIAL) AND                  MUL11640
!C TEMCSU (COMPLETE), WHICH MUST INITIALLY BE SET TO ZERO                MUL11650
!C                                                                       MUL11660
!C                                                                       MUL11670
!C RESET NTEMP2                                                          MUL11680
!C                                                                       MUL11690
                                NTEMP2=NOM2
!C                                                                       MUL11710
!C STEP THROUGH NOMEG2 ELEMENTS FOR SECOND RING                          MUL11720
!C                                                                       MUL11730
                                DO K2=1,NOMEG2
                                    STORE=STORSQ(K2)
                                    STORM=STORSM(K2)
!C                                                                       MUL11770
!C GENERATE THE Z DEPENDENT ARRAY ZMULL, WHICH INDEPENDENT OF            MUL11780
!C K1 AND IANG                                                           MUL11790
!C                                                                       MUL11800
                                    IF(NZ.gt.1) then
                                        DO IZ2=2,NZ
                                            ZSQ=STORE+ZVALSQ(IZ2)
                                            ZSQRT=SQRT(ZSQ)
                                            ZMUL=ZSQRT*STORM
                                            ZSQ=AVL2(EQRAD(NR1),EQRAD(NR2),ZSQRT)/ZSQ 
                                            ZMULL(IZ2)=EXP(-ZMUL)*ZSQ
                                        end do
                                    end if
!C                                                                       MUL11890
!C FOR Z=0 WE USE ZERO FOR 1/L**2 FOR ELEMENTS WHOSE                     MUL11900
!C SEPARATION IS ZERO.                                                   MUL11910
!C                                                                       MUL11920
                                    IZ2=1
                                    IF(STORE.ge.ALIMIT) then
                                        ZSQRT=SQRT(STORE)
                                        ZMUL=ZSQRT*STORM
                                        ZSQ=AVL2(EQRAD(NR1),EQRAD(NR2),ZSQRT)/STORE
                                    else
                                        ZSQ=ZERO(NR1)
                                        ZMUL=MUS(I1)/SQRT(ZSQ)
                                    end if
                                    ZMULL(IZ2)=EXP(-ZMUL)*ZSQ
                                    ZMULA=0.
                                    ZMULD=0.
                                    ZMULE=0.
                                    ZMULB=0.
!                                    write(6,*) nrad(1),eqrad(nr1),eqrad(nr2)
                                    DO IZ2=1,NZ
!                                        write(6,*) iz2,zmull(iz2)
                                        ZMULA=ZMULA+ZMULL(IZ2)*AFACT(IZ2)
                                        ZMULD=ZMULD+ZMULL(IZ2)*DFACT(IZ2)
                                        ZMULE=ZMULE+ZMULL(IZ2)*EFACT(IZ2)
                                        ZMULB=ZMULB+ZMULL(IZ2)*BFACT(IZ2)
                                    end do
!                                    stop
!C                                                                       MUL12130
!C STEP THROUGH THE SCATTERING ANGLES                                    MUL12140
!C                                                                       MUL12150
                                    DO IAZI=1,nazimul
!C                                                                       MUL12170
!C RESET NTEMP1, AND DEFINE CONSTANTS FOR REFERENCE TO ARRAY             MUL12180
!C ALEN2                                                                 MUL12190
!C                                                                       MUL12200
                                        NTEMP1=NOM1
                                        NALEN1=(nazimul-1)*(NOM1-1)+(IAZI-1)*NOMEG1                       
                                        NALEN2=(nazimul-1)*(NOM2-1)+(IAZI-1)*NOMEG2                       
                                        do iang=1,nangmul
                                            TEMA1(iang)=0.0
                                            TEMA2(iang)=0.0
                                            TEMD1(iang)=0.0
                                            TEMD2(iang)=0.0
                                            TEME1(iang)=0.0
                                            TEME2(iang)=0.0
                                            TEMB1(iang)=0.0
                                            TEMB2(iang)=0.0
                                        END DO
!C                                                                       MUL12320
!C STEP THROUGH NOMEG1 ELEMENTS IN FIRST RING                            MUL12330
!C                                                                       MUL12340
                                        DO K1=1,NOMEG1
!C THE VALUE OF NTEMP2 MUST REMAIN BETWEEN THE VALUES OF NOM2            MUL12360
!C AND NOM2+NOMEG2-1 , IN ORDER FOR IT TO ALWAYS REPRESENT ELEMENTS IN THMUL12370
!C 2ND RING                                                              MUL12380
!                                                                       MUL12390
                                            NA=NTEMP2-NOM2
                                            NA=NA/NOMEG2
                                            NTEMP2=NTEMP2-NA*NOMEG2
!C                                                                       MUL12430
!C DEFINE THE REFERENCE  INTEGERS FOR ALEN2 FOR 1ST AND 2ND RINGS        MUL12440
!C                                                                       MUL12450
                                            NSERC1=NALEN1+NTEMP1
                                            NSERC2=NALEN2+NTEMP2
!C                                                                       MUL12480
!C DEFINE PROBABILITIES OF ELEMENT BEING IN IN BEAM AND OF BEING SEEN    MUL12490
!C BY DETECTOR                                                           MUL12500
!C                                                                       MUL12510
                                            RT1T1=NTEST(NTEMP1)
                                            RT1T2=NTEST(NTEMP2)
                                            NT2T1=NT2(NSERC1)
                                            NT2T2=NT2(NSERC2)
!C                                                                       MUL12560
!C PRODUCTS WITH FIRST ELEMENT AS PRIMARY                                MUL12570
!C                                                                       MUL12580
                                            do iang=1,nangmul
                                                PROD1A=ALEN1(NTEMP1)*ALEN2(iang,NSERC2)
                                                PROD1D=PROD1A*NT2T2
                                                PROD1E=PROD1A*RT1T1
                                                PROD1B=PROD1A*NT2T2*RT1T1
!C                                                                       MUL12630
!C ADD TO TEMPORARY SUMS                                                 MUL12640
!C                                                                       MUL12650
                                                TEMA1(iang)=TEMA1(iang)+PROD1A                                 
                                                TEMD1(iang)=TEMD1(iang)+PROD1D                                 
                                                TEME1(iang)=TEME1(iang)+PROD1E                                  
                                                TEMB1(iang)=TEMB1(iang)+PROD1B                                 
                                            end do
!C                                                                       MUL12700
!C ADD ARPROD TO APPROPRIATE AREA PRODUCT SUM                            MUL12710
!C                                                                           MUL12720
                                            AREAA(IAZI,1)=AREAA(IAZI,1)+ARPROD
                                            AREAD(IAZI,1)=AREAD(IAZI,1)+ARPROD*NT2T2
                                            AREAE(IAZI,1)=AREAE(IAZI,1)+ARPROD*RT1T1
                                            AREAB(IAZI,1)=AREAB(IAZI,1)+ARPROD*NT2T2*RT1T1
!C                                                                       MUL12770
!C IF NR1.EQ.NR2 THERE IS NO NEED TO ADD TO SECOND TEMPORARY             MUL12780
!C SUMS AS WELL SINCE THESE ADDITIONS WOULD BE INDISINGIUSHABLE           MUL12790
!C FROM THOSE TO THE FIRST SUMS                                          MUL12800
!C                                                                       MUL12810
                                            IF(NR1.ne.NR2) then
!C                                                                       MUL12830
!C PRODUCTS WITH SECOND ELEMENT AS PRIMARY                               MUL12840
!C                                                                       MUL12850
                                                do iang=1,nangmul
                                                    PROD2A=ALEN1(NTEMP2)*ALEN2(iang,NSERC1)
                                                    PROD2D=PROD2A*NT2T1
                                                    PROD2E=PROD2A*RT1T2
                                                    PROD2B=PROD2A*NT2T1*RT1T2
                                                    TEMA2(iang)=TEMA2(iang)+PROD2A
                                                    TEMD2(iang)=TEMD2(iang)+PROD2D
                                                    TEME2(iang)=TEME2(iang)+PROD2E
                                                    TEMB2(iang)=TEMB2(iang)+PROD2B
                                                end do
!C                                                                       MUL12940
!C ADD ARPROD TO APPROPRIATE AREA PRODUCTS                               MUL12950
!C                                                                       MUL12960
                                                AREAA(IAZI,2)=AREAA(IAZI,2)+ARPROD
                                                AREAD(IAZI,2)=AREAD(IAZI,2)+ARPROD*NT2T1
                                                AREAE(IAZI,2)=AREAE(IAZI,2)+ARPROD*RT1T2
                                                AREAB(IAZI,2)=AREAB(IAZI,2)+ARPROD*NT2T1*RT1T2
                                            end if
!C                                                                       MUL13020
!C INCREMENT NTEMP1 AND INCREMENT NTEMP2 SO THAT THE CORRESPONDING       MUL13030
!C VALUE OF OM2 STARTS AT THE SAME VALUE AS OM1                          MUL13040
!C                                                                       MUL13050
                                            NTEMP1=NTEMP1+1
                                            NTEMP2=NTEMP2+NRAT
!  130 CONTINUE
                                        end do
!C MULTIPLY THE TEMPORARY SUMS BY THE AREA PRODUCTS                      MUL13100
!C AND SUM OVER Z VALUES                                                 MUL13110
!C                                                                       MUL13120
                                        do iang=1,nangmul
                                            TEMA1(IANG)=TEMA1(IANG)*ARPROD                                 
                                            TEMA2(IANG)=TEMA2(IANG)*ARPROD                                  
                                            TEMD1(IANG)=TEMD1(IANG)*ARPROD                                  
                                            TEMD2(IANG)=TEMD2(IANG)*ARPROD                                 
                                            TEME1(IANG)=TEME1(IANG)*ARPROD                                 
                                            TEME2(IANG)=TEME2(IANG)*ARPROD                                 
                                            TEMB1(IANG)=TEMB1(IANG)*ARPROD                                  
                                            TEMB2(IANG)=TEMB2(IANG)*ARPROD                                  
!		write(6,*) nr1,nr2,iazi,iang,tema1(iang),tema2(iang),temb1(iang),temb2(iang)
                                            SUMHA(IANG,IAZI,1)=SUMHA(IANG,IAZI,1)+ZMULA*TEMA1(IANG)        
                                            SUMHD(IANG,IAZI,1)=SUMHD(IANG,IAZI,1)+ZMULD*TEMD1(IANG)        
                                            SUMHE(IANG,IAZI,1)=SUMHE(IANG,IAZI,1)+ZMULE*TEME1(IANG)        
                                            SUMHB(IANG,IAZI,1)=SUMHB(IANG,IAZI,1)+ZMULB*TEMB1(IANG)        
                                            SUMHA(IANG,IAZI,2)=SUMHA(IANG,IAZI,2)+ZMULA*TEMA2(IANG)        
                                            SUMHD(IANG,IAZI,2)=SUMHD(IANG,IAZI,2)+ZMULD*TEMD2(IANG)       
                                            SUMHE(IANG,IAZI,2)=SUMHE(IANG,IAZI,2)+ZMULE*TEME2(IANG)       
                                            SUMHB(IANG,IAZI,2)=SUMHB(IANG,IAZI,2)+ZMULB*TEMB2(IANG)       
!		write(6,*) nr1,nr2,iazi,iang,SUMHA(IANG,IAZI,1),SUMHA(IANG,IAZI,2),SUMHB(IANG,IAZI,1),SUMHB(IANG,IAZI,2)
                                        end do
!  200 CONTINUE
                                    end do
!C                                                                       MUL13300
!C INCREMENT NTEMP2                                                      MUL13310
!C                                                                       MUL13320
                                    NTEMP2=NTEMP2+1
!  140 CONTINUE                                                          MUL13340
                                end do
!C                                                                       MUL13350
!C HAVING COMPLETED ALL INTEGRATIONS FOR THE CURRENT 2ND                 MUL13360
!C ELEMENT, WE CAN LOOK AT NEXT 2ND AND INCREMENT NR2 AND NOM2           MUL13370
!C                                                                       MUL13380
                                NR2=NR2+1
                                NOM2=NOM2+NOMEG2
!   40 CONTINUE
                            end do
!C                                                                       MUL13420
!C HAVING NOW EXHAUSTED ALL VALUES OF R2 FOR THIS VALUE OF R1            MUL13430
!C WE CAN INCREMENT NOM1 AND PROCEED TO A NEW VALUE OF R1                MUL13440
!C                                                                       MUL13450
                            NOM1=NOM1+NOMEG1
!   20 CONTINUE
                        end do
!C ADD SUMH>S                                                            MUL13480
!C                                                                       MUL13490
!                        write(6,*) 'summu1> ',areaa(1,1),areaa(1,2),areab(1,1),areab(1,2)
!                        write(6,*) 'summu1> ',zmula,zmuld,zmule,zmulb
                        do iang=1,nangmul
!		write(6,*) tema1(iang),tema2(iang),temb1(iang),temb2(iang)
!                write(6,*) highcc,highcp,highpc,highpp
                            DO IAZI=1,nazimul
                                SUMA1=SUMHA(IANG,IAZI,1)*ZSTEPS/HIGHCC                           
                                SUMA2=SUMHA(IANG,IAZI,2)*ZSTEPS/HIGHCC                           
                                SUMD1=SUMHD(IANG,IAZI,1)*ZSTEPS/HIGHCP                           
                                SUMD2=SUMHD(IANG,IAZI,2)*ZSTEPS/HIGHCP                           
                                SUME1=SUMHE(IANG,IAZI,1)*ZSTEPS/HIGHPC                           
                                SUME2=SUMHE(IANG,IAZI,2)*ZSTEPS/HIGHPC                           
                                SUMB1=SUMHB(IANG,IAZI,1)*ZSTEPS/HIGHPP                           
                                SUMB2=SUMHB(IANG,IAZI,2)*ZSTEPS/HIGHPP                            
 !  		write(6,*) iang,iazi,i1,i2,suma1,suma2,sumb1,sumb2
500                   FORMAT(1X,5(2X,E13.6))
                                IF(I1.ne.I2) then
                                    IF(AREAa(IAZI,1).ne.0.) then
                                        SUMA(IANG,IAZI,I1,I2)=SUMA(IANG,IAZI,I1,I2)+SUMA1/AREAA(IAZI,1)
                                    end if
                                    if(areaa(IAZI,2).ne.0.) then
                                        SUMA(IANG,IAZI,I2,I1)=SUMA(IANG,IAZI,I2,I1)+SUMA2/AREAA(IAZI,2)
                                    end if
                                    IF(AREAD(IAZI,1).ne.0.) then
                                        SUMD(IANG,IAZI,I1,I2)=SUMD(IANG,IAZI,I1,I2)+SUMD1/AREAD(IAZI,1)
                                    end if
                                    IF(AREAD(IAZI,2).ne.0.) then
                                        SUMD(IANG,IAZI,I2,I1)=SUMD(IANG,IAZI,I2,I1)+SUMD2/AREAD(IAZI,2)
                                    end if
                                    IF(AREAE(IAZI,1).ne.0.) then
                                        SUME(IANG,IAZI,I1,I2)=SUME(IANG,IAZI,I1,I2)+SUME1/AREAE(IAZI,1)
                                    end if
                                    if(AREAE(IAZI,2).ne.0.) then
                                        SUME(IANG,IAZI,I2,I1)=SUME(IANG,IAZI,I2,I1)+SUME2/AREAE(IAZI,2)
                                    end if
                                    IF(AREAB(IAZI,1).ne.0.) then 
                                        SUMB(IANG,IAZI,I1,I2)=SUMB(IANG,IAZI,I1,I2)+SUMB1/AREAB(IAZI,1)
                                    end if
                                    IF(AREAB(IAZI,2).ne.0.) then
                                        SUMB(IANG,IAZI,I2,I1)=SUMB(IANG,IAZI,I2,I1)+SUMB2/AREAB(IAZI,2)
                                    end if
                                else
                                    SUMAR=AREAA(IAZI,1)+AREAA(IAZI,2)                             
                                    IF(SUMAR.NE.0.) then
                                        SUMA(IANG,IAZI,I1,I2)=(SUMA1+SUMA2)/sumar
                                    end if
                                    SUMAR=AREAD(IAZI,1)+AREAD(IAZI,2)                             
                                    IF(SUMAR.NE.0.) then
                                        SUMD(IANG,IAZI,I1,I2)=(SUMD1+SUMD2)/SUMAR                  
                                    end if
                                    SUMAR=AREAE(IAZI,1)+AREAE(IAZI,2)                              
                                    IF(SUMAR.NE.0.) then
                                        SUME(IANG,IAZI,I1,I2)=(SUME1+SUME2)/SUMAR                    	
                                    end if
                                    SUMAR=AREAB(IAZI,1)+AREAB(IAZI,2)                            
                                    IF(SUMAR.NE.0.) then
                                        SUMB(IANG,IAZI,I1,I2)=(SUMB1+SUMB2)/SUMAR
                                    end if
                                endif
                            end do
                        end do
                    end if
                end do !End of second annulus
            end if
        end do !End of first annulus
!        stop
        RETURN
    END subroutine summu1
    
!C                                                                       MUL06340
!C       ***********************                                         MUL06350
!C       *                     *                                         MUL06360
!C       *  SUBROUTINE DATAOP  *                                         MUL06370
!C       *                     *                                         MUL06380
!C       ***********************                                         MUL06390
!C                                                                       MUL06400
    SUBROUTINE DATAOP(il)            

	use beam_routines
        use mul_corr_azi
        
        integer                                :: il
        
        real, dimension(:), allocatable         :: TSUMC,TSUMS,TSUMT,TSUMP,DIVCC,DIVCP!(mcorrang)
        real, dimension(:), allocatable         :: AREAA,AREAB,AREAS,ARBS
        real, dimension(:), allocatable         :: T
        real, dimension(:), allocatable         :: TSUMA,TSUMD,TSUME,TSUMB
        real, dimension(:,:), allocatable       :: AREABS!(mcorrang,mcont)
        
        real pi,fact,sumpri,thirsc,tomulsca,tosca,arproa,arprod,arproe,arprob
        integer i,kk,iang,iazi,ian,i1,i2
        
        PI=4.0*atan(1.0)
!        write(6,*) 'dataop> HIGHT,highb,highs,highbs: ',HIGHT,highb,highs,highbs
!C                                                                       MUL06530
!C GENERATE THE EQUIVALENT SAMPLE AREAS IN THE CALCULATION - THIS        MUL06540
!C MUST BE DONE NUMERICALLY SO AS TO INCORPORATE THE BEAM PROFILE 
!C                                                                       MUL06560
        call reallocate1d_r(areaa,nan)
        call reallocate1d_r(areab,nan)
        call reallocate1d_r(areas,nan)
        call reallocate2d_r(areabs,nazimul,nan)
      
        DO I=1,NAN
            FACT=SIGS(I)*DEN(I)/(4*PI)
            IF(FACT.GT.0.) then
                CALL EQUIV(RAD1(I),RAD2(I),100,AREAA(I),AREAB(I),AREAS(I),ARBS)
                AREAA(I)=AREAA(I)*FACT
                AREAB(I)=AREAB(I)*FACT
                AREAS(I)=AREAS(I)*FACT
                DO KK=1,nazimul
                    AREABS(KK,I)=ARBS(KK)*FACT
                end do
            else
                areaa(i)=0.0
                areab(i)=0.0
                areas(i)=0.0
                DO KK=1,nazimul
                    AREABS(KK,I)=0.0
                end do
            end if
!            write(6,*) 'dataop> i,areaa,areab,areas,areabs(1): ',i,areaa(i),areab(i),areas(i),areabs(1,i)
!            write(6,*) 'dataop> i,sumc(1,1,i),sums(1,1,i),sumt(1,1,i),sump(1,1,i)',i,sumc(1,1,i),sums(1,1,i),sumt(1,1,i) &
 !           ,sump(1,1,i)
        end do
!c	write(6,*) (areaBS(kk,1),kk=1,nazimul)
!c	write(6,*) (sumP1(1,kk,1),kk=1,nazimul)
!c	write(6,*) (sumB1(1,kk,1,1),kk=1,nazimul)
!C                                                                       MUL06680
!C STEP THROUGH THE ANGLES
!C
	call reallocate1d_r(tsumc,nazimul)
        call reallocate1d_r(tsums,nazimul)
        call reallocate1d_r(tsumt,nazimul)
        call reallocate1d_r(tsump,nazimul)
	call reallocate1d_r(tsuma,nazimul)
        call reallocate1d_r(tsumd,nazimul)
        call reallocate1d_r(tsume,nazimul)
        call reallocate1d_r(tsumb,nazimul)
        call reallocate1d_r(divcc,nazimul)
        call reallocate1d_r(divcp,nazimul)
        call reallocate1d_r(tsumt,nazimul)
        call reallocate1d_r(tsump,nazimul)
        
        DO IANG=1,NANGMUL
!C                                                                       MUL06690
!C A) PRIMARY SCATTERING.                                                MUL06700
!C                                                                       MUL06710
            DO IAZI=1,nazimul
                TSUMC(IAZI)=0.
                TSUMS(IAZI)=0.
                TSUMT(IAZI)=0.
                TSUMP(IAZI)=0.
                DO IAN=1,NAN
                    TSUMC(IAZI)=TSUMC(IAZI)+SUMC(IANG,IAZI,IAN)*HIGHT*AREAA(IAN)
                    TSUMS(IAZI)=TSUMS(IAZI)+SUMS(IANG,IAZI,IAN)*HIGHS*AREAS(IAN)
                    TSUMT(IAZI)=TSUMT(IAZI)+SUMT(IANG,IAZI,IAN)*HIGHB*AREAB(IAN)
                    TSUMP(IAZI)=TSUMP(IAZI)+SUMP(IANG,IAZI,IAN)*HIGHBS*AREABS(IAZI,IAN)
                end do
            end do
!            write(6,*) 'dataop> tsumc(1),tsums(1),tsumt(1),tsump(1) ',tsumc(1),tsums(1),tsumt(1),tsump(1)
!C                                                                       MUL06840
!C B) SECONDARY SCATTERING.                                              MUL06850
!C                                                                       MUL06860
            DO IAZI=1,nazimul
                TSUMA(IAZI)=0.
                TSUMD(IAZI)=0.
                TSUME(IAZI)=0.
                TSUMB(IAZI)=0.
            end do
!C                                                                       MUL06930
!C STEP THROUGH ANNULI AS PRIMARY SCATTER                                MUL06940
!C                                                                       MUL06950
            DO I1=1,NAN
!                                                                       MUL06970
!C STEP THROUGH ANNULI AS SECONDARY SCATTERER                            MUL06980
!C                                                                       MUL06990
 !               write(6,*) areaa(i1),areab(i1),areas(i1)
                DO I2=1,NAN
                    ARPROA=HIGHT*HIGHT*AREAA(I1)*AREAA(I2)
                    ARPROD=HIGHT*HIGHS*AREAA(I1)*AREAS(I2)
                    ARPROE=HIGHB*HIGHT*AREAB(I1)*AREAA(I2)
                    ARPROB=HIGHB*HIGHS*AREAB(I1)*AREAS(I2)
                    DO IAZI=1,nazimul
                        TSUMA(IAZI)=TSUMA(IAZI)+ARPROA*SUMA(IANG,IAZI,I1,I2)               
                        TSUMD(IAZI)=TSUMD(IAZI)+ARPROD*SUMD(IANG,IAZI,I1,I2)               
                        TSUME(IAZI)=TSUME(IAZI)+ARPROE*SUME(IANG,IAZI,I1,I2)               
                        TSUMB(IAZI)=TSUMB(IAZI)+ARPROB*SUMB(IANG,IAZI,I1,I2)
                    end do
                end do
            end do
!            write(6,*) arproa,arprod,arproe,arprob,(tsumb(kk),kk=1,nazimul)
!C                                                                       MUL07130
!C OUTPUT RESULTS AND FORM SOME RATIOS - DIVCP AND DIVCC                 MUL07140
!C                                                                       MUL07150
            DO I=1,nazimul
                DIVCC(I)=TSUMA(I)/TSUMC(I)
                DIVCP(I)=TSUMD(I)/TSUMS(I)
            END DO
            DO IAZI=1,nazimul
                SUMPRI=TSUMP(IAZI)
!C                                                                       MUL07260
!C CALCULATE TOTAL MULTIPLE SCATTERING                                   MUL07270
!C                                                                       MUL07280
                THIRSC=TSUME(IAZI)*DIVCP(IAZI)/(1.-DIVCC(IAZI))
                TOMULSCA=TSUMB(IAZI)+THIRSC
                TOSCA=TOMULSCA+SUMPRI
                onescat(il,iang,iazi)=SUMPRI
                mulscat(il,iang,iazi)=TOMULSCA
            end do
        END DO
        RETURN
    END subroutine dataop    
    
!C                                                                       MUL15380
!C       **********************                                          MUL15390
!C       *                    *                                          MUL15400
!C       *  SUBROUTINE EQUIV  *                                          MUL15410
!C       *                    *                                          MUL15420
!C       **********************                                          MUL15430
!C                                                                       MUL15440
    SUBROUTINE EQUIV(R1,R2,MS1,AR,ARB,ARS,ARBS)

        use beam_routines
        use mul_corr_azi
        
        real r1,r2,ar,arb,ars
        real, dimension(:), allocatable, optional            :: ARBS!(mcorrang)
        integer ms1
        real pi,piconv,pi2,rstep,radd,r,omst,arel,omeg,d,p1,p2,d1,ang
        integer ms,i,ir,io,nom
        
        PI=4.0*atan(1.0)
        piconv=pi/180.0
        PI2=2*PI
        MS=(1.-R1/R2)*MS1
        IF(MS.LT.1) MS=1
        RSTEP=(R2-R1)/MS
        RADD=R1-0.5*RSTEP
        AR=0.
        ARB=0.
        ARS=0.
        if(present(arbs)) then
           call reallocate1d_r(arbs,nazimul)
           DO I=1,nazimul
              ARBS(I)=0.
           end do
        end if
        DO IR=1,MS
            R=RADD+IR*RSTEP
            NOM=PI2*R/RSTEP
            OMST=PI2/NOM
            AREL=RSTEP*R*OMST
            DO IO=1,NOM
                OMEG=(IO-1)*OMST
                D=R*SIN(OMEG)
                P1=PROBE(D,PROFIL,NPROF,PRSTEP,A,B)
                P2=1.
                IF(D.GT.A1.OR.D.LT.B1) P2=0.
                AR=AR+AREL
                ARB=ARB+P1*AREL
                ARS=ARS+P2*AREL
                if(present(arbs)) then
                   DO I=1,nazimul
                     ang=azimul(i)*piconv
                     D1=R*SIN(OMEG-ANG)
                     P2=1.
                     IF(D1.GT.A1.OR.D1.LT.B1) P2=0.
                     ARBS(I)=ARBS(I)+P1*P2*AREL
                   end do
                end if
            end do
        end do
        RETURN
    END subroutine equiv
        
END MODULE cylmultof_routines
