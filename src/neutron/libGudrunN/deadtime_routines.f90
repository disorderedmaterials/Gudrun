!     
! File:   get_deadtimes.f90
! Author: aks45
!
! Created on 27 April 2012, 11:46
!

MODULE deadtime_routines

    implicit none

    real                                :: ddeadtime,mdeadtime,adeadtime    ! default detector, module and acquisition deadtimes
    real, dimension(:), allocatable     :: ddeadtimem,mdeadtimem    !detector, module, deadtimes for specific modules
    character(len=256) :: deadtimefile    !name of file containing detector deadtimes

    integer, parameter :: mdeadsolve=100,ndeadsolve=100
    real, dimension(mdeadsolve) :: mdatarate,lambdaplus,lambdaminus,gradplus,gradminus

CONTAINS

    subroutine get_deadtimes()
         
        use reallocation_routines
        use run_par
    
!***********************************************************************************
!*
!*    get_deadtimes.for
!*
!*    Copyright A K Soper, March 2004
!*
!*    This program reads the default detector, module and data acquisition deadtimes
!* from a file called "deadtime.cor". Individual detector and module deadtimes can be
!* specified for individual modules if required.
!*
!***********************************************************************************
!c
!c internal variables
!c
        integer :: im,it,nchard,ierr
        real :: ddeadtimer,mdeadtimer    !used as dummy variables to read deadtimes
        character(len=256) fullfilename

!c Get the solutions for the paralysable deadtime equation
!c This is not used at this time 6/3/06
!c    call solve_deadtime()
!c    do it=1,ndeadsolve
!c       write(6,*) it,mdatarate(it),lambdaminus(it),lambdaplus(it)
!c    end do
        ! Stop if no modules are defined, otherwise allocate the arrays
        if(nmod.gt.0) then
            call reallocate1d_r(ddeadtimem,nmod)
            call reallocate1d_r(mdeadtimem,nmod)
        else
            write(6,*) 'get_deadtimes> No modules defined!'
            stop
        end if
!c
!c initialise the default deadtimes
!c
        ddeadtime=0.0
        mdeadtime=0.0
        adeadtime=0.0
!c
!c counter for number of entries read from dead time file
!c
        it=0
!c
!c determine number of characters in deadtimes file - if specified
!c
        nchard=index(deadtimefile,'.')
        write(6,*) 'Deadtime filename: ',deadtimefile(1:len_trim(deadtimefile)),' ',nchard
!c
!c access the deadtimes file if it is specified.
!c
        if(nchard.gt.0) then
            open(21,file=deadtimefile,status='old',iostat=ierr)

!c If there is a problem with this file, get the revised filename, based on the 
!c default startup folder.

            if(ierr.ne.0) then
                deadtimefile=convertfilename(deadtimefile)
                open(21,file=deadtimefile,status='old',iostat=ierr)
            endif
!c
!c if ierr.ne.0 there is a problem with this file
!c
            do while(ierr.eq.0) 
!c
!c first line is different from others
!c
                if(it.eq.0) then
                    read(21,*,iostat=ierr) ddeadtime,mdeadtime,adeadtime
                    if(ierr.eq.0) then
!c
!c valid line read - increment number of good lines read
!c
                        it=it+1
!c
!c set all the detector deadtimes to these values initially
!c
                        do im=1,nmod
                            ddeadtimem(im)=ddeadtime
                            mdeadtimem(im)=mdeadtime
                        end do
                    endif
                    write(6,*) ierr,it,ddeadtime,mdeadtime,adeadtime
                else
!c
!c read separate values for specific detectors
!c
                    read(21,*,iostat=ierr) im,ddeadtimer,mdeadtimer
                    write(6,*) im,ddeadtimer,mdeadtimer,ierr
                    if(ierr.eq.0) then
!c
!c valid line read - increment number of good lines read
!c
                        it=it+1
!c
!c and set deadtime values to their default values if they are zero
!c
                        if(im.gt.0.and.im.le.nmod) then
                            if(ddeadtimer.ge.0.0) ddeadtimem(im)=ddeadtimer
                            if(mdeadtimer.ge.0.0) mdeadtimem(im)=mdeadtimer
                        endif
                    endif
                endif
            end do
            close(21)
!c
!c if a deadtime file was specified, but produced no valid entries, print out
!c a warning message and stop
!c
            if(it.eq.0) then
                write(6,*) 'get_deadtimes> ' &
                ,'Specified deadtime file either does not exist' &
                ,', or contains invalid entries. Please check'
            endif
        endif
!c
!c no dead time file specified or there is an error with the file - set the default values
!c
        if(it.eq.0) then
            do im=1,nmod
                ddeadtimem(im)=ddeadtime
                mdeadtimem(im)=mdeadtime
            end do
        endif
!c
!c write the deadtime values to the default file 'deadtime.cor'
!c
        open(21,file='deadtime.cor',status='unknown')
        write(21,*) ddeadtime,mdeadtime,adeadtime
        do im=1,nmod
            write(21,*) im,ddeadtimem(im),mdeadtimem(im)
        end do
        close(21)
        return

    end subroutine get_deadtimes

END MODULE deadtime_routines
    
