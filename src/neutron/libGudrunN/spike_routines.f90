!     
! File:   find_spike.f90
! Author: aks45
!
! Created on 02 November 2013, 09:04
!

MODULE spike_routines
    
    implicit none
    
    CONTAINS
    
    subroutine find_spike()
        
        use run_par
        use local_data
        use bad_detectors
        
!***********************************************************************************
!*
!*    find_spike.FOR
!*
!*    A K Soper, October 1999
!*
!*    Searches the data channel by channel for noise spikes
!*
!***********************************************************************************

!c internal variables
!c
        integer                 :: is,isref        !internal indices
        integer                 :: j,jf,jf1,jl,jl1    !index range limits
        integer                 :: iflag1,iflag2        !flags to identify noisy channels
        integer                 :: jref1,jref2        !channel reference integers
        real                    :: rat,fir,sec,tcw1,tcw2    !temporary real values
        real                    :: old,firdev,secdev    !temporary real values
        real                    :: tsum            !temporary ratio sum
        real                    :: err            !temporary error 
        real                    :: rat1,rat2        !temporary real gradients
!c
!c ensure the requested channel range is at least between 2 and nchan-1 to avoid binning errors
!c
        if(nchfir.lt.2) nchfir=2
        if(nchlas.lt.nchfir.or.nchlas.gt.nchan-1) nchlas=nchan-1
!c
!c Start from channel 1
!c
        jf=1
        do is=nspecf,nspecl
!c
!c define final reference integer for counts array
!c
            jf=1
            jl=jf+nchan-1
            isref=is-nspecf+1
!c
!c ignore the first and last channels, which may contain binning errors
!c
            jf1=jf+nchfir-1
            jl1=jl-nchan+nchlas
            jref1=jf1-jf+1
!c
!c search for unphysical step between adjacent channels
!c
!c normalise counts to time width of channel
!c
            tcw1=tcw(jref1)
            fir=tcounts(jf1,isref)/tcw1
            firdev=fir/tcw1
            old=fir
            iflag1=0
            rat1=1
            do j=jf1+1,jl1
!c
!c define actual channel number relative to first
!c
                jref2=j-jf+1
!c
!c time channel width
!c
                tcw2=tcw(jref2)
!c
!c normalise to channel width
!c
                sec=tcounts(j,isref)/tcw2
!c
!c standard deviation of this value
!c
                secdev=sec/tcw2
!c
!c fluctuation around the mean
!c
                rat=sec-fir
                rat2=rat
!c
!c square of fluctuation around the mean
!
                err=rat*rat
!c
!c likely standard deviations of the fluctuation - use the smallest std of the two.
!
                if(firdev.le.secdev) then
                    tsum=firdev+firdev
                else
                    tsum=secdev+secdev
                endif
!c
!c The rms fluctuation should not be more than spikedev times the expected value.
!c However this test only makes sense if tsum is .gt. zero, otherwise the test is
!c indeterminate.
!c
                if(tsum.gt.0.0.and.err.gt.spikedev*tsum) then
                    iflag2=sign(1.0,rat)
                else
                    iflag2=0
                endif
!c
!c if iflag1 .ne.0 and iflag1+iflag2 = 0 then this is a spike and needs to
!c be corrected
!c
                if(iflag1.ne.0) then
!c
!c if a single step has occurred it is likely the spectrum is bad and should not
!c be used under any circumstances
!c
!c but we have to eliminate the possibility that this step is not the first
!c half of a spike, in which case the problem can be corrected
!c
                    if((iflag1+iflag2).eq.0) then
!c
!c correct the spike by taking the average of the channels immediately adjacent
!c to the bad one
!c
                        if(flagsp.ne.0) then
                            fir=0.5*(old+sec)
                            tcounts(j-1,isref)=fir*tcw1
                            iflag1=0
                            iflag2=0
                            rat=sec-fir
                            rat2=rat
                        endif
!c
!c increment the spike count - this will ensure this spectrum does not appear 
! in a group file
!c
!                        spike(is)=spike(is)+1
                    endif
!c
!c if iflag2 is 0, while iflag1 is non-zero (i.e. rat2 is outside range of 
!c statistical test) it probably indicates an unphysical step has occurred
!c hence treat this detector as bad
!c
!c                 else if(iflag2.eq.0) then
!c                    if(ibad(is).eq.0.and.flagsp.ne.0) 
!c     1ibad(is)=iflag1*runno
!c                    write(6,100) jref1,is
!c     *,runno(1:index(runno,' ')),fir,sec
!c100    format(1x,'Step found in channel ',i4,', spectrum ',i5,', File ',a
!c     1/1x,2(1x,e13.6))
!c                endif
                endif
                old=fir
                fir=sec
                firdev=secdev
                iflag1=iflag2
                rat1=rat2
                jref1=jref2
                tcw1=tcw2
            end do
            jf=jl+1
        end do
        return
        
    end subroutine find_spike

END MODULE spike_routines

