!     
! File:   fltmultof_routines.f90
! Author: aks45
!
! Created on 15 November 2013, 09:43
!

MODULE fltmultof_routines
    
    implicit none
    
    CONTAINS
    

    function prime(sec0,sec1,amu,t)
       
        real             :: prime,sec0,sec1,amu,t
        real             :: sum,amut,s,f,amuts
        integer          :: nterms
      
        sum=1.
        amut=amu*t
        s=amut*(sec0-sec1)
	if(abs(s).lt.0.01)then
!c use a simple series expansion
            nterms=5
            f=1-s/real(nterms)
            do while (nterms.gt.2)
               nterms=nterms-1
               f=1.0-s*f/real(nterms)
            end do
            sum=f               
	else
            sum=(1.0-exp(-s))/s
	endif
        if(sec1.lt.0.0) then
            prime=sum
        else
            amuts=amut*sec1
            prime=0.
            if(amuts.lt.50.) prime=sum*exp(-amut*sec1)
        endif
        return
    end function prime

    function aintex(al,sec,amu,del)
        
        real aintex,al,sec,amu,del
        real amux1,amux2,sum
        
        amux1=0.5*del*amu*abs(sec)
        amux2=al*amu*abs(sec)
        sum=0.0
        if(amux2.le.40.0) then
            sum=exp(-amux2)
            if(amux1.le.40.0) then
                sum=sum*2.0*sinh(amux1)
            end if
        endif
        aintex=sum/(sec*amu*del)
!cwrite(6,100) al,z1,z2,sec,aintex
!  100 format(1x,6e13.6)
        return
        
    end function aintex

    function binte1(amu,x0,delx)
        
        real binte1,amu,x0,delx
        real*8 sum,x1,x2,x0d
        
        x0d=abs(x0)
        x1=x0-delx
        x2=x0+delx
        sum=ainte1(amu,x1)+ainte1(amu,x2)-2.*ainte1(amu,x0d)
        binte1=sum/(delx*delx)
        return
        
    end function binte1

    function ainte1(amu,x0)
      
        real*8 ainte1,sum,damux,x0,x,damu
        real amu,amux
        
        damu=amu
        x=abs(x0)
        if(x.eq.0.0) then
            ainte1=1./(damu*damu)
            return
        endif
        amux=amu*x
        damux=amux
        sum=damux*damux*expint(damux)
        sum=sum+1.0+2.0*amux
        if(amux.lt.60.0) sum=sum+(1.0-damux)*exp(-damux)
        sum=sum*0.5/(damu*amu)
        ainte1=sum
        return
        
    end function ainte1
 
      
    function e1int1(b,al)
        
        real e1int1,b,al
        real*8 da1,dal
        real sum,a1,a2,gam,b1
        
        dal=al
        if(b.gt.1.) then
            b1=b-1
            sum=0.0
            a1=b1*al
            a2=b*al
            if(a2.lt.40.0) then
                sum=log(b1)*exp(-a2)
            endif
            sum=sum-expintiexp(a1,a2)-expint(dal)
            e1int1=-sum
            return
        else
            gam=0.5772156649
            sum=-gam-alog(al)
            a2=b*al
            sum=sum*exp(-a2)-expint(dal)
            e1int1=-sum
            return
	endif
    end function e1int1
 
    function e1int2(b,al)
        
        real e1int2,b,al
        real*8 da1,dal
        real b1,sum,a1,a2
      
        b1=b+1
        sum=alog(b1)
        a1=b1*al
        da1=a1
        dal=al
        if(a1.le.100.0) then
            sum=sum+expint(da1)
            a2=b*al
            sum=sum-exp(-a2)*expint(dal)
        end if
        e1int2=sum
        return
      
    end function e1int2


    function expint(x)
!c Requires x > 0
   
        real*8 expint,x,rx,a(4),b(4),c(5),sum1,sum2,sum,ratio,gam,an,prod
        integer n
        data a/8.5733287401,18.0590169730,8.6347608935,0.2677737343/
        data b/9.5733223454,25.6329561486,21.0996530827,3.9584969228/
        data c/0.99999193,-0.24991055,0.05519968,-0.00976004,0.00107857/
        data gam/0.5772156649/

        x=abs(x)
        expint=0.0
        if(x.le.0.0) return

        if (x.ge.1.0) then
!c Expansion from A + S p 231, 5.1.56
            rx=1.0/x
            sum1=a(4)
            sum2=b(4)
            n=4
            do while (n.gt.1)
                n=n-1
                sum1=rx*sum1+a(n)
                sum2=rx*sum2+b(n)
            end do
            sum1=1.0+sum1
            sum2=1.0+sum2
            ratio=sum1/sum2
!c Form the log of arguments to prevent overflows
            sum=x+log(x)-log(ratio)
!c      write(6,*) sum1,sum2,ratio,sum,exp(-sum)
            if(abs(sum).lt.50.0) then
                expint=exp(-sum)
            endif
        else
!c A + S p231 5.1.53
            n=5
            sum=x*c(n)
            do while (n.gt.1)
                n=n-1
                sum=x*(c(n)+sum)
            end do
            expint=sum-gam-log(x)

        endif
        return
        
    end function expint         

    function expintiexp(a1,a2)
      
        real a1,a2
        real*8 expintiexp,da1,da2,gam,sum,t,t2,prod,an
        integer n
!c calculates the value of E_i(a1)*exp(-a2) using the series expansion. Requires a1,a2 > 0.0
      
        da1=a1
        da2=a2
        gam=0.5772156649
        n=int(sqrt(a1))+1
        an=real(n)
        t=exp(-da2/an)
        prod=da1*t
        sum=prod/(an*an)
!c      write(6,*) n,prod,t,an,sum
        t2=1
        do while (n.gt.1)
            n=n-1
            an=real(n)
            t2=t2*t
            sum=prod*(t2/an+sum)/an
!c         write(6,*) n,prod,t,an,sum
        end do
        if(da2.lt.40.) then
            sum=sum+(gam+log(da1))*exp(-da2)
        endif
        expintiexp=sum
        return
        
    end function expintiexp
    
END MODULE fltmultof_routines
