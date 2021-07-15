      program test

      real a1,a2
      real*8 expint,x
      write(6,100,advance='no')
100   format(1x,'Type x: ')
      read(5,*) x
      do while(x.gt.0.0)

         write(6,*) expint(x)
         write(6,100,advance='no')
         read(5,*) x

      end do
      stop
      end

      function expint(x)
c Requires x > 0
      real*8 x,rx,a(4),b(4),c(5),expint,sum1,sum2,sum,ratio,gam,an
     1,prod
      integer n
      data a/8.5733287401,18.0590169730,8.6347608935,0.2677737343/
      data b/9.5733223454,25.6329561486,21.0996530827,3.9584969228/
      data c/0.99999193,-0.24991055,0.05519968,-0.00976004,0.00107857/
      data gam/0.5772156649/

      x=abs(x)
      rx=1.0/x
      expint=0.0
      if(x.le.0.0) return

      if (x.ge.1.0) then
c Expansion from A + S p 231, 5.1.56
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
c Form the log of arguments to prevent overflows
         sum=x+log(x)-log(ratio)
c      write(6,*) sum1,sum2,ratio,sum,exp(-sum)
         if(abs(sum).lt.50.0) then
            expint=exp(-sum)
         else
            expint=0.0
         endif

      else

c A + S p231 5.1.53
         
         n=5
         sum=x*c(n)
         do while (n.gt.1)
            n=n-1
            sum=x*(c(n)+sum)
         end do
         expint=sum-gam-log(x)

      endif

      return
      end         


