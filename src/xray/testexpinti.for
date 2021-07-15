      program test

      real a1,a2
      real*8 expintiexp
      a1=1
      do while(a1.gt.0.0)
         write(6,100,advance='no')
100   format(1x,'Type a1, a2: ')
         read(5,*) a1,a2

         write(6,*) expintiexp(a1,a2)

      end do
      stop
      end

      function expintiexp(a1,a2)
      real a1,a2
      real*8 expintiexp,da1,da2,gam,sum,t,t2,prod,an
      integer n
c calculates the value of E_i(a1)*exp(-a2) using the series expansion. Requires a1,a2 > 0.0
      da1=a1
      da2=a2
      gam=0.5772156649
      n=int(sqrt(a1))+1
      an=real(n)
      t=exp(-da2/an)
      prod=da1*t
      sum=prod/(an*an)
      write(6,*) n,prod,t,an,sum
      t2=1
      do while (n.gt.1)
         n=n-1
         an=real(n)
         t2=t2*t
         sum=prod*(t2/an+sum)/an
c         write(6,*) n,prod,t,an,sum
      end do
      if(da2.lt.40.) then
         sum=sum+(gam+log(da1))*exp(-da2)
      endif
      expintiexp=sum
      return
      end      


