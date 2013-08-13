c
c   routine TSADD adds two time series s1,s2.
c   Result is placed in s1.
c
      subroutine tsadd(s1, s2, ns)
      real s1(*), s2(*)
c
      if (ns .lt. 1) pause '**error in TSADD - ns < 1.'
c
      do i = 1, ns
         s1(i) = s1(i) + s2(i)
      end do
      return 
      end
