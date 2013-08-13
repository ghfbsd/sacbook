c
c   routine TSMULT multiplies two time series s1,s2.
c   Result is placed in s1.
c
      subroutine tsmult(s1, s2, ns)
      real s1(*), s2(*)
c
      if (ns .lt. 1) pause '**error in TSMULT - ns < 1.'
      do i = 1, ns
	 s1(i) = s1(i) * s2(i)
      end do
      return 
      end
