c
c   routine TSERR calculates the L(nnorm) norm distance between sa & sb,
c   for the portion between n1 & n2. It is returned as err.
c
      subroutine tserr(sa, sb, n1, n2, nnorm, err)
c
      real sa(*), sb(*)
      snorm = 0.0
      if ((n1 .lt. 1) .or. (n2 .lt. 1)) then
      write(unit=6, fmt=*) ' error in TSERR - n1,n2 = ', n1, n2
      return 
c	
      end if
      do i = n1, n2
      s = abs(sa(i) - sb(i))
      snorm = snorm + (s ** nnorm)
c
      end do
      if (snorm .eq. 0.0) return 
      snorm = alog(snorm)
      snorm = snorm / float(nnorm)
      snorm = exp(snorm)
c
      err = snorm
      return 
      end
