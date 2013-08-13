      subroutine tsnorm(ss, n1, n2, nnorm, snorm)
c
      real ss(*)
      snorm = 0.0
      if ((n1 .lt. 1) .or. (n2 .lt. 1)) then
      write(unit=6, fmt=*) ' error in TSNORM - n1,n2 = ', n1, n2
      return 
c	
      end if
      if (nnorm .lt. 0) then
      write(unit=6, fmt=*) ' Error in TSNORM - nnorm = 0 '
c
      return 
      else if (nnorm .eq. 0) then
      do i = n1, n2
      s = abs(ss(i))
      if (s .gt. snorm) snorm = s
      end do
c
      return 
      else
      do i = n1, n2
      s = abs(ss(i))
      snorm = snorm + (s ** nnorm)
c
      end do
      if (snorm .eq. 0.0) return 
      snorm = alog(snorm)
      snorm = snorm / float(nnorm)
      snorm = exp(snorm)
c
      end if
      return 
      end
