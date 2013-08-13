      subroutine scdcm(a, n, a1, a2)
      dimension a(1)
      amin = a(1)
      amax = a(1)
      imin = 1
      imax = 1
      if (n .eq. 1) return 
      do 10 i = 2, n
      if (amax .ge. a(i)) goto 9
      amax = a(i)
      imax = i
    9 if (amin .le. a(i)) goto 10
      amin = a(i)
      imin = i
   10 continue
      if (imax - imin) 11, 11, 12
   11 a1 = amax
      a2 = amin
      return 
   12 a1 = amin
      a2 = amax
      return 
      end
