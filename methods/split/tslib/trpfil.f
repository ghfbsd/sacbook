      subroutine trpfil(a, lenf, df, f0, f1, f2, f3, f4)
      complex a(1)
      i1 = ((f1 - f0) / df) + 1.5
      i2 = ((f2 - f0) / df) + 1.5
      i3 = ((f3 - f0) / df) + 1.5
      i4 = ((f4 - f0) / df) + 1.5
      write(unit=6, fmt=100) i1, i2, i3, i4
  100 format(4i10)
      f21 = float(i2 - i1)
      f43 = float(i4 - i3)
      if (f21 .gt. 0.) then
      slope1 = 1. / f21
      do 10 i = i1, i2
      fil = float(i - i1) / f21
   10 a(i) = a(i) * fil
      end if
      if (f43 .gt. 0.) then
      slope2 = 1. / f43
      do 20 i = i3, i4
      fil = 1. - (float(i - i3) / f43)
   20 a(i) = a(i) * fil
      end if
      call zero(a, 1, i1)
      call zero(a, i4, lenf)
      return 
      end
