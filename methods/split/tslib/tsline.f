      subroutine tsline(a, i1, i2)
      dimension a(1)
      slope = (a(i2) - a(i1)) / float(i2 - i1)
      do 10 i = i1, i2
      y = (slope * float(i - i1)) + a(i1)
      a(i) = y
   10 continue
      return 
      end
