      subroutine demean(a, n)
      dimension a(1)
      ave = 0.
      do 10 i = 1, n
   10 ave = ave + a(i)
      ave = ave / float(n)
      do 20 i = 1, n
   20 a(i) = a(i) - ave
      return 
      end
