      subroutine zeroarray(a, n1, n2)
c	zeros out a real array from n1 to n2
      dimension a(1)
      do i = n1, n2
      a(i) = 0.
      end do
      return 
      end
