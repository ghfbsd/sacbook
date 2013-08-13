      subroutine hamwin(x, n)
      dimension x(1)
      tupi = 2. * 3.141592654
      win = tupi / float(n - 1)
      do i = 1, n
      w1 = i - 1
      sg = sign(0.5,x(i))
      x(i) = (x(i) * (0.08 + (0.46 * (1. - cos(w1 * win))))) + sg
      end do
      return 
      end
