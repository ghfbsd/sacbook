      subroutine tscopy(s1, s2, ns)
c
      real s1(*), s2(*)
      do i = 1, ns
      s2(i) = s1(i)
c
      end do
      return 
      end
