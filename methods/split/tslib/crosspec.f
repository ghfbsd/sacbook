c
c   routine CROSSPEC calculates the cross-spectrum (sum of AB*)for s1,s2
c   (elements n1 to n2), and returns it as CSP (complex).
c
      subroutine crosspec(s1, s2, n1, n2, csp)
c
      complex s1(*), s2(*), csp, z
      csp = cmplx(0.0,0.0)
      do i = n1, n2
      z = s1(i) * conjg(s2(i))
      csp = csp + z
c
      end do
      return 
      end
