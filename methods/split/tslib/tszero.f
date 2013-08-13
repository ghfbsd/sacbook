c
c   routine TSZERO sets the first ns elements of ss to zero
c
      subroutine tszero(ss, ns)
c
      real ss(*)
      do i = 1, ns
      ss(i) = 0.0
c
      end do
      return 
      end
