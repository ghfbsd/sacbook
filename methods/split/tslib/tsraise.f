c
c   routine TSRAISE takes ss & raises the first ns elements to the n'th 
cpower
c
      subroutine tsraise(ss, ns, n)
c
      real ss(*)
      do i = 1, ns
      ss(i) = ss(i) ** n
c
      end do
      return 
      end
