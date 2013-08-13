c
c   routine SCALAR multiplies the first ns elements of ss by scal
c
      subroutine scalar(ss, ns, scal)
c
      real ss(*)
      do i = 1, ns
      ss(i) = scal * ss(i)
c
      end do
      return 
      end
