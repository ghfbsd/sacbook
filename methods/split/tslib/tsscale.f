c
c   routine TSSCALE scales first n elements of
c   array ss such that Max(ss(i)))=1.0 
c   sfac is factor that array was multiplied by (output)
c   **	calls SCALAR	**
c
      subroutine tsscale(ns, ss, sfac)
      real ss(*)
      amax = 0.0
      do i = 1, ns
      as = abs(ss(i))
      if (as .gt. amax) amax = as
c
      end do
      if (amax .eq. 0.0) return 
      as = 1.0 / amax
      call scalar(ss, ns, as)
c
      sfac = as
      return 
      end
