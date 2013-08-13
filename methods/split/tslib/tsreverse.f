c
c   routine TSREVERSE reverses the real time series ss,
c   and places the output in st
c
      subroutine tsreverse(ns, ss, st)
      real ss(ns),st(ns)
c
      do i = 1, ns
	 st(i) = ss((ns - i) + 1)
      end do
      return 
      end
