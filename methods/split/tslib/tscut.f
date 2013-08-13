c
c   routine TSCUT takes time series s, length ns, & extracts wavelet
c   with start point n1, end point n2. This is output as ss, starting
c   at point n3 and padded with zeros to total length ns.
c
      subroutine tscut(s, n1, n2, ss, n3, ns)
c
      real s(*), ss(*)
      do i = 1, ns
      ss(i) = 0.0
c
      end do
      nup = (n2 - n1) + n3
      if (((((((n1 .lt. 0) .or. (n1 .gt. ns)) .or. (n2 .lt. 0)) .or. (n2
     & .gt. ns)) .or. (n1 .gt. n2)) .or. (n3 .lt. 1)) .or. (nup .gt. ns)
     &) then
      write(unit=6, fmt=*) ' error in TSCUT - ns,n1,n2,n3 = ', ns, n1, 
     &n2, n3
      call exit
c
      end if
      do i = n1, n2
      ii = (i - n1) + n3
      ss(ii) = s(i)
c
      end do
      return 
      end
