c
c   routine looks at array , and if all elements are zero it replaces
c   the first by +1 & the second by -1
c
      subroutine zerochk(array, npts)
c
      real array(*)
      if (npts .lt. 2) return 
      do i = 1, npts
      if (array(i) .ne. 0.0) return 
c
      end do
      array(1) = 1.0
c
      array(2) = -1.0
      return 
      end
