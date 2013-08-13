c
c   routine TSDEMEAN takes a time series S1 of NP points, 
c   & removes the mean
c
      subroutine tsdemean(np, s1)
c
      real s1(*)
      smean = 0.0
      do i = 1, np
      smean = smean + s1(i)
      end do
c
      smean = smean / float(np)
      do i = 1, np
      s1(i) = s1(i) - smean
c
      end do
      return 
      end
