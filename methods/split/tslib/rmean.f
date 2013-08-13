      subroutine rmean(c, a, npts, nr)
c IT MIGHT BE EASIER IF NR, THE NUMBER OF FILTER POINTS BE EVEN
      dimension c(1), a(1)
      nh = int(nr / 2)
      nb = (npts - nh) + 1
      do i = 1, npts
c
c FIRST NR POINTS
      a(i) = 0.
      if (i .le. nh) then
      is = i - 1
      do k = i - is, i + is
      a(i) = a(i) + c(k)
      end do
      nds = (2 * is) + 1
c
c FROM NR TO RIGHT BEFORE LAST NR POINTS
      a(i) = a(i) / nds
      else if ((i .gt. nh) .and. (i .lt. nb)) then
      do k = i - nh, i + nh
      a(i) = a(i) + c(k)
      end do
c
c LAST NR POINTS
      a(i) = a(i) / (nr + 1)
      else if ((i .ge. nb) .and. (i .lt. npts)) then
      idf = npts - i
      do k = i - idf, i + idf
      a(i) = a(i) + c(k)
      end do
      ndb = (2 * idf) + 1
c LST POINT
      a(i) = a(i) / ndb
      else if (i .eq. npts) then
      a(i) = a(i - 1)
      end if
      end do
      return 
      end
