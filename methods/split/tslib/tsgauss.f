c
c   routine TSGAUSS sets up a Gaussuan wavelet of length NPTS,
c   centered on N0 & with standard deviation SD
c
      subroutine tsgauss(ss, npts, n0, sd)
c
      real ss(*)
      if ((n0 .lt. 1) .or. (n0 .gt. npts)) then
      write(unit=6, fmt=*) ' Error in TSGAUSS - n0 = ', n0
      call exit
      end if
      if (sd .eq. 0.0) then
      write(unit=6, fmt=*) ' Error in TSGAUSS - sd = 0.0 '
      call exit
c
      end if
      do i = 1, npts
      arg = float(i - n0)
      arg = - ((arg ** 2) / (2.0 * sd))
      if (arg .gt. (-30.0)) then
      ss(i) = exp(arg)
      else
      ss(i) = 0.0
      end if
c
      end do
      call tsnorm(ss, 1, npts, 2, snorm)
      sfac = 1.0 / snorm
c
      call scalar(ss, npts, sfac)
      return 
      end
