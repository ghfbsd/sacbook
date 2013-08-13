c
c   routine TSWRITE writes ss(i),i=nmin to nmax
c   to an unformatted file open on unit iunit. The time values
c   are also output.
c   
      subroutine tswrite(ss, nmin, nmax, tzero, dt, iunit)
c
      real ss(*), tt(20000)
      if (dt .eq. 0.0) then
      write(unit=6, fmt=*) ' dt is zero in TSWRITE - stop '
      call exit
c
      end if
      t = tzero
      do i = nmin, nmax
      tt(i) = t
      t = t + dt
c
      end do
      nn = (nmax - nmin) + 1
c
      write(unit=iunit) nn, (tt(i), ss(i),i = nmin, nmax)
      return 
      end
