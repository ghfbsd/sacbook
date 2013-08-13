c***	linearly interpolates an unequally spaced time series
c	to produce an equally spaced one with spacing dx
c	n is # pts in original series,npts is the number of
c	equally spaced time points, and nmax is the maximum
c	allowable.  if npts>nmax,npts is set to nmax and a warning
c	is issued. Note: npts is set by the subroutine, not by you.
c	Output is in ybuf. Paul Silver 2/13/85
      subroutine linpol(x, y, n, ybuf, dx, nmax, npts)
c	xlen=x(n)-x(1)
c	npts=xlen/dx +1.5
c	write(6,*)xlen,npts
c	if(npts.gt.nmax)then
c	  write(6,*)'# pts too large, pgm aborted'
c	  call exit
c	endif
      dimension x(n), y(n), ybuf(*)
      x1 = x(1)
      x2 = x(2)
      xnow = x1
      y1 = y(1)
      y2 = y(2)
      indx = 2
c***	fill ybuf with linearly interpolated time series
      ifl = 0
      i = 0
      do while (indx .ge. 0)
    1 if (xnow .gt. x2) then
      indx = indx + 1
      if (indx .gt. n) goto 2
      y1 = y2
      x1 = x2
      y2 = y(indx)
      x2 = x(indx)
      ifl = 0
      goto 1
      end if
      if (ifl .eq. 0) then
      slope = (y2 - y1) / (x2 - x1)
      ifl = 1
      end if
      i = i + 1
      if (i .gt. nmax) then
      write(unit=6, fmt=*) 
     &'warning: npts is greater than maximum.only max is given'
      i = i - 1
      goto 2
      end if
      ybuf(i) = y1 + (slope * (xnow - x1))
c	  write(6,*)ybuf(i),xnow
      xnow = xnow + dx
      end do
    2 npts = i
      return 
      end
