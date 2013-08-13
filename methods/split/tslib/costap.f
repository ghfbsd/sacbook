c
c   APPLIES COSINE TAPER (10% OR 50%)
c   IWINDO=3  50% (CORRESPONDS TO HANNING WINDOW)
c   IWINDO=4  5%
c
      subroutine costap(x, npts, iwindo)
      dimension x(1)
      pi = 3.141592654
      percnt = float(iwindo - 3) * 5.
      if (percnt .le. 1.0) percnt = 50.0
      icos1 = npts * (percnt / 100.0)
      icos2 = (npts - icos1) + 1
      ia = 1
      ib = icos1
      iend = 1
      j = 1
      f = pi / float(icos1 - 1)
  110 do i = ia, ib
      fi = float(i - 1)
      if (iend .eq. 1) then
      x(j) = (x(j) * (1. - cos(fi * f))) * 0.5
      j = j + 1
      else if (iend .eq. 2) then
c		X(J)=X(J)*(1.+COS(FI*F))*0.5
      x(j) = (x(j) * (1. - cos(fi * f))) * 0.5
      j = j - 1
c	   X(J)=X(J)*(1.-COS(FI*F))*0.5 + SIGN(0.5,X(J))
      end if
      end do
      if (iend .eq. 2) return 
      iend = 2
      j = npts
      goto 110
      end
