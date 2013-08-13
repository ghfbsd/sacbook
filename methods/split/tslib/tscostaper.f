c
c   routine TSCOSTAPER applies a cosine taper of length LEN to time seri
ces A.
c   Taper decays or grows as ISW is 1 or -1. To apply other than at star
ct of
c   time series A, pass to this routine as e.g. A(ipoint).
c   Adapted from FTAPER
c
      subroutine tscostaper(a, len, isw)
      real a(*)
c
      pi = 4.0 * atan(1.0)
      if (len .eq. 1) return 
      arg = pi / float(len - 1)
      if (isw .eq. 1) then
      do 10 i = 1, len
      a(i) = (.5 * (1. + cos(arg * float(i - 1)))) * a(i)
   10 continue
      else if (isw .eq. (-1)) then
      do 20 i = 1, len
      a(i) = (.5 * (1. - cos(arg * float(i - 1)))) * a(i)
   20 continue
      else
      write(unit=*, fmt='('' INCORRECT ISW IN FTAPER:PGM ABORTED'')') 
      call exit
      end if
      return 
      end
