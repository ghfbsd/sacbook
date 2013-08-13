c
c   routine FTAPER applies a cosine taper to the frequency domain values
c A.
c   The taper decays or grows as ISW is +1 or -1.
c
      subroutine ftaper(a, len, isw)
      complex a(1)
      data pi / 3.141592654 /
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
