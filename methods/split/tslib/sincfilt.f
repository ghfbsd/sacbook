c   routine adds (dt>0.0) or takes out (dt<0.0) a sinc**4 function
c
      subroutine sincfilt(a, len, dt)
      logical ladd
      complex a(*)
c
      pi = 4.0 * atan(1.0)
      ladd = dt .gt. 0.0
      dt0 = abs(dt)
      nt2 = (len / 2) + 1
      fny = 0.5 / dt0
      df = 1.0 / (len * dt0)
c
      fr = df
      do i = 2, nt2
      x = (pi * fr) / fny
      sinc = sin(x) / x
      sinc = sinc ** 4
      if (ladd) then
      a(i) = sinc * a(i)
      else
      a(i) = a(i) / sinc
      end if
      fr = fr + df
c
      end do
      return 
      end
