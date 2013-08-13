c
c   routine CONVOLV forms the convolution of s1 & s2 using an FFT.
c   s1,s2 are unchanged, & output is ss.
c   Normalization is such that unit delta * unit delta = unit delta
c   Time series is padded to double length before 
c   padding to power of two, to avoid wrap around.
c
      subroutine convolv(s1, s2, ns, ss)
      parameter (ndmax = 65536)
      real s1(*), s2(*)
c
      real ss(ndmax), st1(ndmax), st2(ndmax)
      nss = 2 * ns
      if (nss .gt. ndmax) then
      write(unit=6, fmt=*) ' Error in CONVOLV - maximum dimension', 
     &' exceeded - abandon '
      call exit
c
      end if
      call tszero(st1, nss)
      call tszero(st2, nss)
      call tscopy(s1, st1, ns)
c
      call tscopy(s2, st2, ns)
c
      lens = nss
      call two(st1, lens)
      call fftl(st1, lens, 1, ierr)
c
      call fftl(st2, lens, 1, ierr)
      len2 = lens / 2
      do i = 1, len2 + 1
      ir = (2 * i) - 1
      ii = ir + 1
      sr = st1(ir)
      si = st1(ii)
      tr = st2(ir)
c	write (6,*) 'i,ir,ii,sr,si,tr,ti =',i,ir,ii,sr,si,tr,ti
      ti = st2(ii)
      st1(ir) = (0.5 * ((sr * tr) - (si * ti))) * float(lens)
      st1(ii) = (0.5 * ((sr * ti) + (si * tr))) * float(lens)
c
      end do
      call fftl(st1, lens, 2, ierr)
c
      call tscopy(st1, ss, nss)
      return 
      end
