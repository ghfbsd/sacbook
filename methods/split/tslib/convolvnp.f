c
c   routine CONVOLVNP forms the convolution of s1 & s2 using an FFT.
c   s1,s2 are unchanged, & output is ss.
c   all arrays should be dimensioned ns<65536
c   Normalization is such that unit delta * unit delta = unit delta
c
      subroutine convolvnp(s1, s2, ns, ss)
      parameter (ndmax = 65538)
      real s1(*), s2(*), ss(*)
c
      real st(ndmax)
      if ((ns + 2) .gt. ndmax) then
      write(unit=6, fmt=*) ' Error 1 in CONVOLVNP - maximum dimension '
     &, 'exceeded - abandon '
      call exit
c
      end if
      do i = 1, ns
      ss(i) = s1(i)
      st(i) = s2(i)
      end do
c
      lens = ns
      call two(ss, lens)
      if ((lens + 2) .gt. ndmax) then
      write(unit=6, fmt=*) ' Error 2 in CONVOLVNP - maximum dimension '
     &, 'exceeded - abandon '
      call exit
c
      end if
      call fftl(ss, lens, 1, ierr)
c
      call fftl(st, lens, 1, ierr)
      len2 = lens / 2
      do i = 1, len2 + 1
      ir = (2 * i) - 1
      ii = ir + 1
      sr = ss(ir)
      si = ss(ii)
      tr = st(ir)
c	write (6,*) 'i,ir,ii,sr,si,tr,ti =',i,ir,ii,sr,si,tr,ti
      ti = st(ii)
      ss(ir) = (0.5 * ((sr * tr) - (si * ti))) * float(lens)
      ss(ii) = (0.5 * ((sr * ti) + (si * tr))) * float(lens)
c
      end do
c
      call fftl(ss, lens, 2, ierr)
      return 
c
      end
