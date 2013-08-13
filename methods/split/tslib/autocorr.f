c Changed by G. Bokelmann, Dec. 1989
c
c   routine AUTOCORR calculates the autocorrelation for ss,
c   using an FFT. Output is auto, and ss is unchanged.
c   Normalized so that autocorr(unit delta) = unit delta
c
      subroutine autocorr(ss, ns, auto)
      parameter (ndmax = 65536)
      real ss(*)
c
      real auto(ndmax), st(ndmax)
      nss = 2 * ns
      if (nss .gt. ndmax) then
      write(unit=6, fmt=*) ' Error in AUTOCORR - maximum dimension ', 
     &'exceeded - abandon '
      call exit
c
      end if
      call tsreverse(nss, ss, st)
      call tszero(st, ns)
c
      call convolvnp(ss, st, nss, auto)
      call tsreverse(nss, auto, st)
c
      call tscopy(st, auto, nss)
      return 
      end
