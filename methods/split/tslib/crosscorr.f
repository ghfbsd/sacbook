c
c   routine CROSSCORR forms the cross-correlation of two time
c   series s1,s2 of length ns, and places the result in ss.
c
      subroutine crosscorr(s1, s2, ns, ss)
      parameter (ndmax = 65536)
      real s1(ns), s2(ns), ss(2*ns)
c
c   pad series to double length
      real st(ndmax), su(ndmax)
      nss = 2 * ns
      if (nss .gt. ndmax) pause '**CROSSCORR - too many points.'

c    copy double length series S1 in reverse order to ST
c   make sure first NS entries of ST are zero
      call tsreverse(ns, s1, st(ns+1))
c
      call tszero(st, ns)
c   copy S2 to first NS entries of SU
      call tszero(su, nss)
c
      call tscopy(s2, su, ns)
c
      call convolvnp(st, su, nss, ss)
      call tsreverse(nss, ss, st)
c
      call tscopy(st, ss, nss)
      return 
      end
