      subroutine dcm(a, adec, n, icnt, xlen, idec)
      common /ppt/ ppimax
      dimension a(1), adec(1)
      ppi = float(n) / xlen
      idec = ((ppi - 1) / ppimax) + 1.
      if (idec .eq. 1) then
      do 10 i = 1, n
      adec(i) = a(i)
   10 continue
      icnt = n
      return 
      end if
      id2 = idec * 2
      icnt = 0
      do 20 i = 1, n, id2
      j = min(id2,(n - i) + 1)
      call scdcm(a(i), j, a1, a2)
      icnt = icnt + 1
      adec(icnt) = a1
      icnt = icnt + 1
      adec(icnt) = a2
   20 continue
      return 
      end
