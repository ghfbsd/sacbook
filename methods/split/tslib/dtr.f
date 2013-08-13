      subroutine dtr(a, n, ifl)
      dimension a(1)
      if (ifl .eq. 1) then
      call demean(a, n)
      else if (ifl .eq. 2) then
      call demean(a, n)
      call dtrend(a, n)
      else
      write(unit=*, fmt='('' INVALID VALUE OF IFL IN DTR.PGM ABORTED'')'
     &) 
      call exit
      end if
      return 
      end
