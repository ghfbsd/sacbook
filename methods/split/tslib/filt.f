      subroutine filt(a, lenf, df, f0, f1, f2, f3, f4, itype)
      complex a(1)
      if (itype .eq. 0) then
      write(unit=*, fmt=*) 'option not active, pgm aborted'
c
c	subroutine filt(a,lenf,df,f0,f1,f2,f3,f4,itype)
c
c
c	************   BOXCAR FILTER   ********************
c
c fmin and fmax are the min and max frequencies within the boxcar
c filter limits.  Data (a(i)) are zeroed outside of those limits.
c 
c
c
c	  CALL ZERO(A,1,INT((FMIN-F0)/DF+1.5))
c	  CALL ZERO(A,INT((FMAX-F0)/DF+1.5),LENF)
c	  WRITE(6,*)'BOX CAR FILTER'
      call exit
c***	COSINE TAPER WITH LOW FREQUENCY TRUNCATION
      else if (itype .eq. 1) then
      if1 = ((f1 - f0) / df) + 1.5
      if2 = ((f2 - f0) / df) + 1.5
      if3 = ((f3 - f0) / df) + 1.5
      if4 = ((f4 - f0) / df) + 1.5
      call ftaper(a(if1), (if2 - if1) + 1, -1)
      call ftaper(a(if3), (if4 - if3) + 1, 1)
      call zero(a, 1, if1 - 1)
c	  WRITE(6,*)'COSINE TAPER'
      call zero(a, if4 + 1, lenf)
c***	TRAPAZOID FILTER
      else if (itype .eq. 2) then
      call trpfil(a, lenf, df, f0, f1, f2, f3, f4)
      write(unit=6, fmt=*) 'TRAPAZOID FILTER'
      else
      write(unit=*, fmt='('' INVALID VALUE OF ITYPE,PGM ABORTED'')') 
      call exit
      end if
      return 
      end
