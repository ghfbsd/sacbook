      subroutine cosfilt(a, lenf, df, f1, f2, f3, f4)
c***	cosine taper with low frequency truncation
      complex a(1)
      if1 = (f1 / df) + 1.5
      if2 = (f2 / df) + 1.5
      if3 = (f3 / df) + 1.5
      if4 = (f4 / df) + 1.5
      call ftaper(a(if1), (if2 - if1) + 1, -1)
      call ftaper(a(if3), (if4 - if3) + 1, 1)
      call zero(a, 1, if1 - 1)
c	  write(6,*)'cosine taper'
      call zero(a, if4 + 1, lenf)
      return 
      end
