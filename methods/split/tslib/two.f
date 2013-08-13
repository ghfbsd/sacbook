c SEE SUBROUTINE NTWO -- DOES NOT CLOBBER LEN (BUG SOURCE FOR USERS OF TWO)
c PAD THE REST OF THE ARRAY WITH ZERO'S
c WHEN LEN IS NOT A POWER OF TWO
      subroutine two(a, len)
      dimension a(1)
      itwo = 1
    1 itwo = itwo * 2
      if (len - itwo) 2, 4, 1
    2 continue
      do 3 i = len + 1, itwo
    3 a(i) = 0.
      len = itwo
    4 continue
      return 
      end
