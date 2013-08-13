      subroutine zero(a, len1, len)
      complex a(1)
      do 10 i = len1, len
   10 a(i) = (0.,0.)
      return 
      end
