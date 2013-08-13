C     ntwo -- Pad an array with zeroes when it isn't a power of 2.
C
C     Assumes:
C        a - array to be padded
C        len - number of elements in a
C
C     Returns:
C        a - optionally padded with zeroes
C        function result - padded length (power of two)
      function ntwo(a, len)
      dimension a(*)
      itwo = 1
    1 itwo = itwo * 2
      if (len - itwo) 2, 4, 1
    2 continue
      do 3 i = len + 1, itwo
    3 a(i) = 0.
    4 continue
      ntwo = itwo
      return 
      end
