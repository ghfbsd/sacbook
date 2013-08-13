C     Common declaration for our version of the SAC file header.

      parameter (lencom = 158, lenchr=24)
      character sachdr(lencom)*4, khdr(lenchr)*8, kcom(lenchr)*8
      common /cmhdr/ sachdr
      common /kmhdr/ khdr
      equivalence (nvhdr,sachdr(77)),(sachdr(111),kcom(1))
