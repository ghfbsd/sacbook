c
c $$$$$ CALLS FFT AND REALTR $$$$$
c
c   IF IABS(NDIR).EQ.1 FFTL FOURIER TRANSFORMS THE N POINT REAL TIME
c   SERIES IN ARRAY X.  THE RESULT OVERWRITES X STORED AS
c   (N+2)/2 COMPLEX NUMBERS (NON-NEGATIVE FREQUENCIES ONLY).  IF
c   IABS(NDIR).EQ.2 FFTL FOURIER TRANSFORMS THE (N+2)/2 COMPLEX FOURIER
c   COEFFICIENTS (NON-NEGATIVE FREQUENCIES ONLY) IN ARRAY X
c   (ASSUMING THE SERIES IS HERMITIAN).  THE RESULTING N POINT REAL
c   TIME SERIES OVERWRITES X.  IF NDIR.GT.0 THE FOREWARD TRANSFORM USES
c   THE SIGN CONVENTION EXP(I*W*T).  IF NDIR.LT.0 THE FOREWARD TRANSFORM
c   USES THE SIGN CONVENTION EXP(-I*W*T).  THE FOREWARD TRANSFORM IS
c   NORMALIZED SUCH THAT A SINE WAVE OF UNIT AMPLITUDE IS TRANSFORMED
c   INTO DELTA FUNCTIONS OF UNIT AMPLITUDE.  THE BACKWARDS TRANSFORM IS
c   NORMALIZED SUCH THAT TRANSFORMING FOREWARD AND THEN BACK RECOVERS
c   THE ORIGINAL SERIES.  IERR IS NORMALLY ZERO.  IF IERR.EQ.1 THEN FFT
c   HAS NOT BEEN ABLE TO FACTOR THE SERIES.  HOWEVER, X HAS BEEN
c   SCRAMBLED BY REALTR.  NOTE THAT IF N IS ODD THE LAST POINT WILL NOT
c   BE USED IN THE TRANSFORM.
c
c                                                     -RPB
      subroutine fftl(x, n, ndir, ierr)
      integer*4 n, n2, n1, i
      dimension x(*)
      n2 = n / 2
      idir = iabs(ndir)
c   DO FOREWARD TRANSFORM (IE. TIME TO FREQUENCY).
      goto (1, 2), idir
    1 call fft(x, x(2), n2, n2, n2, 2, ierr)
      call realtr(x, x(2), n2, 2)
      n1 = (2 * n2) + 2
      scale = 1. / n
      if (ndir .gt. 0) goto 3
      do 5 i = 4, n, 2
    5 x(i) = - x(i)
c   DO BACKWARD TRANSFORM (IE. FREQUENCY TO TIME).
      goto 3
    2 if (ndir .gt. 0) goto 6
      do 7 i = 4, n, 2
    7 x(i) = - x(i)
    6 x(2) = 0.
      x((2 * n2) + 2) = 0.
      call realtr(x, x(2), n2, -2)
      call fft(x, x(2), n2, n2, n2, -2, ierr)
      n1 = 2 * n2
c   NORMALIZE THE TRANSFORM.
      scale = .5
    3 do 4 i = 1, n1
    4 x(i) = scale * x(i)
      return 
      end
