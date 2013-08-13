c
c     subroutine rspln(i1,i2,x,y,q,f)
c $$$$$ calls only library routines $$$$$
c
c   subroutine rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  the
c   interpolation is continuous with continuous first and second
c   derivitives.  it agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  the arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -rpb
c
      subroutine rspln(i1, i2, x, y, q, f)
      dimension x(1), y(1), q(3, 1), f(3, 1), yy(3)
      equivalence (y0, yy(1))
      data yy / 3*0. /
      data tol / 1.e-13 /
# 21 "rspln.for"
      j1 = i1 + 1
c   bail out if there are less than two points total.
# 22 "rspln.for"
      y0 = 0.
# 24 "rspln.for"
      if (i2 - i1) 13, 17, 8
c   search for discontinuities.
# 25 "rspln.for"
    8 a0 = x(j1 - 1)
# 27 "rspln.for"
      do 3 i = j1, i2
      b0 = a0
      a0 = x(i)
      if (abs(a0 - b0) - tol) 4, 4, 3
    3 continue
   17 j1 = j1 - 1
      j2 = i2 - 2
      goto 5
    4 j1 = j1 - 1
c   see if there are enough points to interpolate (at least three).
# 36 "rspln.for"
      j2 = i - 3
c   only two points.  use linear interpolation.
# 38 "rspln.for"
    5 if ((j2 + 1) - j1) 9, 10, 11
# 40 "rspln.for"
   10 j2 = j2 + 2
      y0 = (y(j2) - y(j1)) / (x(j2) - x(j1))
      do 15 j = 1, 3
      q(j,j1) = yy(j)
   15 q(j,j2) = yy(j)
c   more than two points.  do spline interpolation.
# 45 "rspln.for"
      goto 12
# 47 "rspln.for"
   11 a0 = 0.
      h = x(j1 + 1) - x(j1)
      h2 = x(j1 + 2) - x(j1)
      y0 = (h * h2) * (h2 - h)
      h = h * h
c   calculate derivitive at near end.
# 52 "rspln.for"
      h2 = h2 * h2
# 54 "rspln.for"
      b0 = (((y(j1) * (h - h2)) + (y(j1 + 1) * h2)) - (y(j1 + 2) * h))
     & / y0
c   explicitly reduce banded matrix to an upper banded matrix.
# 55 "rspln.for"
      b1 = b0
# 57 "rspln.for"
      do 1 i = j1, j2
      h = x(i + 1) - x(i)
      y0 = y(i + 1) - y(i)
      h2 = h * h
      ha = h - a0
      h2a = h - (2. * a0)
      h3a = (2. * h) - (3. * a0)
      h2b = h2 * b0
      q(1,i) = h2 / ha
      q(2,i) = - (ha / (h2a * h2))
      q(3,i) = - ((h * h2a) / h3a)
      f(1,i) = (y0 - (h * b0)) / (h * ha)
      f(2,i) = (h2b - (y0 * ((2. * h) - a0))) / ((h * h2) * h2a)
      f(3,i) = - ((h2b - ((3. * y0) * ha)) / (h * h3a))
      a0 = q(3,i)
c   take care of last two rows.
# 72 "rspln.for"
    1 b0 = f(3,i)
# 74 "rspln.for"
      i = j2 + 1
      h = x(i + 1) - x(i)
      y0 = y(i + 1) - y(i)
      h2 = h * h
      ha = h - a0
      h2a = h * ha
      h2b = (h2 * b0) - (y0 * ((2. * h) - a0))
      q(1,i) = h2 / ha
      f(1,i) = (y0 - (h * b0)) / h2a
      ha = x(j2) - x(i + 1)
      y0 = - ((h * ha) * (ha + h))
c   calculate derivitive at far end.
# 85 "rspln.for"
      ha = ha * ha
# 87 "rspln.for"
      y0 = (((y(i + 1) * (h2 - ha)) + (y(i) * ha)) - (y(j2) * h2)) / y0
      q(3,i) = ((y0 * h2a) + h2b) / ((h * h2) * (h - (2. * a0)))
c   solve upper banded matrix by reverse iteration.
# 89 "rspln.for"
      q(2,i) = f(1,i) - (q(1,i) * q(3,i))
# 91 "rspln.for"
      do 2 j = j1, j2
      k = i - 1
      q(1,i) = f(3,k) - (q(3,k) * q(2,i))
      q(3,k) = f(2,k) - (q(2,k) * q(1,i))
      q(2,k) = f(1,k) - (q(1,k) * q(3,k))
    2 i = k
c   fill in the last point with a linear extrapolation.
# 97 "rspln.for"
      q(1,i) = b1
# 99 "rspln.for"
    9 j2 = j2 + 2
      do 14 j = 1, 3
c   see if this discontinuity is the last.
# 101 "rspln.for"
   14 q(j,j2) = yy(j)
c   no.  go back for more.
# 103 "rspln.for"
   12 if (j2 - i2) 6, 13, 13
# 105 "rspln.for"
    6 j1 = j2 + 2
c   there is only one point left after the latest discontinuity.
# 106 "rspln.for"
      if (j1 - i2) 8, 8, 7
# 108 "rspln.for"
    7 do 16 j = 1, 3
c   fini.
# 109 "rspln.for"
   16 q(j,i2) = yy(j)
# 111 "rspln.for"
   13 return 
      end
