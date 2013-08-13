c
c     function rsple(i1,i2,x,y,q,s)
c $$$$$ calls only library routines $$$$$
c
c   rsple returns the value of the function y(x) evaluated at point s
c   using the cubic spline coefficients computed by rspln and saved i
c   q.  if s is outside the interval (x(i1),x(i2)) rsple extrapolates
c   using the first or last interpolation polynomial.  the arrays mus
c   be dimensioned at least - x(i2), y(i2), and q(3,i2).
c
c                                                     -rpb
      function rsple(i1, i2, x, y, q, s)
      dimension x(1), y(1), q(3, 1)
c   guarantee i within bounds.
      data i / 1 /
# 15 "rsple.for"
      ii = i2 - 1
# 17 "rsple.for"
      i = max0(i,i1)
c   see if x is increasing or decreasing.
# 18 "rsple.for"
      i = min0(i,ii)
c   x is decreasing.  change i as necessary.
# 20 "rsple.for"
      if (x(i2) - x(i1)) 1, 2, 2
# 22 "rsple.for"
    1 if (s - x(i)) 3, 3, 4
    4 i = i - 1
      if (i - i1) 11, 6, 1
    3 if (s - x(i + 1)) 5, 6, 6
    5 i = i + 1
c   x is increasing.  change i as necessary.
# 27 "rsple.for"
      if (i - ii) 3, 6, 7
# 29 "rsple.for"
    2 if (s - x(i + 1)) 8, 8, 9
    9 i = i + 1
      if (i - ii) 2, 6, 7
    8 if (s - x(i)) 10, 6, 6
   10 i = i - 1
      if (i - i1) 11, 6, 8
    7 i = ii
      goto 6
c   calculate rsple using spline coefficients in y and q.
# 37 "rsple.for"
   11 i = i1
# 39 "rsple.for"
    6 h = s - x(i)
      rsple = y(i) + (h * (q(1,i) + (h * (q(2,i) + (h * q(3,i))))))
      return 
      end
