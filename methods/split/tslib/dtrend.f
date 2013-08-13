      subroutine dtrend(a, n)
      dimension a(1)
      xx = 0.
      xa = 0.
      cent = float(n + 1) / 2.
      do 10 i = 1, n
      x = float(i) - cent
      xx = xx + (x ** 2)
   10 xa = xa + (x * a(i))
      slope = xa / xx
      do 20 i = 1, n
      x = float(i) - cent
   20 a(i) = a(i) - (slope * x)
      return 
      end
