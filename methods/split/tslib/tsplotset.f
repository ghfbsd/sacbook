c
c   routine TSPLOTSET sets up a multiple time series plot
c   it returns appropriate values of height, width etc
c   Allows 10 time series per screen
c
      subroutine tsplotset(w, h, hplot, ymin, xmin, ymax)
      w = 10.0
      h = 13.0
      nts = 10
      hplot = w / float(nts)
      ymin = 1.0
      xmin = 1.0
c
      ymax = w
      return 
      end
