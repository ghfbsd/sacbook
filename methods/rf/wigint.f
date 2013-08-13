      function wigint(x,y,npts,dx,eps,t)
*=====================================================================
* PURPOSE:  Interpolates evenly or unevenly spaced data.
*=====================================================================
* INPUT ARGUMENTS:
*    X:       X array if unevenly spaced, first x if evenly spaced. [f]
*    Y:       Y array. [fa]
*    NPTS:    Length of (X and) Y arrays. [i]
*    DX:      Set to 0.0 if unevenly spaced, to sampling interval
*             if evenly spaced. [f]
*    EPS:     Interpolation factor. [f]
*    T:       Time value to interpolate to. [f]
*=====================================================================
* OUTPUT ARGUMENTS:
*    Function result:       Interpolated y value. [f]
*=====================================================================
* REFERENCE: Wiggins, 1976, BSSA, 66, p.2077.
*=====================================================================

      parameter (tiny = 1e-6)
      dimension x(npts),y(npts)
      epsi=.0001
      if(eps.gt.0.) epsi=eps
      if (dx.ne.0.) then
C        Evenly-spaced data - compute index into y array.
	 j=int((t-x(1))/dx)
	 dxj=t-x(1)-float(j)*dx
	 j=j+1
	 if(dxj.eq.0.) go to 99
	 h=dx
	 dxj1=dxj-h
	 dxd=h
	 dxu=h
      else
C        Unevenly spaced data - do binary search for desired value.
	 ilo = 1
	 ihi = npts
         do 20 i=1,npts
	    j = (ilo + ihi)/2
	    if (t .eq. x(j)) go to 99
	    if (t .lt. x(j)) then
	       ihi = j
	    else
	       ilo = j
	    endif
	    if (ihi-ilo .le. 1) go to 30
  20     continue
  30     continue
	 j = ihi-1
	 dxj=t-x(j)
	 if(dxj.eq.0.) go to 99
	 h = max(x(j+1)-x(j), tiny)
	 dxj1=t-x(j+1)
      endif
      hs=h*h
      hc=hs*h
      dxjs=dxj*dxj
      dxj1s=dxj1*dxj1
      dy=y(j+1)-y(j)
      am=dy/h
      amd=am
      amu=am
      if(j.gt.1) then
         if(dx.ne.0.) then
	    dxd = dx
         else
	    dxd = max(x(j)-x(j-1), tiny)
         endif
	 dyd=y(j)-y(j-1)
	 amd=dyd/dxd
      endif
      if(j+1.lt.npts) then
         if(dx.ne.0.) then
	    dxu = dx
         else
	    dxu = max(x(j+2)-x(j+1), tiny)
         endif
	 dyu=y(j+2)-y(j+1)
	 amu=dyu/dxu
      endif
      wd=1./amax1(abs(amd),epsi)
      w=1./amax1(abs(am),epsi)
      wu=1./amax1(abs(amu),epsi)
      sp=(wd*amd+w*am)/(wd+w)
      sp1=(w*am+wu*amu)/(w+wu)
      t1=y(j)*(dxj1s/hs+2.*dxj*dxj1s/hc)
      t2=y(j+1)*(dxjs/hs-2.*dxj1*dxjs/hc)
      t3=sp*dxj*dxj1s/hs
      t4=sp1*dxjs*dxj1/hs
      wigint=t1+t2+t3+t4
      go to 100
  99  wigint=y(j)
 100  return
      end
