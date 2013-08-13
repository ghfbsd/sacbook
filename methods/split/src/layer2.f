c   program layer2
c   this program computes the apparent delay time and fast polarization
c   direction for two layers.  The inputs are the two delay times and two
c   fast polarization directions and the reference frequency.  The output
c   is the apparent fast polarization, delay time the  apparent rotation
c   and change in size of the initial polarization vector.
c   phi1 and dt1 correspond to the first layer the wave passes though
c    with phi2 and dt2 the second.  The solved for parameters are
c   phi0 and dt0, the apparent fast pol dir and delay time, as well as
c   amp, and phase, the amplitude and phase of the complex scalar multiplying
c    the phat vector. Written by Paul Silver 7/25/91
      complex ak
      data pi/3.141592654/
      rad= 180./pi
      write(stderr,*)'give phi1,dt1,phi2,dt2,freq'
      read(*,*) phi1,dt1,phi2,dt2,freq
      t1=pi*freq*dt1
      t2=pi*freq*dt2
      write(*,100) dt1,phi1,dt2,phi2,freq
      do 10 i = 1,181
        phat=float(i-1)
        al1=2.*(phi1-phat)/rad
        al2=2.*(phi2-phat)/rad
        cc  =cos(t1)*sin(t2)*cos(al2) + cos(t2)*sin(t1)*cos(al1)
        cs  =cos(t1)*sin(t2)*sin(al2) + cos(t2)*sin(t1)*sin(al1)
        ap  =cos(t1)*cos(t2) - sin(t1)*sin(t2)*cos(al2-al1)
        app=                 - sin(t1)*sin(t2)*sin(al2-al1)
        al0=atan((app**2 + cs**2)/(app*ap + cs*cc))
        t0= atan(app/(cs*cos(al0) - cc*sin(al0)))
        ak=cmplx(ap,cc)/cmplx(cos(t0),sin(t0)*cos(al0))
        amp=cabs(ak)
        phase=atan2(aimag(ak),real(ak))
        phi0 = .5*al0*rad+phat
        dt0=t0/(pi*freq)
        phase=phase*rad
c     if dt0 is negative, make pos and add 90 to phi0
        if(dt0.lt.0.)then
          dt0=-dt0
          phi0=phi0+90.
          if(phi0.gt.90.)phi0=phi0-180.
        endif
	write(*,*)'apcs ',ap*cs,'appcc',app*cc
        write(*,110) phat, phi0,dt0,amp, phase
110   format(2(1x,f7.2),f10.2,1x,2f8.3)
10    continue
      stop
100   format(2(1x,f7.2,f10.2),1x, f8.3)
      end


