C gcd, dazll  --  Subroutines to compute azimuth and distance between two points
C                 on a sphere and to compute the lat, lon of a point at a given
C                 distance from another along a given azimuth
C gcdd -- Subroutine to compute derivative of great circle distance wrt lat &
C         longitude.

	subroutine gcd(elat,elon,slat,slon,del,delkm,az,baz)
C       gcd -- Return great circle distance between two
C          points.
C
C       Called with:
C          elat, elon - Lat & lon of first point
C          slat, slon - Lat & lon of second point
C       Returns:
C          del - distance (in degrees) separating points along minor arc
C          delkm - great circle distance (km)
C          az - azimuth from point 1 to point 2 (degrees)
C          baz - azimuth from point 2 to point 1 (degrees)
C
C       Taken from code written by Bob Hermann, St. Louis U.
C          hacked by G. Helffrich
C

	real*4 rad, e2, re
	data rad/0.017453292/,e2/0.9932315/, re/6371.003/
	call dircos(elat,elon,selat,celat,selon,celon,dse,dce,ea,eb,ec)
	call dircos(slat,slon,sslat,cslat,sslon,cslon,dss,dcs,sa,sb,sc)

c-----
c	compute distance
c	Choose correct formula for short and large distances
c-----
	cdel = (ea*sa + eb*sb + ec*sc)
c-----
c	if DEL = [0,20)
c-----
	if(cdel .gt. 0.9396)then
		fac = (ea-sa)**2 + (eb-sb)**2 + (ec-sc)**2
		fac = sqrt(fac)/2.0
		del = 2.0*asin(fac)
c-----
c	if DEL = [20,160]
c-----
	else if(cdel. le. 0.9396 .and. cdel. ge. -0.9396)then
		del = acos(cdel)
c-----
c	if DEL = (160,180]
c-----
	else
		fac = (ea+sa)**2 + (eb+sb)**2 + (ec+sc)**2
		fac = sqrt(fac)/2.0
		del = 2.0*acos(fac)
	endif
	
c-----
c	check for station or epicenter at pole
c-----
	if(elat.eq.90.0 .or. slat.eq.-90.0)then
		az = 180.0
		baz = 0.0
		saz =  0.0
		caz = -1.0
	else if(elat.eq.-90.0 .or. slat.eq.90.0)then
		az = 0.0
		baz = 180.0
		saz = 0.0
		caz = 1.0
	else
		saz = celat*(cslat * sin(rad*(slon - elon)))
		caz = (sslat - cdel*selat)
c       	fac = sqrt(saz**2 + caz**2)
c       	saz = saz / fac
c       	caz = caz / fac
		az = atan2(saz,caz)
	
		sbz = - cslat*(celat * sin(rad*(slon - elon)))
		cbz = (selat - cdel*sslat)
		baz = atan2(sbz,cbz)
	
 		az = az / rad
		baz = baz / rad
	endif
	delkm = del * re
	del = del / rad

c-----
c	put az and baz in the range [0,360)
c-----
	if(az .lt. 0.0)az = az + 360.0
	if(baz .lt. 0.0)baz = baz + 360.0
	return
	end
	
	subroutine dircos(lat,lon,
     +     slat,clat,slon,clon,dslat,dclat,aa,bb,cc)
	real*4 lat,lon,slat,slon,clat,clon, dclat, dslat, aa, bb, cc
	real*4 rad, e2
	parameter (e2geo=0.993305615, e2equ=0.99776354)
	data rad/0.017453292/,e2/e2geo/
c-----
c	convert geographic latitude to geocentric
c	Use flattening of Chovitz (1981) f= 1/298.257 adopted by IUGG 1980
c	
c	The relation between geocentric and geographic latitude is
c	tan phi c = ( 1 - f)^2 tan phi g
c
c	To avoid problems at poles, define sin phi c and cos phi c
c	so that the relation holds ans also that s^2 + c^2 = 1
c
c	For geographic to geocentric use e2 = (1-f)^2 = 0.993395615
c	For geographic to equidistant use e2 =(1-f)^1.5=0.99776354
c	Brown, R. J. (1984). On the determination of source-receiver
c		distances using a new equidistant latitude,
c		Geophys. J. R. astr. Soc. 76, 445-459.
c
c       Returns derivative of slat and clat wrt lat as well.
c	
c-----
	c = cos(rad*lat)
	s = sin(rad*lat)
	e4 = e2**2
	arg = e4 + (1.0-e4)*c*c
	fac = sqrt(arg)
	slat = e2 * s /fac
	clat =      c /fac
	slon = sin(rad*lon)
	clon = cos(rad*lon)

	dclat = s/fac*(c*c*(1.0-e4)/arg - 1.0)
	dslat = e2*c/fac*(s*s*(1.0-e4)/arg + 1.0)

	aa = clat * clon
	bb = clat * slon
	cc = slat
	
	return
	end

	subroutine gcdd(elat,elon,slat,slon,del,ddella,ddello)
C       gcdd -- Return derivative of great circle distance between two
C          points.
C
C       Called with:
C          elat, elon - Lat & lon of first point
C          slat, slon - Lat & lon of second point
C       Returns:
C          del - distance (in degrees) separating points along minor arc
C          ddella - derivative of delta wrt latitude (unitless)
C          ddello - derivative of delta wrt longitude (unitless)
C
C       G. Helffrich/Bristol 4 Nov. 1993
C

	real*4 rad, e2, re
	data rad/0.017453292/,e2/0.9932315/, re/6371.003/
	call dircos(elat,elon,selat,celat,selon,celon,dse,dce,ea,eb,ec)
	call dircos(slat,slon,sslat,cslat,sslon,cslon,dss,dcs,sa,sb,sc)

c-----  compute derivatives of ea, eb, ec wrt lat & lon
        deala = celon * dce
        dealo =-celat * selon
        debla = selon * dce
        deblo = celat * celon
        decla = dse
        declo = 0.0

c-----
c	compute distance
c	Choose correct formula for short and large distances
c-----
	cdel = (ea*sa + eb*sb + ec*sc)
c-----
c	if DEL = [0,20)
c-----
	if(cdel .gt. 0.9396)then
		arg = (ea-sa)**2 + (eb-sb)**2 + (ec-sc)**2
		fac = sqrt(arg)/2.0
		del = 2.0*asin(fac)

		dargla = 2.0*(
     +             (ea-sa)*deala + (eb-sb)*debla + (ec-sc)*decla
     +          )
		darglo = 2.0*(
     +             (ea-sa)*dealo + (eb-sb)*deblo + (ec-sc)*declo
     +          )
		dfacla = dargla/(8.*fac)
		dfaclo = darglo/(8.*fac)
		ddella = 4.0/sqrt(4.0 - arg)*dfacla
		ddello = 4.0/sqrt(4.0 - arg)*dfaclo
c-----
c	if DEL = [20,160]
c-----
	else if(cdel. le. 0.9396 .and. cdel. ge. -0.9396)then
		del = acos(cdel)

		dcdla = sa*deala + sb*debla + sc*decla
		dcdlo = sa*dealo + sb*deblo + sc*declo
		ddella = -dcdla/sqrt(1.0 - cdel**2)
		ddello = -dcdlo/sqrt(1.0 - cdel**2)
c-----
c	if DEL = (160,180]
c-----
	else
		arg = (ea+sa)**2 + (eb+sb)**2 + (ec+sc)**2
		fac = sqrt(arg)/2.0
		del = 2.0*acos(fac)

		dargla = 2.0*(
     +             (ea+sa)*deala + (eb+sb)*debla + (ec+sc)*decla
     +          )
		darglo = 2.0*(
     +             (ea+sa)*dealo + (eb+sb)*deblo + (ec+sc)*declo
     +          )
		dfacla = dargla/(8.*fac)
		dfaclo = darglo/(8.*fac)
		ddella = -4.0/sqrt(4.0 - arg)*dfacla
		ddello = -4.0/sqrt(4.0 - arg)*dfaclo
	endif
	del = del / rad
	end

      SUBROUTINE GCDSE(SLAT,SLON,RLAT,RLON,DELTA,RANGE,SRAZ,RSAZ)
C     GCDSE -- Compute great circle distance in a spherical earth
C
C     Assumes:
C        SLAT, SLON, RLAT, RLON - Receiver lat & lon (degrees)
C
C     Returns:
C        DELTA - distance (degrees)
C        RANGE - distance (km)
C        SRAZ - source->receiver azimuth
C        RSAZ - receiver->source azimuth

      DEGRAD = ATAN(1.0)/45.0
      THA = DEGRAD*(90.0-SLAT)
      THB = DEGRAD*(90.0-RLAT)
      COSDPH = COS(DEGRAD*(SLON-RLON))
      SINDPH = SIN(DEGRAD*(SLON-RLON))
      COSDEL = MAX(-1.0,MIN(1.0,
     +   COS(THA)*COS(THB)+SIN(THA)*SIN(THB)*COSDPH
     +))
      DELTA = ACOS(COSDEL)/DEGRAD
      RANGE = DELTA * 111.1195
      SINCOSZ = SIN(THA)*COS(THB)-SIN(THB)*COS(THA)*COSDPH
      SINSINZ = -SIN(THB)*SINDPH
      SRAZ = ATAN2(SINSINZ,SINCOSZ)/DEGRAD
      SINCOSZ = SIN(THB)*COS(THA)-SIN(THA)*COS(THB)*COSDPH
      SINSINZ = SIN(THA)*SINDPH
      RSAZ = ATAN2(SINSINZ,SINCOSZ)/DEGRAD
      IF (SRAZ .LT. 0.0) SRAZ = SRAZ + 360.0
      IF (RSAZ .LT. 0.0) RSAZ = RSAZ + 360.0
      END

      subroutine dazell(lat0,lon0,dist,az,lat1,lon1)
C     dazell  --  Compute coordinates of a point at a given distance and
C                  azimuth from another point on an elliptical earth.
C
C     Called via:
C        call dazll(lat0,lon0,dist,az,lat1,lon1)
C
C     Assumes:
C        lat0, lon0 - lat, lon (in degrees) of reference point
C        dist - distance (in degrees) of desired coordinates
C        az - azimuth (clockwise from north in degrees) of direction to
C           desired point
C
C     Returns:
C        lat1, lon1 - lat, lon (in degrees) of desired point in elliptical
C           earth.

      real lat0,lon0,lat1,lon1

C        Use Newton-Raphson to locate exact elliptical earth distance to yield
C           specified distance in kilometers.  Only two iterations are usually
C           required to find distance within 0.1%.
      parameter (tol=1e-3)

      if (abs(dist).lt.tol) then
	 lat1 = lat0
	 lon1 = lon0
	 delta = 0.0
	 delkm = 0.0
	 saz = 0.0
	 baz = 0.0
	 go to 10
      endif
      gcarc = dist
      do i=1,5
	 gclo=0.95*gcarc
	 gchi=1.05*gcarc
	 call dazll(lat0,lon0,gclo,az,lat1,lon1)
	 call gcd(lat0,lon0,lat1,lon1,dello,delkm,saz,baz)
	 dello = sign(dello,dist)
	 call dazll(lat0,lon0,gchi,az,lat1,lon1)
	 call gcd(lat0,lon0,lat1,lon1,delhi,delkm,saz,baz)
	 delhi = sign(delhi,dist)
	 call dazll(lat0,lon0,gcarc,az,lat1,lon1)
	 call gcd(lat0,lon0,lat1,lon1,delta,delkm,saz,baz)
	 delta = sign(delta,dist)
	 ddddel = (delhi-dello)/(gchi-gclo)
	 gcnew = gcarc - (delta-dist)/ddddel
	 if (abs(gcnew-gcarc) .lt. tol) go to 10
	 gcarc = gcnew
      enddo
      write(0,*) '**DAZELL:  Strange - could not find proper distance.'
10    continue
      end

        subroutine dazll(lat0,lon0,dist,az,lat1,lon1)
C       dazll  --  Compute coordinates of a point at a given distance and
C                  azimuth from another point.
C
C       Called via:
C          call dazll(lat0,lon0,dist,az,lat1,lon1)
C
C       Assumes:
C          lat0, lon0 - lat, lon (in degrees) of reference point
C          dist - distance (in degrees) of desired coordinates
C          az - azimuth (clockwise from north in degrees) of direction to
C             desired point
C
C       Returns:
C          lat1, lon1 - lat, lon (in degrees) of desired point.
C
C       See formulas in Appendix A, "Seismology and Plate Tectonics",
C       D. Gubbins, Cambridge, 1990

        parameter (pi=3.1415926535897, degrad=pi/180.)
	real lat0, lon0, lat1, lon1
        
	rlat0 = (90.0 - lat0)*degrad
	raz = az*degrad
	rdist = dist*degrad
	cosla1 = sin(rlat0)*cos(raz)*sin(rdist) + cos(rdist)*cos(rlat0)

	lat1 = 90.0 - acos(cosla1)/degrad
	lon1 = lon0 + atan2( sin(raz)*sin(rdist)*sin(rlat0),
     +                       cos(rdist) - cos(rlat0)*cosla1)
     +     / degrad
	end
