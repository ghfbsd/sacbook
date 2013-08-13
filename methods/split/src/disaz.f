c
c in:
c	epip 	epicenter latitude
c	epil 	epicenter longitude
c	stnp 	station latitude
c	stnl 	station longitude
c out:
c	az 	azimuth in degree
c	baz	back azimuth in degree
c	range	geocentric distance in km
c	angi	geocentric distance in degree
c
c     north and east are positive
      subroutine disaz(epip, epil, stnp, stnl, az, baz, range, angi)
# 15 "disaz.for"
      epil = - epil
      stnl = - stnl
      pi = 3.14159265358
      esq1 = 0.9932315
      rad = 0.01745329
      epipr = epip * rad
      epipr = atan((esq1 * sin(epipr)) / cos(epipr))
      epilc = 360.0 - epil
      epilr = epilc * rad
      epipc = 90.0 - (epipr * 57.295780)
      epipr = epipc * rad
      num = 1
      do 999 k = 1, num
      stpc = stnp * rad
      stpc = 90. - (57.295780 * atan(esq1 * (sin(stpc) / cos(stpc))))
      stlc = 360. - stnl
      if (stpc - 180.) 69, 68, 68
   68 angi = 180. - epipc
      az = 180.
      baz = 0.
      goto 998
   69 stnpr = stpc * rad
      stnlr = stlc * rad
      pang = abs(stlc - epilc)
      if (pang - 180.) 71, 71, 70
   70 pang = 360. - pang
   71 pang = pang * rad
      isig = 0
      isig = 1
      angi = (cos(epipr) * cos(stnpr)) + ((sin(epipr) * sin(stnpr)) * 
     &cos(pang))
# 45 "disaz.for"
      if (abs(angi) - 1.0001) 600, 200, 300
  600 if (abs(angi) - 1.0) 201, 201, 200
  200 angi = 1.0
  201 an1 = angi
      asang = asin(an1)
      angi = (pi / 2.) - asang
      isig = 2
      az = (cos(stnpr) - (cos(epipr) * cos(angi))) / (sin(epipr) * sin(
     &angi))
# 53 "disaz.for"
      if (abs(az) - 1.0001) 610, 210, 300
  610 if (abs(az) - 1.0) 211, 211, 210
  210 az = 1.0
  211 an2 = az
      asaz = asin(an2)
      az = 57.29578 * ((pi / 2.) - asaz)
      isig = 3
      baz = (cos(epipr) - (cos(stnpr) * cos(angi))) / (sin(stnpr) * sin(
     &angi))
# 61 "disaz.for"
      if (abs(baz) - 1.0001) 620, 220, 300
  620 if (abs(baz) - 1.0) 221, 221, 220
  220 baz = 1.0
  221 an3 = baz
      asbaz = asin(an3)
      baz = 57.29578 * ((pi / 2.) - asbaz)
      goto 303
  300 a = angi
      b = az
      c = baz
      goto (201, 211, 221), isig
  303 angi = angi * 57.295780
      if (stlc - epilc) 72, 74, 73
   72 if ((stlc - epilc) + 180.) 80, 77, 81
   73 if ((stlc - epilc) - 180.) 80, 77, 81
   74 if (stpc - epipc) 75, 78, 76
   75 az = 0.
      baz = 180.
      goto 998
   76 az = 180.
      baz = 0.
      goto 998
   77 if ((stpc + epipc) - 180.) 78, 78, 79
   78 az = 0.
      baz = 0.
      goto 998
   79 az = 180.
      baz = 180.
      goto 998
   80 baz = 360. - baz
      goto 998
   81 az = 360. - az
  998 range = angi * 111.195
  999 continue
      epil = - epil
      stnl = - stnl
      return 
      end
