 	SUBROUTINE ROTSUB(AN,AE,AR,AT,NPTS,BAZ)
c+
c	SUBROUTINE ROTSUB(AN,AE,AR,AT,NPTS,BAZ)
c
c	AN and AE are northsouth and eastwest time series respectively;
c    	NPTS = number of points; AR and AT are radial and transverse
c  	time series returned from this routine; BAZ is backazimuth
c
c	Subroutine to rotate into radial and transverse components.
c	BAZ is the back-azimuth, that is, the clockwise angle 
c	measured from north to the earthquake with the station as the
c	vertex. The positive radial direction points in the direction of
c	the earthquake from the station. The positive transverse direction
c	is to the left of the radial direction.  This assures a right-handed
c	coordinate system with the z component pointing up. 
c	Note that this convention is not the same as AKI and RICHARDS
c	page 114-115.  They use a coordinate system where Z is down.
c	At some point it may be advantageous to switch the signs of
c	both the radial and transverse components so that they will
c	have the same convention as the program WKBJ.  In that program,
c	x points to the station from the earthquake, y is in transverse
c	direction, with positive y to the left of x.
c
c	Has been modified by pgs 6/29/89 so that the same two arrays
c	can be used for input and output
c
c-
	DIMENSION AN(NPTS),AE(NPTS),AR(NPTS),AT(NPTS)
	SI=SIN(BAZ)
	CO=COS(BAZ)
	DO 10 I=1,NPTS
	AR_temp=AE(I)*SI+AN(I)*CO 
c	  !pgs
	AT_temp=-(AE(I)*CO-AN(I)*SI)
c	  !pgs
	AR(I)=AR_temp 
c	  !pgs
10	AT(I)=AT_temp 
c	  !pgs
	RETURN
	END
