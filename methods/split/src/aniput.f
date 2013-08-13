	subroutine aniput(r,t,rshift,tshift,npts,dt,tlag,azimuth,baz)
c	routine to put in anisotropy with fast azimuth angle (in degrees
c	clockwise from north) and timeshift deltat. 
c	and reconstruct r and t seismograms.	
c	baz is back azimuth in degrees.
c	the rotation is done so that the same arrays can be used in
c	r,rshift and t, tshift  
	parameter (nmax=50000)
	parameter (pi=3.141592654)
	dimension r(1),t(1),rshift(1),tshift(1),ue(nmax),un(nmax),f(nmax),
     +	s(nmax)
	if(npts.gt.nmax)then
	  write(*,*)'npts gt nmax',npts,nmax,' pgm aborted'
	  call exit
	endif
	rad=180./pi
	bangle=baz/rad
c	note: azimuth is clockwise from north. angle is counterclockwise
c	angle from east, taken to be x coord (n is y)
	angle=(90.-azimuth)/rad	
c	derotate to n,e
	call rotsub(r,t,un,ue,npts,bangle)
c	rotate into fast and slow direction. Remember, azimuth is clockwise
c	from north, baz is clockwise from north.	
	call rotcomp(ue,un,f,s,npts,angle)
c	shift f with respect to s
c	shift seismograms by tlag to add anisotropy and rotate into
c	e-w n-s coord	
     	jshift0=tlag/dt          
	call rotcomp(f(1+jshift0),s,
     +	  ue,un,npts,-angle)
c	then to radial,transverse
	call rotsub(un,ue,rshift,tshift,npts,bangle)
	return
	end
