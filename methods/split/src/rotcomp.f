	subroutine rotcomp(x,y,xrot,yrot,npts,angle)
c	rotate seismograms into coord system xrot,yrot from x, y. Angle
c	assumed to be in radians.  Angle is counter-clockwise angle from
c	the positive x axis in radians
c	the rotation is done so that the same arrays can be used in
c	x,xrot and y,yrot. ei x,y overwritten with xrot,yrot
	dimension x(npts),y(npts),xrot(npts),yrot(npts)
	si=sin(angle)
	co=cos(angle)
	do i=1,npts
	  xrot_temp= x(i)*co+y(i)*si
	  yrot_temp=-x(i)*si+y(i)*co 
	  xrot(i)=xrot_temp
	  yrot(i)=yrot_temp
	enddo
	return
	end
