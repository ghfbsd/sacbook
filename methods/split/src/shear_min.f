	subroutine shear_min(var_array,na,nlag,na1,nlag1,
     +	dt,anglemin,tlagmin)    
	dimension var_array(nlag1,na1) 
	common/indxcm/lagmin,ianglemin
c	finds anglemin,tlagmin corresponding to minimum value of transverse
c	component.
c	assumes that nlag is odd =2*jshift+1, and that the zero-lag point
c	is in the jshift+1 position.  This is not checked so be careful
c	Also assumes that first angle is 0 deg and last is 180 (1 degree
c	increments

	jshift=(nlag-1)/2+1
c	look for minimum point of var_array
	varmin=1e20 
	do i=1,nlag
	  do j=1,na
	    if(var_array(i,j).lt.varmin)then
	      varmin=var_array(i,j)
	      lagmin=i
	      ianglemin=j		
	    endif
	  enddo
	enddo     
	tlagmin=(lagmin-jshift)*dt
	anglemin=float(ianglemin-1)
	return
	end
