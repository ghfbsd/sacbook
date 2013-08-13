	subroutine spline_int(y,y_int,fac,npts,npts_int)
c	subroutine to cubic-spline-interpolate y by a factor of fac.
c	y assumed to be equally spaced in x.
c	calls rspln,rsple
	parameter (ndim=25000)
	dimension y(npts),y_int(*)
	dimension q(3,ndim),f(3,ndim),x(ndim)
	if(npts.gt.ndim)then
	  write(*,*)'error,you have reached limit: ndim,npts',ndim,npts
	  call exit
	endif
c	fill up x array
	do i=1,npts
	  x(i)=float(i-1)
	enddo
	call rspln(1,npts,x,y,q,f)
	dx=1./fac	
	npts_int=float(npts-1)*fac+1.
	do i=1,npts_int
	  s=float(i-1)*dx
	  y_int(i)=rsple(1,npts,x,y,q,s)
	enddo
	return
	end
	subroutine unspline(y,fac,npts_int)
c	subroutine to decimate data array y by nfac. Put back in y
	dimension y(npts_int)
	indx=0
	do i=1,npts_int,nfac
	  indx=indx+1
	  y(indx)= y(i)
	enddo
	return
	end
