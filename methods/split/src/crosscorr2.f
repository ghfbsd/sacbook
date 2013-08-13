	subroutine crosscorr2(s1,s2,n,ss)
c	routine to compute the cross correlation function
c	between two time series of length n.  It is normalized
c	by sqrt(acf1(0)*acf2(0)) if norm =1, or not if norm=0.
c	It is normalized by default.
c	acf1(0) and acf2(0) are stored in common
c  	The result ss is of length 2n+1.  The zero-lag point
c	is in position n+1. ss(1) and ss(2n+1) are set to zero.
c	This is a test of shearsp by doing the cc in time domain
c	Assumes the signal is zero outside the boundaries of the array
c	
	dimension s1(n),s2(n),ss(2*n+1)
	common/acfcom/acfs1,acfs2,norm
	data norm/1/
	n2=2*n
	indx=1
	ss(1)=0.
	ss(n2)=0.
	do lag=-(n-1),n-1,1
	  indx=indx+1
	  sum = 0.0
c	  find overlap
	  iover=n-iabs(lag)
	  if(lag.le.0)then
	    do i=1,iover
	      sum = sum + s1(i)*s2(i-lag)
	    enddo
	  else
	    do i=1,iover
	      sum = sum + s1(i+lag)*s2(i)
	    enddo
	  endif
	  ss(indx) = sum
	enddo
	acfs1=0.
	acfs2=0.
	do i=1,n	
	  acfs1=acfs1+s1(i)**2
	  acfs2=acfs2+s2(i)**2
	enddo
	fac=sqrt(acfs1*acfs2)
c	normalize cc
	if(norm.eq.1)then
	  write(*,*)'cross correlation normalized'	  
	  call scalar(ss,n2,1./fac)	
	else
	  write(*,*)'cross correlation not normalized'
	endif
	return
	end
                               
