	subroutine e_t(ss_ee,ss_nn,ss_en,n0,
     +	  e_t_array,na,nlag,na1,nlag1,ebaz)
c	subroutine to compute energy on shear component for angles
c	0 through 180 and for lags of +-jshift (defined below)
c	ss... are the cross and auto correlations and n0 is
c	the index of the zero-lag point.
	dimension e_t_array(nlag1,na1),ss_ee(n0+nlag),
     +	ss_nn(n0+nlag),ss_en(n0+nlag)
	data pi/3.141592654/
	rad=180./pi
c	check if nlag is odd
	if(mod(nlag,2).eq.0)then
	  write(*,*)'nlag must be odd, pgm aborted: nlag',nlag
	  call exit
	endif
	jshift=(nlag-1)/2
	do iangle=0,180,1
	  theta=float(iangle)/rad
	  baz=ebaz/rad
	  co_t =cos(theta)
	  si_t =sin(theta)
	  co_tb=cos(theta+baz)
	  si_tb=sin(theta+baz)
	  co_2tb=cos(2.*theta+baz)  
c	  !for a test
	  si_2tb=sin(2.*theta+baz)  
c	  !for a test
	  e_t_0=
     +	    co_tb**2*(ss_ee(n0)*co_t**2   +
     +	              ss_nn(n0)*si_t**2   +				
     +	           2.*ss_en(n0)*si_t*co_t)    +
     +	    si_tb**2*(ss_ee(n0)*si_t**2   +
     +	              ss_nn(n0)*co_t**2   -				
     +	           2.*ss_en(n0)*si_t*co_t) 
	  ilag=0
	  do lag=n0-jshift,n0+jshift,1	  
	    ilag=ilag+1	
	    e_t_array(ilag,iangle+1)=e_t_0-
     +	    2.*si_tb*co_tb*((-ss_ee(lag)+ss_nn(lag))*si_t*co_t   +
     +	                  ss_en(lag)*co_t**2   -	
     +	                  ss_en(2*n0-lag)*si_t**2)   
	    if(e_t_array(ilag,iangle+1).lt.0.)then
	      write(*,*)'e_t_array lt 0:',e_t_array(ilag,iangle+1)
	      write(*,*)'lag,iangle,e_t_0',lag,iangle,e_t_0
	      if(lag.eq.n0)then
	        e_t_test=ss_ee(n0)*co_2tb**2 +
     +	                 ss_nn(n0)*si_2tb**2 +
     +		      2.*ss_en(n0)*si_2tb*co_2tb
	        write(*,*)'lag=n0, recalculating using simple formula'
		write(*,*)'e_t_test=',e_t_test
	        write(*,*)'ss_ee,ss_nn,ss_en',
     +	          ss_ee(n0),ss_nn(n0),ss_en(n0)
	      endif
	    endif
	  enddo
	enddo
	return
	end
