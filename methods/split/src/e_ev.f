	subroutine e_ev(ss_ee,ss_nn,ss_en,n0,
     +	  e_ev_array,pol_array,na,nlag,na1,nlag1,ratio)
c	subroutine to compute energy on minimum-eigenvalue component for angles
c	0 through 180 and for lags of +-jshift (defined below)
c	ss... are the cross and auto correlations of the original
c	e-w, n-s seismograms (rotated form the radial and transverse)and n0 is
c	the index of the zero-lag point.
	dimension e_ev_array(nlag1,na1),pol_array(nlag1,na1),
     +  ss_ee(n0+nlag),ss_nn(n0+nlag),ss_en(n0+nlag),c(2,2)
	data pi/3.141592654/
	rad=180./pi
c	check if nlag is odd
	if(mod(nlag,2).eq.0)then
	  write(*,*)'nlag must be odd, pgm aborted: nlag',nlag
	  call exit
	endif
	jshift=(nlag-1)/2
	alammin=1.e+22
	do iangle=0,180,1
	  theta=float(iangle)/rad
	  c2 =cos(theta)**2
	  s2 =sin(theta)**2
	  sc =sin(theta)*cos(theta) 
c	  calculate lag-independent quantities
	  ss_epep = ss_ee(n0)*c2 + 2.*ss_en(n0)*sc + ss_nn(n0)*s2
	  ss_npnp = ss_ee(n0)*s2 - 2.*ss_en(n0)*sc + ss_nn(n0)*c2
	  ilag=0
	  do lag=n0-jshift,n0+jshift,1	  
	    ilag=ilag+1	
c	    calculate lag dependent term
	    ss_epnp = -ss_ee(lag)*sc - ss_en(2*n0-lag)*s2 +
     +                 ss_en(lag)*c2+ ss_nn(lag)*sc
c	    matrix elements
	    c(1,1) = ss_epep 
	    c(2,2) = ss_npnp 
     	    c(1,2) = ss_epnp 
	    c(2,1)=c(1,2)
	    trace= c(1,1) + c(2,2)
	    euc  = c(1,2)**2 - c(1,1)*c(2,2)
	    alam1 = .5*(trace+sqrt(trace**2+4.*euc))	    	    
	    alam2 = .5*(trace-sqrt(trace**2+4.*euc))	    	    
c	    write(*,*)'alam1, alam2 = ',alam1,alam2
	    trace_check=ss_ee(n0)+ss_nn(n0)
	    e_ev_array(ilag,iangle+1) = alam2
c	    polarization in azimuth degrees(clockwise from north)
c	    note, vectors are in 'primed' coord frame and must add
c	    theta to the angle.
	    pol_array(ilag,iangle+1)=
     +	      90.-rad*(atan2((alam1-c(1,1)),c(1,2))+theta)		    
	    if(alam2.lt.alammin)then
	      alammin=alam2
	      rmin=alam2/alam1
	      minlag = lag-n0
	      minang = iangle
c           else if(alam2.eq.alammin)then
c             write(*,9)minlag,minang,lag-n0,iangle
9             format(' same lambda min found - was ',2(i4,1x),
     +               ', now is ',2(i4,1x))
	    endif
	    if(lag.eq.n0.and.iangle.eq.0)rorig=alam2/alam1
	    if(e_ev_array(ilag,iangle+1).lt.0.)then
	      write(*,*)'e_ev_array lt 0:',e_ev_array(ilag,iangle+1)
	      write(*,*)'lag,iangle,e_ev_0',lag,iangle,e_ev_0
	      write(*,*)'alam1,alam2,trace,euc',
     +	        alam1,alam2,trace,euc
	    endif
	  enddo
	enddo
c       prevents zero division if both (and probably all) are zero
	if (rmin .eq. rorig) then
	  ratio = 1.0
        else
	  ratio=rmin/rorig
        endif
c	write(*,*)'ratio in e_ev',rmin,rorig,ratio
c       show pol_array value at complementary position to minimum
	mclag = -minlag
	if (minang .gt. 90) then
	   mcang = minang - 90
	else
	   mcang = minang + 90
	endif
c	write(*,*) 'Min. lambda2 ',
c    +     rmin,e_ev_array(minlag+jshift+1,minang+1),
c    +     ' complement ',e_ev_array(mclag+jshift+1,mcang+1)
c	if (alammin .ne. e_ev_array(mclag+jshift+1,mcang+1))
c    +     write(*,10) alammin,e_ev_array(mclag+jshift+1,mcang+1)
10      format('**Note:  Min. phi, dt lambda isn''t the same as its ',
     +     'complement: ',2g12.6)
	return
	end
