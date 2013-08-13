	subroutine shear_bars_ev(var_array,varmax,na,nlag,na1,nlag1,
     +	dt,anglemin,tlagmin,angle_bars,tlag_bars,varmin)    
c	version of shear_bars to use with shearmain_ev
	dimension var_array(nlag1,na1),tlag_bars(2),angle_bars(2)
c	common/varplt/tlagpl(15000),anglepl(15000),npoints
c	common to return array indices of minimum variance location
	common/indxcm/lagmin,ianglemin
	logical first_time
	data eps/1e-10/
c	assumes that nlag is odd =2*jshift+1, and that the zero-lag point
c	is in the jshift+1 position.  This is not checked so be careful
c	Also assumes that first angle is 0 deg and last is 180 (1 degree
c	increments
c	note, the error bars correspond to 1 standard deviation in that
c	twice the error bars corresponds to the 95% confidence interval.

c	initialize npoints
	npoints=1
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
c	normalize var_array by varmin for ftest
c	check for varmin=0 (arises in synthetic cases) and set to eps
	if(varmin.eq.0.)varmin=eps
	do i=1,na 
	  call scalar(var_array(1,i),nlag,1./varmin)
	enddo
c	tlagpl(npoints)=tlagmin 
c	anglepl(npoints)=anglemin 
	first_time=.true.
	jshift=(nlag-1)/2+1
	iangle1=ianglemin-45
	iangle2=ianglemin+45
	do i=iangle1,iangle2
c	  check for end effects
	  if(i.lt.1)then
	    ii=180+i
	  elseif(i.gt.180)then
	    ii=i-180
	  else
	    ii=i
	  endif
	  angle0=i-1
	  do j=1,nlag 
	    if(var_array(j,ii).lt.varmax)then
	      if(first_time)then
	        amin=angle0
	        amax=amin
	        nlagmin=j 
	        nlagmax=nlagmin
	        first_time=.false.
	      else
		if(angle0.lt.amin)amin=angle0
		if(angle0.gt.amax)amax=angle0
	        nlagmin=min0(nlagmin,j)
	        nlagmax=max0(nlagmax,j)
	      endif
c	      npoints=npoints+1
c	      if(npoints.gt.15000)then
c	        write(*,*)'npoints has gotten too large, abort'
c	        call exit
c	      endif
c	      tlagpl(npoints)=dt*(j-jshift-1) 
c	      anglepl(npoints)=angle0
	    endif    
	  enddo                  
	enddo
	tlag_bars(1)= .5*((nlagmin-jshift-1)*dt-tlagmin)
	tlag_bars(2)= .5*(tlagmin-(nlagmax-jshift-1)*dt)
	angle_bars(1)=.5*( amin-anglemin)
	angle_bars(2)=.5*(anglemin-amax)
	write(*,*)'tlagmin,anglemin',tlagmin,anglemin
	write(*,*) 'tlag_bars,angle_bars',tlag_bars,angle_bars
	return
	end
