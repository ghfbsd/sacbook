	function f_test2(ndf)
c	this routine calculates the critical value of the
c	variance ratio of two statistical distributions for the
c	95% confidence level by cubic spline interpolation of 
c	table values.  It is specialized to the case of 2 degrees of freedom.
c	returns the critical degrees of freedom up to ndf=999. for
c	ndf gt 999, it prints a caution and then returns the value
c	for 999.
c	It extrapolates using cubic spline interpolation.
c	call spline routines rspln and rsple in utilib
	logical first_time
	parameter (ndim=29,ndfmax=999)
	dimension var_rat(ndim),x(ndim),q(3,ndim),f(3,ndim)
	save first_time,x,var_rat,q
	data first_time/.true./
	data (x(n),var_rat(n),n=1,ndim) /
     +    1.,199.50, 2.,019.00,  3.,009.55,   4.,006.94, 5.,005.79,
     +    6.,005.14, 7.,004.74,  8.,004.46,   9.,004.26,10.,004.10,
     +   11.,003.98,12.,003.89, 13.,003.81,  14.,003.74,15.,003.68,
     +   16.,003.63,17.,003.59, 18.,003.55,  19.,003.52,20.,003.49,
     +   22.,003.44,24.,003.40, 26.,003.37,  28.,003.34,30.,003.32,
     +   40.,003.23,60.,003.15,120.,003.07, 999.,003.00
     +  /

	if(ndf.gt.ndfmax)then
	  write(*,*)'ndf exceeds ndfmax of table',ndf,ndfmax
	  write(*,*)'ndf set to ndfmax'
	  ndf1=ndfmax
	else
	  ndf1=ndf	
	endif
	if(first_time)then
	  call rspln(1,ndim,x,var_rat,q,f)
	  first_time=.false.
	endif
	f_test2=rsple(1,ndim,x,var_rat,q,float(ndf1))
	return
	end
