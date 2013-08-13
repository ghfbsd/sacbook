       	function ndf_fun2(a,n,norig) 
c+	A array is assumed to contain a possibly interpolated time series
c	of length n.  The length of the original uninterpolated time series
c	is norig.  ndf_spect computes the effective
c	number of degrees of freedom in a, which for a gaussian white noise
c	process should be equal to norig.  The true ndf should be less than
c	norig.
	parameter (ndim=10000)
	dimension a(*), b(ndim)
c***	following two lines for test, remove when done
	complex b_comp(ndim/2)
	equivalence (b,b_comp)
c***    end of test
	if(n.gt.ndim-2)then
	  write(*,*)'n bigger than nmax, abort', n,ndim-2
	  call exit
	endif
	nadd=n
c	zero out work array
1	call zeroarray(b,1,ndim)
	call tscopy(a,b,n)
	call fftl(b,nadd,-1,ierr)
	if(ierr.eq.1)then
	  write(*,*)'fft not factorable, add a point'
	  nadd=nadd+1
	  go to 1 
	endif
c	determine number of frequency points for original time series.
c
	nf=(norig+2)/2
c*****************************
c*******test block
c	as a test, apply a filter to the spectrum
c 	write(*,*)'in ndf_fun2, testing with a filter. Please remove'
c	do i=1,nf
c	  b_comp(i)=b_comp(i)*float(i-1)
c 	enddo
c*******end of test block
	ndf_fun2=ndf_spect(b,nf)
	write(*,*)'orig data pts, ndf:', norig,ndf_fun2
	if(ndf_fun2.gt.norig)then
	  write(*,*)'caution:ndf gt norig'
	endif
	return
	end
