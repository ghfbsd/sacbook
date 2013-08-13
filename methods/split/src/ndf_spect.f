	Function ndf_spect(a,nf)
c+	This function calculates the number of degrees of freedom
c	appropriate for a summed and squared time series, whose spectrum
c	is array A. A is assumed to be a complex function of frequency with 
c	A(1) corresponding to zero frequency and A(nf) corresponding to
c	the Nyquist frequency.  If A(t) is gaussian distributed
c	then ndf_spect should be the points in the original
c	timer series 2*(nf-1).
c
c       Original code by P. Silver modified to accumulate double precision
c       sums to prevent underflow and overflow.  (G. Helffrich/U. Bristol)
c-
	complex a(*)
	double precision temp, f2, f4

c	for zero frequency and for Nyquist:
c***	note the following line is from Silver and Chan 1990.
c	but this appears to give a slightly too large value
c	for the number of df for a gaussian distributed time series.
c	using 1/2 instead of 1/3 seems to give the right answer
c	but I dont know why.  It depends on the fourth moment of
c	the zero frequency and nyquist frequency points.

	temp = dble(cabs(a(1)))**2 + dble(cabs(a(nf)))**2
	f2=0.5d0*temp
	temp = dble(cabs(a(1)))**4 + dble(cabs(a(nf)))**4
	f4=0.5d0*temp
	do i=2,nf-1
	  temp=dble(cabs(a(i)))
	  f2=f2+temp**2
	  f4=f4+temp**4
 	enddo
c	based on theory, the following expression should yield
c	the correct number of degrees of freedom.(see appendix of silver
c	and chan, 1990.  In practise, it gives a value that gives
c	a df that is slightly too large.  eg 101 for an original time
c	series of 100.
	ndf_spect=2.d0*(2.d0*f2**2/f4 -1.d0)
	return
	end
