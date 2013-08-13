	function ftest_mse_2(n)
c	a version of ftest_mse that is specialized to the case of 2 parameters
c	the 95% confidence level.  It computes the value of the ftest by
c	interpolation of a table.  Not recommended above n=120.
c	this function provides the critical value of the ratio of mean
c	squared misfit function at the minimum to any other value, such
c	the values inside the region defined by ratio<ftest_mse_2 corresponds
c	to a confidence interval at the 95% confidence level.  
c	'n' is the "effective" number
c	of data, ie the number of degrees of freedom, assumed to be an integer.
	ftest_mse_2= 1. + 2.*f_test2(n-2)/float(n-2)
	return
	end
