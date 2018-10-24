C     Replacement routine for Numerical Recipes realft() routine, modified
C     from routines by A. Shakal to take/return data in storage format
C     specified in Numerical Recipes.  The difference is that the first
C     complex element (c(1); two real elements (x(1),x(2))) of the FFT contains
C     the DC component (x(1) or real(c(1))) and the Nyquist component
C     (x(2) or imag(c(1))).  Because the imaginary component of both complex
C     FFT values are zero for the FFT of a real series, the two imaginary zero
C     values are omitted and the real values are packed into c(1).
C
C     By G. Helffrich/U. Bristol
C        12 Aug. 2013, 24 Oct. 2018

      subroutine realft(data,n,isign)
      dimension data(n)
      if (isign.eq.1) then
         call dfftr(data,n,'forward',1.0)
      else
         call dfftr(data,n,'inverse',0.5)
      endif
      end

      subroutine dfftr (x,nft,dirctn,delta)
c                                              a.shakal, 1/78, 15 jul 80
c           this subroutine does a fast fourier transform on a real
c        time series.  it requires 1/2 the storage and e1/2 the time
c        required by a complex fft.
c
c     forward transform, "call dfftr(x,nft,'forward',dt)":
c           input = x(1),x(2),..,x(nft) = real time series of nft points
c          output = x(1),x(2),..,x(nft) = nft/2+1 complex spectral points
c        these spectral points are identical to the first nft/2+1 returned
c        by subroutine fft (i.e., pos freq terms).  thus, the coefficient
c        at fj, the j-th frequency point (where fj = (j-1)*delf, j=1,nft
c        and delf = 1/(nft*dt)), is in x(i-1),x(i), where i=2j.  x(1) is
c        dc term (no imag part because real time series), x(2) is real part
c        of nyquist coef (no imaginary part because real series).
c
c     inverse transform, "call dfftr(x,nft,'inverse',delf)":
c        input and output are interchanged.
c
c           if this subroutine is called with 'forward', and then with 'inverse'
c        and delf of 1/(nft*dt), the original time series is recovered.
c        identical results (but for scaling) can be obtained by calling
c        fft(x,nft,isign), but in fft a real time series must be stored
c        complex array with zero imaginary parts, which requires 2*nft p
c        of array x.  also, the coefs returned by the fft will differ by
c        n-scaling, since fft's leave out the dt,delf of the approximate
c        integrations.  this subroutine calls fft.
c           this subroutine is a modification of the subroutine 'fftr',
c        written by c.frasier.  the principal modifications are:
c             1) the delt,delf of the integrations are included to make
c                a discrete approximation to the fourier transform.
c
      logical forwrd, invrse
      character dirctn*7
      complex  csign, c1, c2, c3, speci, specj
      real x(nft)
      pi = 3.1415927
c
      call locast(dirctn,invrse,forwrd)
c
      nftby2 = nft/2
      if (forwrd) then
c            forward transform..
         call fft (x,nftby2,+1)
         sign = +1.
      else if (invrse) then
         sign = -1.
      else
         stop 'dirctn bad to dfftr'
      endif
c
c           manipulate elements as appropropriate for a 1/2 length
c        complex fft, after the forward fft, or before the inverse.
      piovrn = pi*sign/float(nftby2)
      csign = cmplx(0.,sign)
      do 10 i = 3,nftby2,2
      j = nft-i+2
      c1 = cmplx(x(i)+x(j), x(i+1)-x(j+1))
      c2 = cmplx(x(i)-x(j), x(i+1)+x(j+1))
      w = piovrn*float(i/2)
      c3 = cmplx(cos(w),sin(w))*c2
      speci = c1 + csign*c3
      x(j) = real(speci)/2.
      x(j+1) = -aimag(speci)/2.
      specj = conjg(c1) + csign*conjg(c3)
      x(i) = real(specj)/2.
      x(i+1) = -aimag(specj)/2.
   10 continue
      if (forwrd) then
c            include dt of integration, for forward transform...
         tmp=x(1)
	 x(1)=tmp+x(2)
	 x(2)=tmp-x(2)
         dt = delta
         do i = 1,nft
            x(i) = x(i)*dt
	 enddo
      else
         tmp=x(1)
	 x(1)=(tmp+x(2))/2
	 x(2)=(tmp-x(2))/2
c            do the inverse transform...
         call fft (x,nftby2,-1)
c            in the inverse transform, include the df of the integration
c            and a factor of 2 because only doing half the integration
c            (i.e., just over the positive freqs).
         twodf = 2.*delta
         do i = 1,nft
            x(i) = x(i)*twodf
	 enddo
      endif
      end

      subroutine locast(dirctn,invrse,forwrd)
      character dirctn*7
      logical forwrd,invrse
      if(dirctn.eq.'forward') go to 1
      if(dirctn.eq.'inverse') go to 2
      write(0,100)dirctn
  100 format(1x,a7,2x,'is meaningless to dfftr, use forward or inverse
     *only')
      invrse=.false.
      forwrd=.false.
      return
    1 invrse=.false.
      forwrd=.true.
      return
    2 invrse=.true.
      forwrd=.false.
      return
      end

      subroutine fft(data,nn,isign)
c                                              a.shakal, 1/78, 10 jul 80
c        cooley-tukey 'fast fourier trnasform' in ansi fortran 77.
c
c           transform(j) = sum {data(i)*w**u(i-1)*(j-1)e}, where i and
c        j run from 1 to nn, and w = exp(sign*twopi*sqrtu-1e/nn).
c        data is a one-dimensional complex array (i.e., the real and
c        imaginary parts of the data are located immediately adjacent
c        in storage, such as fortran places them) whose length nn is
c        a power of two.  isign is +1 or -1, giving the sign of the
c        transform.  transform values are returned in array data,
c        replacing the input data.  the time is proportional to
c        n*log2(n), rather than the non-fft n**2.  modified from the
c        fortran ii coding from n.brenner's mit-ll tech rept.
c
      real data(2*nn)
      pi = 3.1415926
c
      n = 2*nn
      j = 1
      do 5 i = 1,n,2
      if (j.gt.i) then
         tempr = data(j)
         tempi = data(j+1)
         data(j) = data(i)
         data(j+1) = data(i+1)
         data(i) = tempr
         data(i+1) = tempi
      endif
      m = n/2
    1 if (m.ge.2 .and. j.gt.m) then
         j = j-m
         m = m/2
	 go to 1
      endif
      j = j+m
    5 continue
c
c
      mmax = 2
    6 if (n .gt. mmax) then
         istep = 2*mmax
         pibymx = pi*float(isign)/float(mmax)
c
         do 8 m = 1,mmax,2
         theta = pibymx*float(m-1)
         wr = cos(theta)
         wi = sin(theta)
         do 8 i = m,n,istep
         j = i + mmax
         tempr = wr*data(j) - wi*data(j+1)
         tempi = wr*data(j+1) + wi*data(j)
         data(j) = data(i) - tempr
         data(j+1) = data(i+1) - tempi
         data(i) = data(i) + tempr
         data(i+1) = data(i+1) + tempi
   8     continue
         mmax = istep
         go to 6
      endif
      return
      end
