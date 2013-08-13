      program mtdecon
C
C     MTDECON -- Deconvolve a wavelet from a seismogram using multi-taper
C         spectral cross-correlation.
C
C     Multitaper version engineered by George Helffrich/U. Bristol
C        24 Jan. 2005.
C
C     Takes as command line parameters:
C        -freq x - Frequency domain tapering limit.  Stop frequency is at x Hz.
C        -shift x - shift deconvolved signal by x in time
C        -window ts te - Window in trace destined for deconvolution, start &
C           end time.
C        -wavelet ts te - Wavelet desired for deconvolution, start & end time.
C        -noise ts te - Noise sample between ts and te in wavelet trace.
C
C     Also on command line:
C        wavelet input file name, decon input file name, output file name.
C
C     Alternatively, if none of -window x y, -wavelet x y, -noise x y
C        supplied, then read standard input for file names and windows.  This
C        feature provided for batch processing of files.
C
C     Sets user2 - average squared coherence inside frequency band
C          user3 - fmax value cos**2 taper cutoff frequency
C
C     Copyright G. Helffrich 2005-2013.
C
C Copyright (c) 2013 by G. Helffrich.
C All rights reserved.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are met:
C   * Redistributions of source code must retain the above copyright
C     notice, this list of conditions and the following disclaimer.
C   * Redistributions in binary form must reproduce the above copyright
C     notice, this list of conditions and the following disclaimer in
C     the documentation and/or other materials provided with the distribution.  
C   * The names of the authors may not be used to endorse or promote products
C     derived from this software without specific prior written permission.
C THERE IS NO WARRANTY EXPRESSED OR IMPLIED, NOR LIABILITY FOR DAMAGES ASSUMED,
C IN PROVIDING THIS SOFTWARE.

      parameter(NPMX=14, MAXPOINTS=2**NPMX, MAXDAT=32000)
      parameter(windt=10.0, npta=2**10, nmtw=3, ovlp=0.75)
      real dfft(MAXPOINTS,nmtw), wfft(MAXPOINTS,nmtw)
      real buf(MAXPOINTS), nfft(MAXPOINTS,nmtw)
      real data(MAXDAT), s0(MAXPOINTS), cohsq(MAXPOINTS)
      real*8 tap(npta,nmtw), el(nmtw)
      complex cbuf(MAXPOINTS/2), cnfft(MAXPOINTS/2,nmtw)
      complex cdfft(MAXPOINTS/2,nmtw), cwfft(MAXPOINTS/2,nmtw)
      equivalence (buf,cbuf), (dfft,cdfft), (wfft,cwfft), (nfft,cnfft)
      complex csum
      parameter(nfmax=3, ntmax=9)
      character arg*80,file(nfmax)*80,inlin*256,tok(ntmax)*80
      logical gnum,debug,stdin
      data tdelay/10.0/, ts,te/2*0.0/, tns,tne/2*0.0/, tds,tde/2*0.0/
      data fmax/5.0/, debug/.false./, stdin/.true./, nline/0/

      pi = 4*atan(1.0)
      iskip = 0
      nf = 0
      do 10 i=1,iargc()
         if (i .le. iskip) go to 10
	 call getarg(i,arg)
	 if (arg(1:5) .eq. '-freq') then
	    call getarg(i+1,arg)
	    if (gnum(arg,fmax,'-freq value')) stop
	    iskip = i+1
	 else if (arg .eq. '-shift') then
	    call getarg(i+1,arg)
	    if (gnum(arg,tdelay,'time shift value')) stop
	    iskip = i+1
	 else if (arg .eq. '-debug') then
	    debug = .true.
	 else if (arg .eq. '-window') then
	    call getarg(i+1,arg)
	    if (gnum(arg,tds,'deconvolution start')) stop
	    call getarg(i+2,arg)
	    if (gnum(arg,tde,'deconvolution end')) stop
	    iskip = i+2
	    stdin = .false.
	 else if (arg .eq. '-wavelet') then
	    call getarg(i+1,arg)
	    if (gnum(arg,ts,'wavelet start')) stop
	    call getarg(i+2,arg)
	    if (gnum(arg,te,'wavelet end')) stop
	    iskip = i+2
	    stdin = .false.
	 else if (arg .eq. '-noise') then
	    call getarg(i+1,arg)
	    if (gnum(arg,tns,'noise window start')) stop
	    call getarg(i+2,arg)
	    if (gnum(arg,tne,'noise window end')) stop
	    iskip = i+2
	    stdin = .false.
	 else
	    if (nf .ge. nfmax) then
	       write(0,*) '**Dont understand ',arg(1:index(arg,' ')-1),
     &           ', skipping.'
	    else
	       nf = nf + 1
	       file(nf) = arg
	    endif
         endif
10    continue

C     Many times I find that people have version 1 of Numerical Recipes
C     REALFT, which is not compatible with version 2 used here.  This
C     small check ensures that the right version gets used.
      do i=1,4
	 data(i) = -1
      enddo
      call realft(data,2,+1)
      if (data(3).ne.-1. .or. data(4).ne.-1.)
     &   stop '**WRONG VERSION OF REALFT!!  READ INSTALL NOTES!!'

C     If reading input from standard input, get file names and windows:
C        wavelet 
 
1000  continue

      if (stdin) then
         nline = nline + 1
         read(*,'(a)',iostat=ios) inlin
	 if (ios .ne. 0) stop
	 call tokens(inlin,ntmax,npts,tok)
	 if (npts.ne.ntmax) then
	    write(0,*) '**Bad input, line ',nline,', skipping.'
	    go to 1000
	 endif
C        write(0,'(a,1x,a)') '(mtdecon)',inlin(1:nbl(inlin))
	 file(1) = tok(1)
	 file(2) = tok(4)
	 file(3) = tok(9)
	 inlin = tok(2)(1:nbl(tok(2))) // ' ' //
     &           tok(3)(1:nbl(tok(3))) // ' ' //
     &           tok(5)(1:nbl(tok(5))) // ' ' //
     &           tok(6)(1:nbl(tok(6))) // ' ' //
     &           tok(7)(1:nbl(tok(7))) // ' ' //
     &           tok(8)(1:nbl(tok(8)))
	 read(inlin,*,iostat=ios) tds,tde,ts,te,tns,tne
	 if (ios.ne.0) then
	    write(0,*) '**Bad input, line ',nline,', skipping.'
	    go to 1000
	 endif
      endif

      if (ts .ge. te) then
         write(0,*) '**Nonsense/no wavelet start/end values.'
	 if (stdin) go to 1000
	 stop
      endif
      if (tns .ge. tne) then
         write(0,*) '**Nonsense/no noise start/end values.'
	 if (stdin) go to 1000
	 stop
      endif
      if (tds .ge. tde) then
         write(0,*) '**Nonsense/no decon start/end values.'
	 if (stdin) go to 1000
	 stop
      endif

C     Check existence of trace to deconvolve and determine fft lengths.

      call rsac1(file(2),data,npts,b,dt,0,nerr)
      if (nerr .eq. 0) then
	 if (1+nint((tds-b)/dt) .gt. MAXDAT .or.
     &       1+nint((tde-b)/dt) .gt. MAXDAT)
     &      stop '**Decon file too large.'
	 npts = 1+nint((tde-tds)/dt)
         if (npts .gt. MAXPOINTS) stop '**Decon window too long.'
	 nft = 2
	 do i=1,32
	    if (nft .ge. npts) go to 15
	    nft = 2*nft
	 enddo
	 stop '**Decon input file too long - shorten.'
15       continue
      endif
      if (nerr .ne. 0) stop '**Trouble reading decon input file.'
      do i=1,nmtw
	 call zero(wfft(1,i),1,nft)
	 call zero(dfft(1,i),1,nft)
	 call zero(nfft(1,i),1,nft)
      enddo
      df = 1/(dt*nft)
      nfc = min(nint(fmax/df),nft/2)
      ni0 = 1+2*nfc
      n00 = 1+nft-ni0
      nyq = 1+nft/2

C     Read wavelet file, window with taper and zero the remaining data.

      call rsac1(file(1),data,npts,b,du,MAXDAT,nerr)
      if (nerr .ne. 0) stop '**Trouble reading wavelet input file.'
      if (abs(du-dt)/dt .gt. 1e-3) stop '**Sample rate mismatch.'
      nend = 1+nint((te-b)/dt)
      if (nend .gt. MAXDAT) stop '**File too long - shorten.'
      if (.not.stdin) then
C        Assume synthetic if stdin; zero level set properly
         call demean(data,npts)
         call dtrend(data,npts)
      endif

C     Know sample rate.  Now set up multitaper window values.
      ntap = 1 + 2*(nint(windt/dt)/2)
      ncen = 1 + ntap/2
      if (npta.le.ntap) stop '**Taper storage too small.'
      call mtaper(2.5,ntap,nmtw,npta,el,tap)
      tsum = 0.0
      do j=1,nmtw
	 do i=1,ntap
	    tsum = tsum + tap(i,j)
	 enddo
      enddo
C     ofac is overlap factor within one taper
      ofac=1/(1-ovlp)

C     Extract deconvolution wavelet and obtain fourier spectrum.  Use tapers
C        to window off extraneous part of waveform. 
C     This relies on experiments showing that for 20 sps data with a window
C        size of 10 seconds, the amplitude of the Fourier transformed, tapered
C        rises from zero to full scale in 3 seconds.  Thus we back off 3 sec.
C        from the start of the window to center the first taper through the
C        data.
C
C     nwave is the number of nonzero points in the wavelet that we want
C     nwin  is the number if taper windows needed to cover the wavelet

      nwav = 1 + nint((te-ts+2*0.3*windt)/dt)
      nwin = ofac*nwav/ntap
      is = nint((ts-0.3*windt-b)/dt)

C     Form fft of wavelet.

      do j=1,nwin
	 io = (j-1)*ntap/4
	 do k=1,nmtw
	    do i=1,nft
	       ix=ncen+(i-1)-io
	       id=is+i
	       if (ix.gt.ntap .or. ix.lt.1 .or.
     &             id.lt.1 .or. id.gt.npts) then
		  buf(i) = 0.0
	       else
		  buf(i) = data(id)*tap(ix,k)
	       endif
	    enddo
	    call realft(buf,nft,+1)
	    do i=1,nft
	       wfft(i,k) = wfft(i,k) + buf(i)
	    enddo
	 enddo
      enddo
C     Boost for possibly short window wrt whole trace.
      fac = nwav/float(nwin*nft)
      do k=1,nmtw
	 do i=1,nft
	    wfft(i,k) = wfft(i,k)*fac
	 enddo
      enddo

C     Dump wavelet if wanted.

      if (debug) then
	 write(*,*) 'Wavelet (nwin,nwav):',nwin,nwav
C        shft=ts-0.3*windt-b
         shft=0
         fsarg = 2*pi*shft/(dt*nft)
	 call zero(buf,1,nft)
	 do k=1,nmtw
	    do i=1,nft
	       buf(i) = buf(i) + wfft(i,k)
	    enddo
	 enddo
	 buf(2) = buf(2)*cos(fsarg*nft/2)
	 do i=2,nft/2
	    cbuf(i) = cbuf(i)*exp(cmplx(0,(i-1)*fsarg))
	 enddo
	 call realft(buf,nft,-1)
         fac = tsum/2
	 do i=1,nwav
	    buf(i) = buf(i)/fac
	 enddo
	 call setnhv('npts',nwav,nerr)
         call setfhv('b',ts-0.3*windt,nerr)
	 call wsac0('/tmp/decon.wav',buf,buf,nerr)
      endif

C     Extract noise sample from wavelet trace too.

      nwav = nint((tne-tns+2*0.3*windt)/dt)
      nwin = ofac*nwav/ntap
      is = nint((tns-0.3*windt-b)/dt)

C     Form fft of noise sample.

      do j=1,nwin
	 io = (j-1)*ntap/4
	 do k=1,nmtw
	    do i=1,nft
	       ix=ncen+(i-1)-io
	       id=is+i
	       if (ix.gt.ntap .or. ix.lt.1 .or.
     &             id.lt.1 .or. id.gt.npts) then
		  buf(i) = 0.0
	       else
		  buf(i) = data(id)*tap(ix,k)
	       endif
	    enddo
	    call realft(buf,nft,+1)
	    do i=1,nft
	       nfft(i,k) = nfft(i,k) + buf(i)
	    enddo
	 enddo
      enddo
C     Boost for possibly short window wrt whole trace.  Find power in noise
C     window for normalization factor
      fac = nwav/float(nwin*nft)
      call zero(s0,1,nyq)
      do k=1,nmtw
	 cnfft(1,k) = cnfft(1,k)*fac
	 s0(1) = s0(1) + nfft(1,k)**2
	 s0(nyq) = s0(nyq) + nfft(2,k)**2
	 do i=2,nft/2
	    cnfft(i,k) = cnfft(i,k)*fac
	    s0(i) = s0(i) +
     &         (real(cnfft(i,k))**2 + imag(cnfft(i,k))**2)/el(k)
	 enddo
      enddo

C     Dump noise sample interval if wanted.

      if (debug) then
	 write(*,*) 'Noise (nwin,nwav):',nwin,nwav
C        shft=tns-0.3*windt-b
         shft=0
         fsarg = 2*pi*shft/(dt*nft)
	 call zero(buf,1,nft)
	 do k=1,nmtw
	    do i=1,nft
	       buf(i) = buf(i) + nfft(i,k)
	    enddo
	 enddo
	 buf(2) = buf(2)*cos(fsarg*nft/2)
	 do i=2,nft/2
	    cbuf(i) = cbuf(i)*exp(cmplx(0,(i-1)*fsarg))
	 enddo
	 call realft(buf,nft,-1)
         fac = 2/tsum
	 do i=1,nwav
	    buf(i) = buf(i)*fac
	 enddo
	 call setnhv('npts',nwav,nerr)
         call setfhv('b',tns-0.3*windt,nerr)
	 call wsac0('/tmp/decon.noi',buf,buf,nerr)
      endif

C     Now form fft of trace for deconvolution.

      call rsac1(file(2),data,nptd,b,dtd,MAXDAT,nerr)
      if (nerr .ne. 0) stop '**Deconvolution trace unreadable.'
      if (.not.stdin) then
C        Assume synthetic if stdin; zero level set properly
         call demean(data,nptd)
         call dtrend(data,npts)
      endif
      npts=1+nint((tde-tds+2*0.3*windt)/dt)
      is = nint((tds-0.3*windt-b)/dt)
      nwin = ofac*npts/ntap

      fac = float(npts)/float(nwin*nft)
      do j=1,nwin
	 io = (j-1)*ntap/4
	 do k=1,nmtw
	    do i=1,nft
	       ix=ncen+(i-1)-io
	       id=is+i
	       if (ix.gt.ntap .or. ix.lt.1 .or.
     &             id.lt.1 .or. id.gt.nptd) then
		  buf(i) = 0.0
	       else
		  buf(i) = data(id)*tap(ix,k)
	       endif
	    enddo
	    call realft(buf,nft,+1)
	    do i=1,nft
	       dfft(i,k) = dfft(i,k) + buf(i)
	    enddo
	 enddo
      enddo
      do k=1,nmtw
	 do i=1,nft
	    dfft(i,k) = dfft(i,k)*fac
	 enddo
      enddo

C     Dump deconvolution trace if wanted.

      if (debug) then
	 write(*,*) 'Trace (nwin,nwav):',nwin,npts
	 call zero(buf,1,nft)
	 do k=1,nmtw
	    do i=1,nft
	       buf(i) = buf(i) + dfft(i,k)
	    enddo
	 enddo
	 call realft(buf,nft,-1)
         fac = 2/tsum
	 do i=1,npts
	    buf(i) = buf(i)*fac
	 enddo
	 call setfhv('b',tds-0.3*windt,nerr)
	 call setnhv('npts',npts,nerr)
	 call wsac0('/tmp/decon.trc',buf,buf,nerr)
      endif

C     Now form estimated deconvolution by cross-correlation in spectral
C     domain.  Also estimate spectral coherence.

      call zero(buf,1,nft)
      call zero(cohsq,1,nfc)
      sumfac = 1.0
      rumn1 = 0.0
      rumd1 = 0.0
      cohd1 = 0.0
C     rumn2 = 0.0
C     rumd2 = 0.0
      do k=1,nmtw
	 rumn1 = rumn1 + wfft(1,k)*dfft(1,k)
	 rumd1 = rumd1 + wfft(1,k)*wfft(1,k)
C        rumn2 = rumn2 + wfft(2,k)*dfft(2,k)
C        rumd2 = rumd2 + wfft(2,k)*wfft(2,k)
	 cohd1 = cohd1 + dfft(1,k)*dfft(1,k)
      enddo
      cohsq(1) = rumn1**2/(cohd1*rumd1)
      buf(1) = rumn1/(rumd1 + s0(1))
C     buf(2) = rumn2/(rumd2 + s0(nyq))
      do i=2,nfc
         fac=cos(pi/2*float(i-1)/nfc)**2
	 csum = 0.0
	 sumdw = 0.0
	 sumdd = 0.0
	 do k=1,nmtw
	    csum = csum + cdfft(i,k)*conjg(cwfft(i,k))
	    sumdw = sumdw + real(cwfft(i,k))**2 + imag(cwfft(i,k))**2
	    sumdd = sumdd + real(cdfft(i,k))**2 + imag(cdfft(i,k))**2
	 enddo
	 cbuf(i) = csum*fac/(sumdw + s0(i))
	 cohsq(i) = (real(csum)**2 + imag(csum)**2)/(sumdw*sumdd)
	 sumfac = sumfac + fac
      enddo
      call zero(buf,ni0,n00)
      cohsum = cohsq(1)
      fsarg = 2*pi*tdelay/(dt*nft)
      buf(2) = buf(2)*cos(fsarg*nft/2)
      do i=2,nfc
         fac=cos(pi/2*float(i-1)/nfc)**2
	 cbuf(i) = cbuf(i)*exp(cmplx(0.,fsarg*(i-1)))
	 cohsum = cohsum + fac*cohsq(i)
      enddo
      cohsum = cohsum/sumfac
	
C     Inverse transform the deconvolved trace

      call realft(buf,nft,-1)
C     Normalization factors:
C        nfc/nft - Cos**2 taper, nfc/nft loss of tapered high frequencies
C        nft/2 - nft/2 FFT normalization
      fac = float(nfc)/float(nft) * nft/2
      do i=1,npts
	 buf(i)=buf(i)/fac
      enddo


C     Modify begin time for any time shift, and save water level & Gaussian
C        window width parameters in user2 and user3
      call setfhv('b',-tdelay,nerr)
      call setnhv('npts',npts,nerr)
      call setfhv('user2',cohsum,nerr)
      call setfhv('user3',fmax,nerr)
      call wsac0(file(3),buf,buf,nerr)
      ix = index(file(3),' ')-1
      if (nerr .ne. 0) then
         write(0,*) '**Trouble writing ',file(3)(1:ix)
      else if (stdin) then
         write(*,'(a)') file(3)(1:ix)
      else
	 write(*,*) file(3)(1:ix),': Avg. squared coherence ',cohsum
      endif

C     Write squared coherence at each frequency if wanted.
      if (debug) call wsac1('/tmp/decon.coh',cohsq,nfc,0.0,df,nerr)

      if (stdin) then
         call flush
         go to 1000
      endif
      end

C************************************************************* 

      function nbl(str)
      character str*(*)

      do i=len(str),1,-1
         if (str(i:i) .ne. ' ') then
	    nbl = i
	    return
	 endif
      enddo
      nbl = 1
      end

      subroutine zero(data,ibeg,n)
      dimension data(*)
      do i=ibeg,ibeg+n-1
	 data(i) = 0.0
      enddo
      end

      subroutine demean(a, n)
      dimension a(n)
      ave = 0.
      do 10 i = 1, n
   10 ave = ave + a(i)
      ave = ave / float(n)
      do 20 i = 1, n
   20 a(i) = a(i) - ave
      return
      end

      subroutine dtrend(a, n)
      dimension a(n)
      xx = 0.
      xa = 0.
      cent = float(n + 1) / 2.
      do 10 i = 1, n
      x = float(i) - cent
      xx = xx + (x ** 2)
   10 xa = xa + (x * a(i))
      slope = xa / xx
      do 20 i = 1, n
      x = float(i) - cent
   20 a(i) = a(i) - (slope * x)
      return
      end
