C  -----------------------------------------------------------------------------
C  Program to stack a sequence of 3-D SAC files of splitting measurement
C  variance.
C
C  Program arguments:
C
C     [-jack] infile outfile
C
C  where infile is a list of .var files, outfile is an output file for the
C  joint estimate, and the optional parameter -jack indicates that the input
C  is a set of files already processed with this program which comprise a
C  jackknife estimates of the collected data (and lack back-azimuth and
C  S/N info for weighting).
C
C  Individual files are weighted according to signal-to-noise ratio and the
C  azimuthal sampling of raypaths at the station.
C  Read filenames from input file one by one. Read data into an array and
C  header variable USER5 (number of degrees of freedom) into real variable fndf.
C  Compute depmin (lambda2.min) of individual files by means of the subroutine
C  "mxmi" (Can't get value of header variable DEPMIN out of SAC).
C  Normalize the data values to each file's depmin and store the progressive
C  sum of jth elements in the jth slot of an auxiliary array.
C  Keep updating the value and stkndf (sum of degrees of freedom from each file
C  - converted to an integer value for purpose of using function ftest - ) to be
C  able to return it to the corresponding header variable in the output file.
C  Calculate also the value of stkmin (lambda2.min of stacked file) by "mxmi".
C  Compute the variance ratio for 95% confidence in the stacked file, as a
C  function of the total number of degrees of freedom, invoking the function
C  ftest_mse_2.
C  Normalize then the stacked values to minimum (stkmin) = 1 and scale them to
C  the 95% confidence limit value = 1.
C  Read this array into a matrix with indexes (ixphi, ixdt) expressing the
C  position on a phi-dt grid (-90/90 deg; 0/4 sec), to retrieve then the
C  corresponding phi and dt associated with the minimum value of lambda2 and to
C  quantify the standard deviation (1 sigma) of both variables around that
C  point. Convert angles to 0/180 deg grid in the case that the 95% confidence
C  area extends across the border of former grid where -90=90.
C  Stacked values are finally written to an output 3-D sac file (together with
C  the new header variables values - stkndf to be re-converted to a real value
C  to be stored in floating point variable USER5 - ) ready to be displayed as a
C  contour plot in SAC.
C
C  Original version by A. Restivo, April 1998
C
C  Modified to handle different sample rates by G. Helffrich, July 1999
C  Modified to fix bugs in azimuthal weight calculation, G. Helffrich, Oct. 2000
C  Modified to add weighting of degrees of freedom, G. Helffrich, Oct. 2000
C  Modified to include misc. routines (mxmi, sort), G. Helffrich, Oct. 2004
C  Modified to do jackknifing, G. Helffrich, 25 Mar. 2005
C  Modified to do without rsach, G. Helffrich, 4 Aug. 2011
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
C  -----------------------------------------------------------------------------

      program sstack

      parameter (maxp=40000)
      parameter (mxnx=201)
      parameter (tol=1e-4)
      real values(maxp),stack(maxp)
      real ebaz(mxnx),dbaz(mxnx)
      character infile*256,sacnam*256,arg*16
      logical ojack
      data ojack/.false./, noarg/0/

      infile = ' '
      sacnam = ' '
      iskip = 0
      do 5 i=1,iargc()
	 if (i .le. iskip) go to 5
	 call getarg(i,arg)
	 if (arg .eq. '-jack') then
	    ojack = .true.
	 else if (infile .eq. ' ') then
	    call getarg(i,infile)
	 else if (sacnam .eq. ' ') then
	    noarg = i
	 else
	    write(0,*) '**',arg(1:index(arg,' ')),'is invalid, ignoring.'
	 endif
5     continue

      do i=1,mxnx
         ebaz(i) = 0
         dbaz(i) = 0
      enddo
      stkndf = 0

      if (infile.eq.' ' .or. noarg.eq.0)
     &   stop '**No input or output file'
      open (1,file=infile)

      j = 0
      do 10 i=1,mxnx
         read (1,'(a)',iostat=ios) sacnam

         if (ios .ne. 0) go to 15
	 if (0 .ne. index('*#',sacnam(1:1))) go to 10
	 j = j + 1

C        call rsach(sacnam,nerr)
C             if (nerr .ne. 0) go to 999
         call rsac1(sacnam,values,nlen,beg,del,0,nerr)
              if (nerr .ne. 0 .and. nerr .ne. -803) go to 999
         call getfhv('BAZ',ebaz(j),nerr)
              if (nerr .ne. 0) go to 999
	 call getnhv('NXSIZE',n,nerr)
              if (nerr .ne. 0) go to 999
	 call getfhv('XMINIMUM',xmin,nerr)
              if (nerr .ne. 0) go to 999
	 call getfhv('XMAXIMUM',xmax,nerr)
              if (nerr .ne. 0) go to 999
	 sumwth = (xmax-xmin)/(n-1)
	 if (j .eq. 1) then
	    sdt = sumwth
	    nmax = n
	 else if (abs(sdt-sumwth) .gt. tol) then
	    sdt = min(sdt,sumwth)
	    nmax = max(n,nmax)
	 endif

         if (ebaz(j) .ge. 180.0) ebaz(j) = ebaz(j)-180.0
10    continue
      write(0,*) '**Too many data files, first ',mxnx,' used.'

15    continue
      close(1)

C  -----------------------------------------------------------------------------
C  End of first reading; all ebaz values read and stored, and minimum sample
C  rate determined (std).
C  Variable nfile keeps record of number of splitting measurements.
C  Now second reading with values as well and updating stack(i).
C  Determination of weights therefore to be performed within this section.
C  -----------------------------------------------------------------------------

      nfile = j

      do i=1,maxp

         stack(i) = 0

      enddo

      sumwth = 0

      open (1,file=infile)

      do n=1,nfile
20       continue
	    read (1,'(a)',iostat=ios) sacnam
	    if (ios .ne. 0) then
	       write(0,*) '**Error, can''t re-read file "',
     +            sacnam(1:index(sacnam,' ')-1),'".'
               stop
	    endif
	 if (0 .ne. index('*#',sacnam(1:1))) go to 20

         call rsac1(sacnam,values,nlen,beg,del,maxp,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('BAZ',baz,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('USER4',snrat,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('USER5',fndf,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('USER6',vnorm,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('XMINIMUM',xmin,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('XMAXIMUM',xmax,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('YMINIMUM',ymin,nerr)
              if (nerr .ne. 0) go to 999
         call getfhv('YMAXIMUM',ymax,nerr)
              if (nerr .ne. 0) go to 999
         call getnhv('NXSIZE',nxsize,nerr)
              if (nerr .ne. 0) go to 999
         call getnhv('NYSIZE',nysize,nerr)
              if (nerr .ne. 0) go to 999
	 call mxmi(nlen,values,depmin,depmax)

         if (nxsize .gt. mxnx) then
            stop '**ERROR: sampling rate exceeds max.!'
         endif
	 if (nysize*nmax .gt. maxp) then
	    stop '**Stack storage too small -- recompile.'
	 endif
	 dt = (xmax-xmin)/(nxsize-1)

C  -----------------------------------------------------------------------------
C  Azimuthal and S/N sampling weight determination.
C  -----------------------------------------------------------------------------
         
	 if (ojack) then
	    azwght = 1.0
	    snwght = 1.0
	 else
            if (baz .ge. 180.0) baz = baz-180.0

            do i=1,nfile
               dbaz(i) = abs(ebaz(i)-baz)
               if (dbaz(i) .gt. 90.0) dbaz(i) = 180.0-dbaz(i)
            enddo

            if (nfile.gt.1) call sort(nfile,dbaz)

            do i=1,nfile
               if (dbaz(i) .gt. 10.0) go to 50
            enddo
50          continue

            azwght = 1.0/(i-1.0)
            snwght = funwth(snrat)
	 endif

C  -----------------------------------------------------------------------------
C  Weighted stack to be next performed.
C  -----------------------------------------------------------------------------

	 k = 0
         do 110 j=1,nysize
	    do 111 i=1,nmax
	       k = k + 1
C              Interpolate to proper intermediate lag or pick matching time
	       if (abs(dt-sdt) .gt. tol) then
	          t = xmin + (i-1)*sdt
		  val = wigint(xmin,values(1+(j-1)*nxsize),nxsize,dt,
     &                         0.0,t)
	       else
		  val = values(k)
	       endif
               stack(k) = stack(k)+val/depmin*azwght*snwght
111         continue
110      continue

         stkndf = stkndf+fndf*snwght*azwght
         sumwth = sumwth+snwght*azwght
      enddo

      close(1)

      nlen = nysize*nmax
      call scalar(stack,nlen,1./sumwth)

C     Calculate weighted degrees of freedom by taking weighted mean and
C        multiplying by the number of files.  If jackknife, use average
C        DF of files, not total DF in all files.
      if (ojack) then
	 stkndf = stkndf/sumwth
      else
         stkndf = stkndf/sumwth*nfile
      endif

C  -----------------------------------------------------------------------------
C  End of stacking procedure; now building of uncertainty plot. 
C  -----------------------------------------------------------------------------

      call mxmi(nlen,stack,stkmin,stkmax)
      varmax = ftest_mse_2(nint(stkndf))

      ixmin = 1
      do j=1,nlen
         stack(j) = stack(j)/(stkmin*varmax)
	 if (stack(j) .lt. stack(ixmin)) ixmin = j
      enddo

C     Find minimum value and calculate best-fit phi and dt
      stkmin = stack(ixmin)
      poldir = ymin + int((ixmin-1)/nmax)*(ymax-ymin)/(nysize-1)
      dtlag =  xmin + mod(ixmin-1,nmax)*sdt
      twosig = 1.0

C     Pass through all data values and find max phi difference with values
C        <= two sigma from minimum value.
C
C     Phidif is 2*max angle difference (in radians); dtmax and dtmin are min
C        and max dt.
      phidif = 0.0
      dtmax = dtlag
      dtmin = dtlag
      rad = 4.0*atan(1.0)/180.0
      sdphi = (ymax-ymin)/(nysize-1)

      do j=1,nlen
         if (stack(j) .le. twosig) then
	    angle = ymin + int((j-1)/nmax)*sdphi
	    phidif = max(phidif,acos(cos(2*rad*(angle-poldir))))

	    delay =  xmin + mod(j-1,nmax)*sdt
	    dtmax = max(dtmax,delay)
	    dtmin = min(dtmin,delay)
	 endif
      enddo

C     We report 1-sigma values, halving the 2-sigma intervals.
      poldev = min(sdphi*nint(phidif/(4*rad*sdphi)),22.5)
      dtdev = nint(max(dtlag-dtmin,dtmax-dtlag)/sdt)*sdt/2

      call scalar(stack,nlen,1./twosig)

C  -----------------------------------------------------------------------------
C  Writing of sac output file.
C  -----------------------------------------------------------------------------

      call getarg(noarg,sacnam)

      call setnhv('NPTS',nmax*nysize,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('DEPMIN',stkmin,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('USER5',stkndf,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('USER4',poldir,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('USER3',poldev,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('USER2',dtlag,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('USER1',dtdev,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('XMINIMUM',xmin,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('XMAXIMUM',xmax,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('YMINIMUM',ymin,nerr)
           if (nerr .ne. 0) go to 999
      call setfhv('YMAXIMUM',ymax,nerr)
           if (nerr .ne. 0) go to 999
      call setnhv('NXSIZE',nmax,nerr)
           if (nerr .ne. 0) go to 999
      call setnhv('NYSIZE',nysize,nerr)
           if (nerr .ne. 0) go to 999 
      call wsac0(sacnam,stack,stack,nerr)
           if (nerr .ne. 0) stop '**Unable to write output file.'

999   continue

      end




C  -----------------------------------------------------------------------------
C  Function funwth defines weight to give to samples with a determined S/N rat.
C  It is defined so that samples below S/N* = 1.0 are given a weight that from
C  0.01 tends asynthotically to 0.0, while samples over S/N* = 21.0 are given a
C  weight that from 0.99 tends to 1.0. Values of S/N* in between return a weight
C  which steadily increases with S/N.
C  S/N* indicates the S/N ratio calculated by the shear program (biased higher).
C  -----------------------------------------------------------------------------

      function funwth(x)

      real k,mu

      epstop = 21.0
      epsbot = 1.0
      width = epstop-epsbot
      mu = width/2
      funetp = 0.99
      funebt = 0.01

      k = (log(funebt**2/funetp**2))/width

      funwth = 1/(exp(k*(x-mu))+1)

      end

      subroutine mxmi(n,vals,vmin,vmax)
      real vals(n)

      if (n .ge. 1) then
	 vmin = vals(1)
	 vmax = vals(1)
      endif
      do i=2,n
	 v = vals(i)
	 if (v .lt. vmin) vmin = v
	 if (v .gt. vmax) vmax = v
      enddo
      end

      SUBROUTINE SORT(N,RA)
      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
