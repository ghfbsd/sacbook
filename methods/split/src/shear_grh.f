c 	program shear_ev_t
c	pgs 11/1/88
c	pgs 7/7/89 modified to change statistics
c	pgs 1/12/90
c       George Helffrich Aug./Oct. 1992
c          Read event information and data from SAC file.
c          Use SAC exclusively for plotting information
c       George Helffrich Aug. 1994
c          Fix bug in filtering code in changing of length from npts to a power
c          of two.  This would not zero the data area constituting the pad.
c       George Helffrich Aug. 1994
c          Fix bug that would cause blow up if measurement window width was
c          smaller than nshift.  This would usually cause array indexing
c          problems in e_t or e_ev.
c       George Helffrich Aug. 1995
c          Normalize gains before processing unless gains are zero
c       George Helffrich May 1996
c          Return true polarization of incoming wave based on:
c          t-option - r-z particle motion correlation
c          e-option - eigenvector of horizontal particle motion matrix
c       George Helffrich Nov. 1997
c          Smarten up interpolation code so that it won't run out of space
c          while interpolating low-sample-rate data.  The strategy is to
c          shorten the interpolated traces so as to fill up the interpolation
c          array while centering the measurement window in the array.
c       George Helffrich Dec. 1997
c          Fix bug reported by Andrea Restivo that left an extra line of data
c          in the .var file (which wasn't described by the NXSIZE & NYSIZE 
c          header variables.  Harmless, but made files larger than needed.
c          Put additional information into the .var header so that station
c          name and location can be transferred to the variance stack.
c       George Helffrich Nov. 2012
c          Fix legacy Fortran I/O runtime errors.
c
c Copyright (c) 2013 by G. Helffrich.
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c   * Redistributions of source code must retain the above copyright
c     notice, this list of conditions and the following disclaimer.
c   * Redistributions in binary form must reproduce the above copyright
c     notice, this list of conditions and the following disclaimer in
c     the documentation and/or other materials provided with the distribution.  
c   * The names of the authors may not be used to endorse or promote products
c     derived from this software without specific prior written permission.
c THERE IS NO WARRANTY EXPRESSED OR IMPLIED, NOR LIABILITY FOR DAMAGES ASSUMED,
c IN PROVIDING THIS SOFTWARE.
c
c	this version modified allow either minimization of the transverse component or
c	smaller eigenvalue of the covariance matrix
c	program to compute shear splitting parameters
c	It reads in 3-component data in the order radial, transverse, vertical
c	and expects the specification of two basis vectors in terms of colat,lon.
c	colat in degrees measured from vertical up.  Azimuthal angle
c	measured counterclockwise (right handed) from the
c	radial direction(x coordinate).
c	The program solves for the angle theta which minimizes the smallest
c	eigenvalue of the 2-d covariance matrix of horizontal particle motion.
c	(appropriate for studying S phase with unknow polarity) 
c	for several events at the same station. The individual variance
c	arrays are summed and can be used for a joint estimate for a particular
c	station.  The
c	two parameters are the azimuth of the fast polarization direction
c	and the timelag tlag.
c	The two basis vectors are assumed to be orthogonal
c	The conventions are the following:
c	if deltat is found to be negative, it means that u1 is the fast
c	direction and u2 is the slow direction.  
c	V1 is north and V2 is east.  The data are 'derotated
c	to be in this frame.  The angle found is with
c	respect to the original v1.  THis is converted to clockwise azimuth
c	of fast direction with respect to north. 
c	note: error bars correspond to 1sd.
c***    INPUT PARAMETERS
C	Program is run by giving the command shear_ev_t
c	Four arguments are required that must be put on command line
c	1) SAC file names (three of them, E, N, Z) containing data and 
c          event parameters.
c	2) Inversion mode: 't' for minimizing transverse component (appropriate
c	for SKS or 'e' for minimizing smaller eigenvalue of covariance matrix.
c	This is appropriate for other S phases but may also be used for
c	SKS.
c	3) hatndf.  This is the number of degrees freedom per data point.
c	if this parameter is greater than zero, it will be used to compute
c	the ndf.  Typical values are .3 for RSTN ip channel, .1 for DWWSSN
c	ip channel, and .05 to .1 for GSN 20s/s data.  If hatndf is less
c	than or equal to zero, then it is estimated from the data. 
c	4) fourth argument is the name of the measurement file where the
c	resulting measurments will be stored.  If this parameter isn't given,
c       a measurement file name is constructed from the name of the SAC file.
        parameter (ndim=25000,nadim=500)
	parameter (nshift=801,tlaglim=4.)
      	parameter (nangle=181)
	parameter (ncsac=201, aztol = 0.001)
	common/acfcom/acfs1,acfs2,norm
	common/indxcm/lagmin,ianglemin
	logical lev,lndf,ex
	dimension se(ndim), sn(ndim),s(ndim,3),syn(ndim,4)
	dimension ss_ee(1+ndim*2),ss_nn(1+ndim*2)
	dimension ss_en(1+ndim*2),ss_rz(1+ndim*2)
	dimension s1(ndim),s2(ndim),s1_int(ndim),s2_int(ndim),
     +	  s_work(ndim)
	dimension var_array(nshift,nadim),var_array_int(nshift,nadim),
     +	pol_array(nshift,nadim),pol_array_int(nshift)
	dimension tlag_bars(2),angle_bars(2)
	dimension varvec(nshift*nadim)
	equivalence (varvec,var_array_int),(ss_rz,varvec)
c	for testing only. remove when done
c	dimension s_tst(ndim,2)
        real caz(3),cin(3),cgain(3)
	integer*4	iyr,iday,ihr,imin,npts
	character*4     dnet,stn,chan,file(3)
	character charev*1,charndf*10,invmode*1
	character labr*16,labt*16
	character*80 buffer,line*255,flname,
     +  logfl,channel*1,comp*4,dot*1
	character*4 scomp(3)
	character*5 bevn,sphase*7 
	character sacpfx*128
	data lev/.true./
	data pi/3.141592654/
	data in/8/,iout/12/,ioutin/10/,iscr/11/
	data ddegree/1./,dt_int/.05/,dt_long/.25/
	data comp/'prtz'/,dot/'.'/
	data varlim/1.5/

c	debugging IEEE errors.  Set environment variable IEEE_TRAP to
c       activate.
c	call ieeeset('environment')

c	get log file name
	read(*,'(a)')logfl
	  open(unit=7,file=logfl,status='unknown',err=201)

c       Check presence of SAC files.  Read Z component header and fill info
c         from this.
	do i=1,3
	  call getarg(i,line)
	  inquire (file=line,exist=ex)
	  if (.not. ex) then
	    write(0,*) '**',line(1:index(line,' ')),'doesn''t exist,',
     +         'quitting.'
	    stop
	  endif
	enddo
c       call rsach(line,nerr)
c       if (nerr .ne. 0) stop '**Can''t read SAC file.'
	call rsac1(line,s1,npts,beg,dt,0,nerr)
        if (nerr .ne. -803 .and. nerr .ne. 0) then
	   write(0,*) 'nerr = ',nerr
	   stop '**Can''t read SAC file.'
        endif
	call getnhv('NZYEAR',iyr,nerr)
	call getnhv('NZJDAY',iday,nerr)
	call getnhv('NZHOUR',ihr,nerr)
	call getnhv('NZMIN',imin,nerr)
	call getnhv('NZSEC',isec,nerr)
	call getnhv('NZMSEC',imsec,nerr)
	sec = isec + float(imsec)*1e-3
c       Accumulate event information.  Header doesn't have 
c         emb, ems, numsta, moments or mechanism information.
	call getkhv('KSTNM',stn,nerr)
	call getkhv('KCMPNM',chan,nerr)
	call getkhv('KNETWK',dnet,nerr)
	if (dnet .eq. '-123') dnet = 'NONE'
	if (chan .eq. '-123') chan = 'NONE'
	call getfhv('O',origin,nerr)
	call addtim(ieyr,ieday,iehr,iemin,esec,
     +    iyr,iday,ihr,imin,sec,origin)
	call getfhv('EVLA',elat,nerr1)
	call getfhv('EVLO',elon,nerr2)
	call getfhv('EVDP',edep,nerr3)
	if (nerr1 .ne. 0 .or. nerr2 .ne. 0 .or. nerr3 .ne. 0) then
	  stop '**Event information isn''t in file header.'
	endif
	call getfhv('GCARC',edist,nerr1)
	call getfhv('AZ',eaz,nerr2)
	call getfhv('BAZ',ebaz,nerr3)
	if (nerr1 .ne. 0 .or. nerr2 .ne. 0 .or. nerr3 .ne. 0) then
	  call getfhv('STLA',stla,nerr1)
	  call getfhv('STLO',stlo,nerr2)
	  if (nerr1 .ne. 0 .or. nerr2 .ne. 0)
     +      stop '**Station information isn''t in file header.'
	  call disaz(elat,elon,stla,stlo,eaz,ebaz,ekm,edist)
	endif
c       Prefix to generated SAC file names is SAC file name stripped of its
c         .xxx suffix.  This prefix becomes the header "file" name.
        isfx = 0
	ipfx = 0
	do i=lenstr(line),1,-1
	  if (line(i:i) .eq. '.' .and. isfx .eq. 0) isfx = i
	  if (line(i:i) .eq. '/' .and. ipfx .eq. 0) ipfx = i
	enddo
	if (isfx .eq. 0) isfx = lenstr(line)+1
	sacpfx = line(1:isfx-1)
	lsacpf = isfx-1
	file(1) = sacpfx(ipfx+1:ipfx+4)
	file(2) = sacpfx(ipfx+5:ipfx+8)
	file(3) = sacpfx(ipfx+9:ipfx+12)

	charev=' '
	call getarg(4,charev)
	if(charev.eq.'t'.or.charev.eq.'T')then
	  lev=.false.
	  write(*,*)'transverse component will be minimized'
	  write(7,*)'transverse component will be minimized'
	elseif(charev.eq.'e'.or.charev.eq.'E'.or.charev.eq.' ')then
	  lev=.true.
	  write(*,*)'smaller eigenvalue will be minimized'
	  write(7,*)'smaller eigenvalue will be minimized'
	elseif(charev.eq.' ')then
	  write(*,*)'inversion option not designated, abort'
	  stop
	else
	  write(*,*)'invalid option ',charev,' abort'
	  stop
	endif
	charndf(1:1)=' '
	ndf_indx=0
	hatndf_sum=0.
	lndf=.false.
	call getarg(5,charndf)
	if(charndf(1:1).ne.' ')then
	  read(charndf,*,iostat=ios)hatndf
	  if(ios.eq.0 .and. hatndf.gt.0.)then
	    write(*,*)'ndf computed from input value of hatndf=',hatndf
	    lndf=.true.	
	  endif
	endif	
	if(.not.lndf)write(*,*)'ndf computed from the data'	
	rad=180./pi
c	zero var_array,var_array_int 
	do i=1,nshift
	  do j=1,nangle
	    var_array(i,j)=0.
	    var_array_int(i,j)=0.
	  enddo
	enddo
c	Bring data in and open measurement file
	itype=0
c	option='r'
	fastmin=999.
	call getkhv('KSTNM',buffer,nerr)
	len=lenstr(buffer)
	call getfhv('DELTA',delta,nerr)
	if (delta .ge. dt_long) channel = 'l'
	if (delta .lt. dt_long) channel = 'i'
	if (delta .lt. dt_int) channel = 's'
	call getarg(6,flname)        
c	  !allows you name measurement file
	if(flname(1:1).eq.' ')then    
	  flname=buffer(1:len)//dot//channel//'mea'
	endif
	inquire(file=flname,exist=ex)
	if(.not.ex)then
	  open(unit=iout,file=flname,status='unknown',err=204)
	  line = buffer(1:len) // ' measurement file'
	  write(iout,'(a)') line(1:lenstr(line))
	else
	  open(unit=iout,file=flname,status='old',position='append',
     &       err=206)
	endif
6	flname=buffer(1:len)//dot//channel//comp
	write(7,'(1x,a)')flname

c	Read in data - no component inversion necessary for SAC azimuth 
c	convention.  Rotate to radial, tangential for processing.
c       There is no constraint on the orientation of the components except
c       that they differ in orientation by 90 degrees.
	do icomp=1,3
	  call getarg(icomp,line)
	  call rsac1(line,s(1,icomp),npts,beg,dt,ndim,nerr)
	  if (nerr .ne. 0) stop '**Trouble reading SAC input file.'
	  if(icomp.eq.1)then
	    write(6,150) file, iyr, iday, ihr, imin, sec, 
     &         dt, npts, dnet, chan, stn, ' ', scale, offset,
     &         ieyr, ieday, iehr, iemin, esec, elat, elon, edep, 
     &         emb, ems, numstn, edist, eaz, ebaz
	    write(7,150) file, iyr, iday, ihr, imin, sec, 
     &         dt, npts, dnet, chan, stn, ' ', scale, offset,
     &         ieyr, ieday, iehr, iemin, esec, elat, elon, edep, 
     &         emb, ems, numstn, edist, eaz, ebaz
	  endif
	  call getkhv('KCMPNM',scomp(icomp),nerr)
	  call getfhv('CMPAZ',caz(icomp),nerr)
	  call getfhv('CMPINC',cin(icomp),nerr)
	  call getfhv('SCALE',cgain(icomp),nerr)
	  if (nerr .ne. 0) cgain(icomp) = 0.0
	  if (cgain(icomp) .lt. 0.0) then
	     write(*,123) icomp
	     write(7,123) icomp
	  endif
	enddo
c       Normalize gains of all components
	if (cgain(1) .ne. 0.0 .and.
     &      cgain(2) .ne. 0.0 .and.
     &      cgain(3) .ne. 0.0) then
	   do icomp=1,3
	      call scalar(s(1,icomp),npts,1./cgain(icomp))
	   enddo
	else
	   write(*,122)
	   write(7,122)
	endif

c       Sort out whether any rotation to N and E is necessary.  By the end
c       of this operation, sn has N time series, se has E time series, and
c       s(n,3) has vertical.  
	if (cin(3) .eq. 0.0) then
	   i1 = 1
	   i2 = 2
	   ih = 3
	else if (cin(2) .eq. 0.0) then
	   i1 = 1
	   i2 = 3
	   ih = 2
	else if (cin(1) .eq. 0.0) then
	   i1 = 2
	   i2 = 3
	   ih = 1
	else
	   stop '**No vertical component given.'
	endif
	if (abs(cin(i1)-90.0) .gt. aztol .or. 
     +      abs(cin(i2)-90.0) .gt. aztol) 
     +     stop '**"Horizontal" components aren''t horizontal.'
        if (abs(dtheta(caz(i1)/rad,caz(i2)/rad)*rad - 90.0)
     +      .lt. aztol) then
c           1 leads 2 by 90.0 - rotate 1 to E
            if (caz(i2) .ne. 0.0) then
               write(*,130) caz(i2),caz(i1)
               write(7,130) caz(i2),caz(i1)
	    endif
            call rotcomp(s(1,i1),s(1,i2),se,sn,npts,caz(i2)/rad)
        else if (abs(dtheta(caz(i2)/rad,caz(i1)/rad)*rad - 90.0)
     +      .lt. aztol) then
c           2 leads 1 by 90.0 - rotate 2 to E
            if (caz(i1) .ne. 0.0) then
               write(*,130) caz(i1),caz(i2)
               write(7,130) caz(i1),caz(i2)
	    endif
            call rotcomp(s(1,i2),s(1,i1),se,sn,npts,caz(i1)/rad)
	else
	    stop '**Horizontal components not 90 degrees apart.'
	endif
	if (ih .ne. 3) call tscopy(s(1,ih),s(1,3),npts)
	call tscopy(se,s(1,1),npts)
	call tscopy(sn,s(1,2),npts)
	
c       Read in start time relative to file zero and name of phase.

	read(*,*) tloc1,tloc2,sphase    
	loc1 = nint((tloc1-beg)/dt)
	loc2 = nint((tloc2-beg)/dt)
	if (loc1 .gt. npts .or. loc2 .gt. npts .or.
     +      loc1 .lt. 1 .or. loc2 .lt. 1) 
     +    stop '**Data file too large, cut smaller.'
c	  !,sstn,iiyr,iiday
	write(7,*) 'loc1,loc2,sphase:', loc1,loc2,' ',sphase
c	  !,sstn,iiyr,iiday
c*******end of io, now compute**********************************************
c*******Filter option***********************************************
        read(*,*)f1,f2,f3,f4
	if(f1.ne.0 .or. f2.ne.0 .or. f3.ne.0 .or. f4.ne.0) then
	  write(*,131)f1,f2,f3,f4
	  write(7,131)f1,f2,f3,f4
	  if(f4.lt.f3.or.f3.lt.f2.or.f2.lt.f1)
     +        stop ' invalid filter frequencies'
	  ifilter = 1
	else
	  ifilter = 0
	endif
c*******end of filter option
c
c	demean data , detrend data  
 9997   call dtr(s(1,1),npts,2)	
	call dtr(s(1,2),npts,2)	
	call dtr(s(1,3),npts,2)
c*******check for long period************************************************
	if(dt.gt.dt_long)then
c	  long period
	  do icomp=1,3
	    call zeroarray(s_work,1,ndim)
	    call spline_int(s(1,icomp),s_work,dt/dt_long,npts,npts_long)
	    call tscopy(s_work,s(1,icomp),npts_long)
	  enddo
c	  redefine loc1,loc2,npts
	  loc1=float(loc1-1)*dt/dt_long + 1.
	  loc2=float(loc2-1)*dt/dt_long + 1.
	  npts=npts_long
	  dt=dt_long
	  write(*,*)' npts, npts_long, dt, dt_long',npts,npts_long,dt,dt_long
	  write(7,*)' npts, npts_long, dt, dt_long',npts,npts_long,dt,dt_long
	endif 
c*******end of long period check
	if(ifilter.ne.0)then
c*******apply filter***************************************************
	  do icomp=1,3
	   write(*,*)'icomp,ifilter,npts,npts_long',icomp,ifilter,npts,npts_long
	   write(7,*)'icomp,ifilter,npts,npts_long',icomp,ifilter,npts,npts_long
	      length=npts
	      call tscopy(s(1,icomp),s_work,npts)
	      df=1./(dt*length)
	      write(*,*)' dt, df, length', dt,df,length
	      call two(s_work,length)
	      call fftl(s_work,length,-1,ierr)
	      df=1./(dt*length)
	      lenf=length/2+1
	      call filt(s_work,lenf,df,0.,f1,f2,f3,f4,1)
	      call fftl(s_work,length,-2,ierr)
	      call tscopy(s_work,s(1,icomp),npts)
	  enddo
	endif
c*******end of filter block

c       Rotate to change n-e to r-t time series.  s(i,1)->r, s(i,2)->t
	call tscopy(s(1,1),se,npts)
	call tscopy(s(1,2),sn,npts)
	call rotsub(sn,se,s(1,1),s(1,2),npts,ebaz/rad)

	npts1=loc2-loc1+1
	if(.not.lev)then  
c	  !compute tenergy
	  renergy0=0.
	  tenergy0=0.
	  do i=loc1,loc2
	    renergy0=renergy0+s(i,1)**2
	    tenergy0=tenergy0+s(i,2)**2
	  enddo
	  toverr0=tenergy0/renergy0
	endif
	nlag = int(tlaglim/dt+.5)
	nshift1 = 2*nlag + 1
c	check if nshift1 is larger than nshift
	if(nshift1.gt.nshift)then
	  write(*,*)'nshift1 greater than nshift, pgm abort'
	  write(7,*)'nshift1 greater than nshift, pgm abort'
	  stop
	endif
c
c****************************************************************************
c*******Begin of test block**************************************************
c*****************************************************************************
c 	for test, zero out s_tst(1,2) (transverse component) and 
c 	set s(1,1) to s_tst(1,1)
c 	call zeroarray(s_tst(1,2),1,ndim)  
c	  !test
c 	call zeroarray(s_tst(1,1),1,ndim)  
c	  !test
c 	call tscopy(s(1,1),s_tst(1,1),ndim)  
c	  !test
c 	call zeroarray(s(1,1),1,ndim)      
c	  !test
c 	call zeroarray(s(1,2),1,ndim)      
c	  !test
c 	tlag_tst= 1.                      
c	  !test
c 	azimuth_tst=65.                    
c	  !test
c  	add noise to both components. 
c  	get rms signal and set noise to .1 of this	
c 	rms=0.
c 	do kk=1,npts
c 	  rms=rms+s_tst(kk,1)**2
c 	enddo
c 	rms=sqrt(rms)
c 	iseed= 84324
c 	do kk=1,npts
c 	  s_tst(kk,1)=s_tst(kk,1)+.025*rms*gauss(iseed)
c 	  s_tst(kk,2)=s_tst(kk,2)+.025*rms*gauss(iseed)
c 	enddo
c 	call aniput(s_tst(1,1),s_tst(1,2),s(1,1),s(1,2),
c	  !test
c    +	  npts,dt,tlag_tst,azimuth_tst, ebaz)            
c	  !test
c	aniput test ended                               
c	  !test
c***************************************************************************
c*******End of test block*************************************************
c***************************************************************************
c	form cross correlations on subinterval
c	zero out ss, cross cor arrays
 	do kk=1,ndim*2+1
	  ss_ee(kk)=0.
	  ss_nn(kk)=0.
	  ss_en(kk)=0.
 	enddo
c       cross-correlation needs to be at least the length of nshift.
c       make sure of this by checking that npts1 is > (nshift-1)/2.  if
c       not, align cross-correlation info into the arrays ss_-- so that 
c       they begin at least after (nshift-1)/2.
        ixcc = npts1 - (nshift-1)/2
	if (ixcc .lt. 0) then
	  ixcc = 1 - ixcc
	  ixcczl = 1 + (nshift-1)/2
	else
	  ixcc = 1
	  ixcczl = npts1+1
	endif
c	don't normalize
	norm=0
	call crosscorr2(se(loc1),se(loc1),npts1,ss_ee(ixcc))
	call crosscorr2(sn(loc1),sn(loc1),npts1,ss_nn(ixcc))
	call crosscorr2(se(loc1),sn(loc1),npts1,ss_en(ixcc))
	call crosscorr2(s(loc1,1),s(loc1,3),npts1,ss_rz(ixcc))
	ppol = ss_rz(ixcczl)
c	write(*,*)' radial-vertical cross-corr. at zero lag ',ppol
	npts2=npts1*2+1   
c	  !length of the cross correlation time series
c	zero lag is in ixcczl position
	if(lev)then  
c	  !calculate covariance matrix
c	  calculate smaller (of the two) eigenvalues for all angles,tlag
	  call e_ev(ss_ee,ss_nn,ss_en,ixcczl,
     +	    var_array,pol_array,nangle,nshift1,nadim,nshift,ratio)
c	  write(*,*)' ratio in main = ', ratio
c	  write(7,*)' ratio in main = ', ratio
	else 
c	  !transverse component
c	  calculate energy on transverse component for all angles,tlag
	  call e_t(ss_ee,ss_nn,ss_en,ixcczl,
     +	    var_array,nangle,nshift1,nadim,nshift,ebaz)  
	endif
c	interpolate var_array in the tlag direction. check if necessary
	if(dt.lt.dt_int)then
	  dt_int=dt
	  write(*,*)'dt_int set to dt=',dt
	  write(7,*)'dt_int set to dt=',dt
	endif
	if(dt.eq.dt_int)then
	  nshift_int=nshift1
	  do i=1,nangle
	    call tscopy(var_array(1,i),var_array_int(1,i),nshift_int)
	  enddo
	else
	  fac=dt/dt_int
  	  do i=1,nangle	
	    call spline_int(var_array(1,i),var_array_int(1,i),fac,
     +	    nshift1,nshift_int) 
	    if(nshift_int.gt.nshift.and.i.eq.1)then
	      write(*,*)'nshift_int gt nshift ',nshift_int,nshift,'abort'
	      write(7,*)'nshift_int gt nshift ',nshift_int,nshift,'abort'
	      stop
	    endif
	  enddo
	endif
c	zero lag is in the (nshift_int-1)/2+1 position
c	find minimum value and confidence interval	
c	find minimum value of angle tlag. we need this now to get reconstructed
c	seismogram to get number of degrees of freedom.  Will be computed
c	again in shear_bars_all
	call shear_min(var_array_int,nangle,nshift_int,nadim,nshift,
     +	  dt_int,anglemin,tlagmin)   
c	convert anglemin and tlag to the azimuth (clockwise from
c	north and a positive tlag.  Before conversion, anglemin is the
c	counterclockwise angle from east of the fast direction 
c	direction if tlag is negative, and the slow direction tlagmin
c	is positive. thus...
	if(tlagmin.le.0.)then
	  sign=-1.
	  anglemin=90.-anglemin
	else
	  sign=1.
	  anglemin=-anglemin
	  if(anglemin.lt.-90.)anglemin=anglemin+180.
	endif
	tlagmin=abs(tlagmin)
c	reconstruct (radial and transverse) or (pol and pol-perp) 
c	seismograms with lag 
c	first, rotate into fast, slow coordinate system.  Since anglemin
c	is clockwise wrt north, and rotcomp need counterclockwise wrt E
c	we have to send 90-anglemin if x=east and y=north
	call rotcomp(se,sn,s1,s2,npts,(90.-anglemin)/rad)
c	now s1 and s2 are in the appropriate frame
c	interpolate data to dt_int sec. If dt=dt_int, then just copy arrays
	if(dt.eq.dt_int)then
	  call tscopy(s1,s1_int,npts)
	  call tscopy(s2,s2_int,npts)
	  npts_int=npts
	  loc1_int=loc1	
	  loc2_int=loc2	
	else
	  if (npts*dt/dt_int .gt. ndim) then
C           Long trace.  Shift the data so that it fits into a window centered
C           on the measurement window, and is of length ndim.  Adjust file begin
C           time accordingly.
	    nintp = ndim*dt_int/dt
	    locint = max(1,(loc1+loc2 - nintp)/2)
	  else
	    nintp = npts
	    locint = 1
	  endif
	  call spline_int(s1(locint),s1_int,dt/dt_int,nintp,npts_int)
	  call spline_int(s2(locint),s2_int,dt/dt_int,nintp,npts_int)
	  loc1_int=float(loc1-locint)*dt/dt_int + 1.
	  loc2_int=float(loc2-locint)*dt/dt_int + 1.
	  beg = beg + (locint-1)*dt
	endif
	npts1_int=loc2_int-loc1_int+1
c	zero out syn 
	do k=1,4
	  call zeroarray(syn(1,k),1,npts_int) 
	enddo
c	shift seismograms by tlag to get rid of anisotropy and rotate into
c	e-w n-s coord	
c	note, at this stage tlagmin is always positive
     	jshift0=tlagmin/dt_int          
	if(jshift0.lt.0)then
	  write(*,*)'jshift0 is less than 0 ',jshift,' abort'
	  write(7,*)'jshift0 is less than 0 ',jshift,' abort'
	  stop
	endif
        call rotcomp(s1_int,s2_int(1+jshift0),
     +  syn(1,1),syn(1,2),npts_int,(anglemin-90.)/rad)
	call rotcomp(s1_int,s2_int,
     +	syn(1,3),syn(1,4),npts_int,(anglemin-90.)/rad)
	if(lev)then
c*******  find polarization of the minumum variance array elt
	  if(dt.eq.dt_int)then
	    polmin=pol_array(lagmin,ianglemin)
	  else
	    call spline_int(pol_array(1,ianglemin),pol_array_int,
     +	    fac,nshift1,nshift_int)
	    polmin=pol_array_int(lagmin)
	  endif
c	  write(*,*)'sign is ',sign
	  write(*,*)'polarization direction=',polmin,
     +       'degrees(clock from n)'
	  write(7,*)'polarization direction=',polmin,
     +       'degrees(clock from n)'
	  call rotsub(syn(1,2),syn(1,1),syn(1,1),syn(1,2),
     +	  npts_int,(180.+polmin)/rad)
	  call rotsub(syn(1,4),syn(1,3),syn(1,3),syn(1,4),
     +	  npts_int,(180.+polmin)/rad)
	else
	  write(*,*)'rotating into radial transverse'
	  write(7,*)'rotating into radial transverse'
	  call rotsub(syn(1,2),syn(1,1),syn(1,1),syn(1,2),
     +	  npts_int,ebaz/rad)
	  call rotsub(syn(1,4),syn(1,3),syn(1,3),syn(1,4),
     +	  npts_int,ebaz/rad)
	endif

C--------------------------------------------------------------------------
C Code added by Andrea Restivo - Jan 21st, 1998.
C
C Compute signal-to-noise ratio from evaluation of max amplitude of signal
C within measurement window on the 'rl' component and 2 sigma of samples
C distribution over the same (lagged) window on the (noise-only / mean=0)
C tangential 'tl' component. Store it in variable 'snrat' to be transferred
C then into .var file header.
C--------------------------------------------------------------------------

        sqrsmp = 0.
        numsmp = 0
        do i=loc1_int,loc2_int
          sqrsmp = sqrsmp+syn(i,2)**2
          numsmp = numsmp+1 
        enddo
        tnoise = 2*sqrt(sqrsmp/numsmp)

        call mxmn(syn(loc1_int,1),numsmp,ampmn,ampmx) 
        signal = sqrt(abs(ampmx)*abs(ampmn))

        snrat = signal/tnoise 

C--------------------------------------------------------------------------

cc	compute error bars for angle,tlag
	if(.not.lndf)then
       	  ndf=ndf_fun2(syn(loc1_int,2),npts1_int,npts1) 
c	  !pgs
	  write(*,*)'ndf calculated from seismogram'
	  ndf_indx=ndf_indx+1
	  hatndf=float(ndf)/float(npts1)
	  hatndf_sum=hatndf_sum+ hatndf 
	else
	  ndf=hatndf*npts1
	  write(*,*)'ndf calculated from input value'
	endif
c	ftest_mse_2 is specialized to two parameters and 95% conf level.
	varmax=ftest_mse_2(ndf)
	write(*,*)'varmax,npts1,ndf=',varmax,npts1,ndf
	write(7,*)'varmax,npts1,ndf=',varmax,npts1,ndf
	call shear_bars_ev(var_array_int,varmax,
     +	  nangle,nshift_int,nadim,nshift,
     +	  dt_int,anglemin,tlagmin,angle_bars,tlag_bars,varmin)   
c	get average bars and round dt error to multiple of the sample rate
	angle_bar=amax1(abs(angle_bars(1)),abs(angle_bars(2)))
	tlag_bar=amax1(abs(tlag_bars(1)),abs(tlag_bars(2)))
	tlag_bar=dt_int*int((tlag_bar+dt_int/2)/dt_int)
c*******the following has to be done again to make sure nothing is screwed***
c*******up since tlagmin and anglemin have been changed above***************
c	convert anglemin and tlag to the azimuth (clockwise from
c	north and a positive tlag.  Before conversion, anglemin is the
c	clockwise angle from east of the fast direction of the fast
c	direction if tlag is negative, and the slow direction tlagmin
c	is positive. thus...
	if(tlagmin.le.0.)then
	  sign=-1.
	  anglemin=90.-anglemin
	else
	  sign=1.
	  anglemin=-anglemin
	  if(anglemin.lt.-90.)anglemin=anglemin+180.
	endif
	tlagmin=abs(tlagmin)
c*******end of repeated section***********************************************
c	write(*,*)'sign is ',sign
	write(*,*)'anglemin,tlagmin',anglemin,tlagmin
	write(7,*)'anglemin,tlagmin',anglemin,tlagmin
	write(bevn,'(i2.2,i3.3)')mod(iyr,100),iday	
	if(.not.lev)then
	  tenergy_min=varmin/tenergy0
	endif
c       adjust polmin to get proper definition based on rad-vert sense of motion
	if(.not.lev)then  
c	  !min eigenvalue
	  ratio=tenergy_min
	  if (ppol .lt. 0) then
	    polmin=ebaz+180.0
	  else
	    polmin=ebaz
	  endif
	  invmode='T'
	else
c         no need to adjust if using e mode
c	  if (ppol .gt. 0.0) polmin=polmin+180.0
	  invmode='E'
	endif
c	put polmin on interval -90 90
c801	if(polmin.lt.-90..or.polmin.gt. 90.)then
c	  if(polmin.lt.-90.)polmin=polmin+180.
c	  if(polmin.gt. 90.)polmin=polmin-180.
c	  go to 801
c	endif
c	put polmin on interval 0 360
        if(polmin .lt. 0.0)
     +     polmin = polmin + 360.0*(1 + int(polmin/360.0))
        if(polmin .gt. 360.0)
     +     polmin = polmin - 360.0*int(polmin/360.0)
	loc1p = mod(loc1,10000)
	loc2p = mod(loc2,10000)
	call getarg(3,line)
	line = line(ipfx+1:)
	k = lenstr(line)
	write(*,100)bevn,stn,channel,sphase,anglemin,angle_bar,
     +    tlagmin,tlag_bar,
     +	  ratio,ebaz,loc1p,loc2p,polmin,arr_angle,varmax,ndf,hatndf,
     +    invmode,f1,f2,f3,f4,line(1:k)
	write(7,100)bevn,stn,channel,sphase,anglemin,angle_bar,
     +	  tlagmin,tlag_bar,
     +	  ratio,ebaz,loc1p,loc2p,polmin,arr_angle,varmax,ndf,hatndf,
     +    invmode,f1,f2,f3,f4,line(1:k)
	write(iout,100)bevn,stn,channel,sphase,anglemin,angle_bar,
     +    tlagmin,tlag_bar,
     +	  ratio,ebaz,loc1p,loc2p,polmin,arr_angle,varmax,ndf,hatndf,
     +    invmode,f1,f2,f3,f4,line(1:k)
	write(*,121)bevn,stn,channel,sphase,anglemin,angle_bar,
     +	  tlagmin,tlag_bar,ratio,polmin,ebaz

c       Write diagnostic output files.
c       If npts_int really changes here, header B must be updated too.
c       We assume that B time hasn't changed, however.
c       Write out the rotated, lagged signals, radial and tangential.
        write(labr,'(a,f6.1)') 'R',polmin
        write(labt,'(a,f6.1)') 'T',polmin-90.0
	call setfhv('B',beg,nerr)
	call setnhv('NPTS',npts_int,nerr)
	call setfhv('DELTA',dt_int,nerr)
	call setfhv('CMPAZ',polmin,nerr)
	call setkhv('KCMPNM',labr,nerr)
        call wsac0(sacpfx(1:lsacpf)//'.ro',syn(1,3),syn(1,3),nerr)
        call wsac0(sacpfx(1:lsacpf)//'.rl',syn(1,1),syn(1,1),nerr)
	call setkhv('KCMPNM',labt,nerr)
	call setfhv('CMPAZ',polmin-90.0,nerr)
        call wsac0(sacpfx(1:lsacpf)//'.to',syn(1,4),syn(1,4),nerr)
        call wsac0(sacpfx(1:lsacpf)//'.tl',syn(1,2),syn(1,2),nerr)

c	If making the final plot, show time series and 
c	horizontal particle motion
c	in fast, slow frame, both with and without time shift delta t.	
2	continue
c	resample vertical to match horizontals if necessary
        if(dt.ne.dt_int)then
	  call spline_int(s(locint,3),syn(1,3),dt/dt_int,
     +	    nintp,npts_int)
	else
	  call tscopy(s(1,3),syn(1,3),npts)
	endif
c       De-mean the shifted and rotated signal.  Copy the shifted data into
c       another work area so that detrending doesn't introduce a discontinuity
c       in the unshifted signal.
	call dtr(s1_int(loc1_int),npts1_int,1)
        call tscopy(s2_int(loc1_int+jshift0),s_work,npts1_int)
	call dtr(s_work,npts1_int,1)

c       Check whether signal needs inversion by checking sign of 
c       cross-correlation.  If negative, invert.
	s12=0.
	do kk=1,npts1_int
	  s12=s12+s1_int(loc1_int+kk)*s_work(kk)
	enddo
	if(s12.lt.0.)then
	  do kk=1,ndim
	    s2_int(kk)=-s2_int(kk)
	    s_work(kk)=-s_work(kk)
	  enddo
	endif
c	write(*,*)'shifted horizontal-vertical cross corr.',s12
c	for now, only plot the original seismograms in the fast-slow frame
c	and the shifted seismograms.
c	Note: the arrays s1_int and s2_int are the fast and slow directions 
c       respectively.  To shift, simply add jshift0 to the starting array 
c       index for s2_int
        call setfhv('B',tloc1,nerr)
	call setnhv('NPTS',npts1_int,nerr)
	call setfhv('DELTA',dt_int,nerr)
	call setfhv('CMPAZ',polmin,nerr)
	call setkhv('KCMPNM',labr,nerr)
	call wsac0(sacpfx(1:lsacpf)//'.scr',
     +             s1_int(loc1_int),s1_int(loc1_int),nerr)

c       Load SAC file header with measurement results.  These are placed
c       in user0-user6, kuser0-kusern variables:
c          user0 - angle
c          user1 - angle error
c          user2 - time lag
c          user3 - time lag error
c          user4 - polarization direction
c          user5 - # degrees of freedom
c          user6 - # degrees of freedom per data point
        call setfhv('USER0',anglemin,nerr)
        call setfhv('USER1',angle_bar,nerr)
        call setfhv('USER2',tlagmin,nerr)
        call setfhv('USER3',tlag_bar,nerr)
        call setfhv('USER4',polmin,nerr)
        call setfhv('USER5',float(ndf),nerr)
        call setfhv('USER6',hatndf,nerr)
	call setkhv('KCMPNM',labt,nerr)
	call setfhv('CMPAZ',polmin-90.0,nerr)
	call wsac0(sacpfx(1:lsacpf)//'.sct',
     +             s2_int(loc1_int),s2_int(loc1_int),nerr)
	call wsac0(sacpfx(1:lsacpf)//'.scl',
     +             s_work,s_work,nerr)
c	determine extrema of time series, load into header.  
	call mxmn(s1_int(loc1_int),npts1_int,s1min,s1max)
	call mxmn(s2_int(loc1_int),npts1_int,s2min,s2max)
	call setlhv('LEVEN',.false.,nerr)
	call setihv('IFTYPE','IXY',nerr)
	call setfhv('XMAXIMUM',s1max,nerr)
	call setfhv('XMINIMUM',s1min,nerr)
	call setfhv('YMAXIMUM',s2max,nerr)
	call setfhv('YMINIMUM',s2min,nerr)
	call wsac0(sacpfx(1:lsacpf)//'.xy1',
     +             s1_int(loc1_int),s2_int(loc1_int),nerr)
	call mxmn(s_work,npts1_int,s2min,s2max)
	call setfhv('YMAXIMUM',s2max,nerr)
	call setfhv('YMINIMUM',s2min,nerr)
	call wsac0(sacpfx(1:lsacpf)//'.xy2',
     +             s1_int(loc1_int),s_work,nerr)

4	continue
	if(.not.lndf)then
	  hatndf_sum=hatndf_sum/float(ndf_indx)
	  write(7,140) ndf_indx, hatndf_sum
	endif
c	Try writing out the variance file as a SAC 3-D file for contouring.
c       Place in proper format for contour plot.  Write out only positive
c       lag values.  vmult is the variance conversion to units of 95% conf.
c       level.  SAC can't really contour very large datasets.  To prevent it
c       from choking on big ones, the contour plot is interpolated to fewer
c       lags if the number is excessive (larger than ncsac).
        vmult = varmin*varmax
	ix = 1
	n = nshift1-nlag
	nc= min(n,ncsac)
	fnc = float(nc)/float(n)
        if (sign .eq. 0) then
c         First angle corresponds to fast direction of east.  Load
c           contouring array in reverse order.
	  do j=nangle,1,-1
c	    do i=1+nlag,nshift1
c	      varvec(ix) = var_array(i,j)/vmult
c	      ix = ix + 1
c	    enddo
            if (n .eq. nc) then
c             straight copy, no need to interpolate
	      call tscopy(var_array(1+nlag,j),varvec(ix),nc)
	    else
c             interpolate to resample to fewer data points
	      call spline_int(var_array(1+nlag,j),varvec(ix),fnc,n,nc)
	    endif
	    ix = ix + nc
	  enddo
	else
c         First angle corresponds to fast direction of north.  Load
c           contouring array from center (west) in reverse down to
c           the beginning (north), then from the end (south) down to
c           the center (east).
 	  do j=(nangle+1)/2,1,-1
            if (n .eq. nc) then
c             straight copy, no need to interpolate
	      call tscopy(var_array(1+nlag,j),varvec(ix),nc)
	    else
c             interpolate to resample to fewer data points
	      call spline_int(var_array(1+nlag,j),varvec(ix),fnc,n,nc)
	    endif
	    ix = ix + nc
	  enddo
	  do j=nangle-1,(nangle+1)/2,-1
            if (n .eq. nc) then
c             straight copy, no need to interpolate
	      call tscopy(var_array(1+nlag,j),varvec(ix),nc)
	    else
c             interpolate to resample to fewer data points
	      call spline_int(var_array(1+nlag,j),varvec(ix),fnc,n,nc)
	    endif
	    ix = ix + nc
	  enddo
	endif
	call scalar(varvec,ix-1,1./vmult)

	call newhdr
	call setkhv('KSTNM',stn,nerr)
	call setfhv('STLA',stla,nerr)
	call setfhv('STLO',stlo,nerr)
        call setfhv('BAZ',ebaz,nerr)
        call setfhv('USER0',anglemin,nerr)
        call setfhv('USER1',angle_bar,nerr)
        call setfhv('USER2',tlagmin,nerr)
        call setfhv('USER3',tlag_bar,nerr)
        call setfhv('USER4',snrat,nerr)
	call setfhv('USER5',float(ndf),nerr)
	call setfhv('USER6',vmult,nerr)
	call setihv('IFTYPE','IXYZ',nerr)
	call setnhv('NPTS',ix-1,nerr)
	call setnhv('NYSIZE',nangle,nerr)
	call setnhv('NXSIZE',nc,nerr)
	call setfhv('YMINIMUM',-90.0,nerr)
	call setfhv('YMAXIMUM',90.0,nerr)
	call setfhv('XMINIMUM',0.0,nerr)
	call setfhv('XMAXIMUM',tlaglim,nerr)
	call wsac0(sacpfx(1:lsacpf)//'.var',varvec,varvec,nerr)
	stop
201	write(*,*)'error in opening log file, abort'
	stop
202	write(*,*)'error in opening vararray file, abort'
	stop
204	write(*,*)'error in opening mea file, abort'
	stop
205	write(*,*)'error in opening outin file, abort'
	stop
206	write(*,*)'error in opening already exist mea file, abort'
	stop
100	format(1x,a5,1x,a4,1x,a1,1x,a5,2x,f4.0,1x,f4.0,3x,f5.2,
     +  5x,f4.2,2x,e8.2,f5.0,2i4,f6.1,f8.1,f6.2,i8,f6.3,5x,a1,
     +  4(1x,f5.2),1x,a)
110	format(2i10,f5.2/(10f5.2))	
120	format(1x,a5,1x,a4,1x,a1,1x,a5,1x,f4.0,1x,f4.0,2x,f5.2,
c     +	1x,f4.2,1x,f6.3,1x,f5.0,'(',f5.0,')',1x,f5.0,2i4,e13.4,
     +  1x,f4.2,1x,f6.3,1x,f5.0,1x,f5.0,2i4,e13.4,
     +  f10.6,f4.0)
121	format(1x,a5,1x,a4,1x,a1,1x,a5,1x,f4.0,1x,f4.0,2x,f5.2,
     +  1x,f4.2,2x,e8.2,f5.0,1x,f5.0)
122     format(1x,'***Zero gain found, skipping normalization.')
123     format(1x,'***NEGATIVE GAIN FOUND FOR COMPONENT ',i2,'***')
130	format(1x,'Note:  Components not originally N and E (',
     +     f6.1,',',f6.2,')')
131	format(' filter freq:'/4(f6.3,1x))
140	format(' calculated hatndf for ', i5, ' records=',e13.3)
150     format(1x,3a4/1x,i4,1h:,i3,2(1h:,i2),1h:,
     +    f5.2,5x,e12.4,3x,i6/1x,2(a4,1h/),a4,4x,a4,2e12.2,/
     +    1x,'event information:',/,1x,i4,1h:,i3,2(1h:,i2),1h:,
     +    f5.2,5x,2f12.3,f10.1/1x,2f6.1,i10,/,
     +    1x,'distance:',f8.3,1x,'azimuth:',f8.3,1x,'back azimuth:',
     +    f8.3)
	end

	subroutine mxmn(a,n,vmin,vmax)
C       mxmn -- Return extremal values of an array
        real a(n)
	vmin = a(1)
	vmax = a(1)
	do i=2,n
	   vmin=min(a(i),vmin)
	   vmax=max(a(i),vmax)
	enddo
	end

      function dtheta(theta1,theta2)
C     angdif - Find difference between angles expressed in radians
      dsin = sin(theta1)*cos(theta2) - cos(theta1)*sin(theta2)
      dcos = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)
      dtheta = atan2(dsin,dcos)
      end

      function lenstr(string)
      character string*(*)
      nch = len(string)
      do i = 0, nch-1
         ich = nch - i
         if (string(ich:ich) .ne. ' ') then
            lenstr = ich
            return 
         end if
      end do
      lenstr = 0
      return 
      end
