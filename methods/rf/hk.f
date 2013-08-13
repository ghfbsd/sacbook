C     Do H-K stacking of receiver functions in a collection of SAC files,
C     aligning to a specified time in each.
C
C     Input:
C        Command line has file name for output.
C        -vp x - Base Vp value.
C        -h lo hi n - range of or H values desired, and number of discrete
C            values.
C        -k lo hi n - range of kappa values desired, and number of discrete
C            values.
C        -weight wt1 wt2, wt3 - Weights for t1 (Ps), t2 (PpPs) and t3 (PsPs).
C        -phaseweight n - phase weight exponent.
C        -p header-name [ s/km | s/deg ] - header variable name where slowness
C            available, and optional units (default s/km).  Needed for
C            calculation of t1, t2 and t3.
C        -amp - preserve amplitude in stack
C        -info - write out distance info for each trace
C  
C     Std input has list of file names.
C
C     Output XYZ file has standard deviation in USER0.
C
C     George Helffrich, U. Bristol
C        Copyright 2003-2013.
C
C     Modifications:
C        - account for 180 degree phase reversal of PsPs.
C        - calculate analytic uncertainty in H and k.
C        - units for slowness value.
C        - Fix problem in t3 calculation.
C        - Report projected uncertainty in H and k.
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
      parameter (np2=17,nmax=2**np2,nsmax=2*32768,clfac=2.0)
      real w(3),t(3),s(3),data(nmax),hdata(nmax),stack(nsmax)
      real hlo, hhi, klo, khi, kappa
      equivalence (t1,t(1)),(t2,t(2)),(t3,t(3))
      integer ndat(nsmax)
      complex ci,csum
      parameter (ntmax=3)
      character sacfile*128,arg*128,token(ntmax)*128,phval*8,trcid*16
      logical gnum, opws, oinfo, goth, gotk, oamp
      data oinfo /.false./, goth, gotk, oamp/3*.false./
      data opws,root /.true.,1.0/, vp/0.0/
      data eps/0.0/, ci/(0.0,1.0)/, s/1.0,1.0,-1.0/
      data wsum, w/1.0, 0.7,0.2,0.1/, pfac/1.0/

C     Read data files and parameters.
C     call ieeeset('environment')
      pi = 4.0*atan(1.0)
      phval = ' '
      sacfile = ' '
      iskip = 0
      do i=1,iargc()
	 if (i .le. iskip) cycle
	 call getarg(i,arg)
	 if (arg(1:1) .ne. '-') then
	    sacfile = arg
	 else if (arg .eq. '-h') then
	    call getarg(i+1,arg)
	    if (gnum(arg,hlo,'low H value')) stop
	    call getarg(i+2,arg)
	    if (gnum(arg,hhi,'high H value')) stop
	    call getarg(i+3,arg)
	    if (gnum(arg,v,'H value count')) stop
	    nh = nint(v)
	    iskip = i+3
	    if (hlo .ge. hhi) stop '**H lo value > high?'
            if (hlo.le.0 .or. nh.le.0) stop '**Nonsense -h range.'
	    goth = .true.
	 else if (arg .eq. '-k') then
	    call getarg(i+1,arg)
	    if (gnum(arg,klo,'low K value')) stop
	    call getarg(i+2,arg)
	    if (gnum(arg,khi,'high K value')) stop
	    call getarg(i+3,arg)
	    if (gnum(arg,v,'K value count')) stop
	    nk = nint(v)
	    iskip = i+3
	    if (klo .ge. khi) stop '**K lo value > high?'
            if (klo.le.0 .or. nk.le.0) stop '**Nonsense -k range.'
	    gotk = .true.
	 else if (arg .eq. '-amp') then
	    oamp = .true.
	 else if (arg .eq. '-p') then
	    call getarg(i+1,phval)
	    call getarg(i+2,arg)
	    if (arg.eq.'s/km') then
	       pfac = 1.0
	       iskip = i+2
	    else if (arg.eq.'s/deg') then
	       pfac = 1.0/111.19493
	       iskip = i+2
	    else
	       iskip = i+1
	    endif
	 else if (arg .eq. '-vp') then
	    call getarg(i+1,arg)
	    if (gnum(arg,vp,'Vp value')) stop
	    iskip = i+1
	 else if (arg(1:3) .eq. '-ph') then
	    call getarg(i+1,arg)
	    if (gnum(arg,root,'phaseweight exponent')) stop
	    iskip = i+1
	    opws = .true.
	 else if (arg(1:2) .eq. '-w') then
	    wsum = 0
	    do j=1,3
	       call getarg(i+j,arg)
	       if (gnum(arg,w(j),'weight value')) stop
	       wsum = wsum + w(j)
	    enddo
	    if (wsum .le. 0) stop '**Bad weights (sum <=0)'
	    iskip = i+3
	 else if (arg .eq. '-info') then
	    oinfo = .true.
	 else
	    write(0,*) '**Unrecognized: ',arg(1:index(arg,' '))
	 endif
      enddo
      if (sacfile .eq. ' ') stop '**No output file given.'
      if (.not.goth .or. .not.gotk) stop '**-h or -k not given.'
      if (nh*nk .gt. nmax) stop '**Too many H and K values.'
      if (vp .le. 0) stop '**No -vp given (or < 0).'
      if (phval .eq. ' ') stop '**No -p given.'

C     Zero stack
      ix = 0
      do i=1,nh
         do j=1,nk
	    ix = ix + 1
	    stack(ix) = 0.0
	 enddo
      enddo

C     Add signs to weights
      psum = 0
      do i=1,3
         w(i) = w(i)*s(i)
	 if (w(i).gt.0) psum = psum + 1
      enddo
      psum = max(1.,psum)

      dk = (khi-klo)/max(1,nk-1)
      dh = (hhi-hlo)/max(1,nh-1)

C     Start reading files.
      n = 0
1000  format(a,1x,f8.5,2(1x,f8.2),1x,1pe8.2)
100   continue
	 read(5,'(a)',iostat=ios) arg
	 if (ios .ne. 0) go to 200
	 if (arg(1:1) .eq. '*') go to 100
	 call tokens(arg,ntmax,nt,token)
	 if (nt .lt. 1) then
	    write(0,*) '**Missing info in input.'
	    go to 100
	 endif
	 ix = nblen(token(1))
         call rsac1(token(1),data,npts,d0,dt,nmax,nerr)
	 if (nerr .ne. 0) then
	    write(0,*) '**Can''t read ',token(1)(1:ix),', skipping.'
	    go to 100
	 endif
	 call getfhv(phval,p,nerr)
	 if (nerr .ne. 0) then
	    write(0,*) '**Header value "',phval(1:nblen(phval)),'" not ',
     &         'in header for ',token(1)(1:ix),', skipping.'
            go to 100
	 endif
	 p = p*pfac
	 n = n + 1
	 call hilbtf(data,hdata,npts,nmax)
	 dn = d0 + (npts-1)*dt

	 if (oinfo) then
	    write(*,1000) token(1)(1:ix),p,d0,dn,dt
	 endif

	 is = 1
	 nx = 0
	 tp = sqrt(1/vp**2 - p**2)
         do i=1,nh
	    h=hlo + float(i-1)*dh
	    do j=1,nk
	       kappa=klo + float(j-1)*dk
	       vs = vp/kappa
	       ts = sqrt(1/vs**2 - p**2)

	       t1 =   h*(ts - tp)
	       t2 =   h*(ts + tp)
	       t3 = 2*h*ts
	       csum = 0
	       sum = 0
	       fac = 1
	       do k=1,3
	          if (t(k).ge.d0 .and. t(k).le.dn) then
		     v = wigint(d0,data,npts,dt,eps,t(k))
		     if (opws) then
			c = wigint(d0,hdata,npts,dt,eps,t(k))
			if(w(k).ne.0)then
			   if (w(k).gt.0) then
			      csum = csum + exp(ci*atan2(c,v))
			   else
			      csum = csum + exp(ci*atan2(-c,-v))
			   endif
			endif
		     endif
		     sum = sum + w(k)*v
		  else
		     if (w(k).ne.0) nx = nx + 1
		  endif
	       enddo
	       if (opws) fac = abs(csum/psum)**root
	       sum = sum/wsum * fac

	       stack(is) = stack(is) + sum
	       is = is + 1
	    enddo
	 enddo
	 if (nx.gt.0) then
	    v = p*6371*pi/180
	    write(0,*) '**',token(1)(1:ix),': ',nx,
     &         ' H-K comb. not in trace, p ',v,'sec/deg.'
	 endif
      go to 100

200   continue
      write(*,*) n,' files processed.'
      n = max(n,1)
      nb = nh*nk

C     Normalize stack by number of seismograms; calculate power.
      sum = 0
      srn = 0
      do i=1,nb
         if (oamp) then
            stack(i) = stack(i)/n
	 else
            stack(i) = abs(stack(i)/n)
	 endif
	 sum = sum + stack(i)**2
	 srn = max(srn,stack(i))
      enddo
      sstd = sqrt(sum/nb)
      if (srn.eq.0) srn = 1

C     Scale stack to 0 - 1, find max.
      ixm = 1
      do i=1,nb
	 stack(i) = stack(i)/srn
	 if (stack(i) .gt. stack(ixm)) ixm=i
      enddo

C     Estimate second deriv. at max. in H and k directions.
      d2sdk2 = -(stack(max(1,ixm-1))+stack(ixm+1) - 2*stack(ixm))/dk**2
      if (ixm-nk.lt.1) then
         ixdm = ixm
      else
         ixdm = ixm-nk
      endif
      if (ixm+nk.gt.nb) then
         ixdp = ixm
      else
         ixdp = ixm+nk
      endif
      d2sdh2 = -(stack(ixdm)+stack(ixdp) - 2*stack(ixm))/dh**2
      sdk = sstd*sqrt(2/d2sdk2)
      sdh = sstd*sqrt(2/d2sdh2)

C     Find projection on H and k axes of 1 std. dev. down from peak value
      ixm = ixm-1
      ikmx = mod(ixm,nk)
      ikmn = ikmx
      ihmx = ixm/nk
      ihmn = ihmx
      srn = 1-sstd
      do i=1,nb
         if (stack(i).gt.srn) then
	    ik = mod(i-1,nk)
	    ih = (i-1)/nk
	    if (ik.gt.ikmx) ikmx=ik
	    if (ik.lt.ikmn) ikmn=ik
	    if (ih.gt.ihmx) ihmx=ih
	    if (ih.lt.ihmn) ihmn=ih
	 endif
      enddo
      h = hlo + (ixm/nk)*dh
      kappa = klo + mod(ixm,nk)*dk
      prh = 0.5*(ihmx-ihmn)*dh
      prk = 0.5*(ikmx-ikmn)*dk
      write(*,*) 'Max value: ',h,kappa
      write(*,*) ' h +/- ',sdh,' k +/- ',sdk,', analytic'
      write(*,*) ' h +/- ',prh,' k +/- ',prk,', projected'
      write(*,*) 'Stack RMS amp.: ',sstd

C     Write stack file
      if (opws) then
	 write(trcid,'(a,i2)') 'PWS',int(root)
      else
         trcid = 'Linear'
      endif
      call newhdr
      call setnhv('NPTS',nh*nk,nerr)
      call setfhv('B',0.0,nerr)
      call setfhv('DELTA',1.0,nerr)
      call setkhv('KSTNM','STACK',nerr)
      call setkhv('KEVNM',trcid,nerr)
      call setihv('IFTYPE','IXYZ',nerr)
      call setnhv('NXSIZE',nk,nerr)
      call setnhv('NYSIZE',nh,nerr)
      call setfhv('XMINIMUM',klo,nerr)
      call setfhv('XMAXIMUM',khi,nerr)
      call setfhv('YMINIMUM',hlo,nerr)
      call setfhv('YMAXIMUM',hhi,nerr)
      call setfhv('USER0',sstd,nerr)
      call wsac0(sacfile,stack,stack,nerr)
      end

      function nblen(string)
      character string*(*)
      do i=len(string),1,-1
	 if (string(i:i) .ne. ' ') go to 1
      enddo
      i = 1
1     continue
      nblen = i
      end
