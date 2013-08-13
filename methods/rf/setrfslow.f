C       Program to read in a series of SAC files, potentially recording
C       near-receiver Ps conversions (receiver functions), and to calculate
C       the slowness, optionally putting it into the file header.
C
C       Alternatively, will read a file name and an optional slowness, and
C       will find a P/PcP/ScP or S/ScS arrival to match the given slowness, and
C       will change KA, EVLA, EVLO, EVDP, and, possibly, STLA and STLO, to
C       arrange for the arrival to be at that slowness for slant stack reasons.
C
C       Command line input:
C       -set ffff - set value into file header field named ffff
C       -model xxx - travel time model to use (iasp91 is default).
C       -units [s/km | s/deg]
C       -phase pppp - override info in file header, use pppp as phase
C       -verbose - verbose output
C
C       Reads list of SAC file names from input.  Lines beginning with * or #
C       are skipped.
C
C       George Helffrich/U. Bristol
C          Copyright 2003-2013.

	parameter (Re=6371.0, tol=1e-2)
	parameter (ntmx=2)
	character inline*64, hdr*8, mphs*8, token(ntmx)*8
	real dlim(2), rlim(2)

	logical ok,verbose,okm,omp
	data verbose/.false./, okm,omp /2*.false./, dlim/2*0.0/

	pi = 4.0*atan(1.0)
	hdr = ' '

	iskip = 0
	do 5 i=1,iargc()
	   if (i .le. iskip) go to 5
	   call getarg(i,inline)
	   if ('-model' .eq. inline) then
	      call getarg(i+1,inline)
	      if (inline .eq. ' ') stop '**No model given.'
	      call tpmod(inline)
	      iskip = i+1
	   else if ('-set' .eq. inline) then
	      call getarg(i+1,hdr)
	      if (hdr .eq. ' ') stop '**No header field for -set.'
	      iskip = i+1
	   else if ('-ph' .eq. inline(1:3)) then
	      call getarg(i+1,mphs)
	      if (mphs .eq. ' ') stop '**No named phase for -phase.'
	      omp = .true.
	      iskip = i+1
	   else if ('-unit' .eq. inline(1:5)) then
	      call getarg(i+1,inline)
	      okm = inline .eq. 's/km'
	      iskip = i+1
	   else if ('-verbose' .eq. inline) then
	      verbose = .true.
	   else
	      write(0,*) '**Don''t understand ',
     &          inline(1:index(inline,' ')-1),', skipping.'
	   endif
5       continue

10      continue
	   read (*,'(a)',iostat=ios) inline
	   if (ios .ne. 0) stop
C          Skip comment in input
	   if (0 .ne. index('*#',inline(1:1))) go to 10
	   call tokens(inline, ntmx, ntok, token)
	   if (ntok.gt.1) then
	      ix = index(inline,' ')-1
	   else
	      ix = lenb(inline)
	   endif
	   inquire(file=inline(1:ix),exist=ok)
	   if (.not.ok) then
	      write(0,*) '**File doesn''t exist, skipping: ',
     &          inline(1:index(inline,' ')-1)
	      go to 10
	   endif

	   call rsach(inline(1:ix),nerr)
	   if (nerr .ne. 0) go to 10

	   call getfhv('GCARC',d,nerr)
	   if (nerr .ne. 0) go to 10

	   call getfhv('EVDP',evdp,nerr)
	   if (nerr .ne. 0) go to 10

	   if (.not.omp) then
	      call getkhv('KA',mphs,nerr)
	      if (mphs .eq. '_') call getkhv('KT0',mphs,nerr)
	      if ((mphs .ne. 'P' .and. mphs .ne. 'PP' .and.
     &             mphs .ne. 'PcP' .and.
     &             mphs(1:3) .ne. 'SKS' .and. mphs .ne. 'S' .and.
     &             mphs .ne. 'Pdiff' .and. mphs(1:3) .ne. 'PKP')
     &            .or. nerr .ne. 0) then
	         write(0,*) '**File ',inline(1:ix),
     &           ' has no named main phase',
     &           ' (P or PP or S etc. as KA value in header)',
     &           ' - skipping processing it.'
                 go to 10
	      endif
	   endif

C          What behavior?  If explicit slowness given, hunt for range of
C             phase with that slowness.
	   if (ntok.gt.1) then
	      read(token(2),*,iostat=ios) p
	      if (ios.ne.0) then
	         write(0,*) '**File ',inline(1:ix),
     &           ' has bad explicit slowness - skipping.'
                 go to 10
	      endif
	      call getfhv('stla',stla,nerr1)
	      call getfhv('stlo',stla,nerr2)
	      if (nerr1 .ne. 0 .or. nerr2 .ne. 0) then
	         write(0,*) '**File ',inline(1:ix),
     &              ' has no station info, setting to zero.'
                 stla = 0.001
		 stlo = 0.001
		 call setfhv('stla',stla,nerr)
		 call setfhv('stlo',stlo,nerr)
	      endif
C             If not found, seek boundary where diffraction starts.
	      ix = index('PS',mphs(1:1))
	      if (dlim(ix) .eq. 0.0) then
	         rlo = 85
		 rhi = 110
		 do i=1,10
		    rng = 0.5*(rlo+rhi)
		    tt = tpttxx(mphs(1:1),rng,0.0,dtdd,dtdh,.true.)
		    if (tt.lt.0) then
		       rhi = rng
		    else
		       rlo = rng
		    endif
	         enddo
		 tt = tpttxx(mphs(1:1),rlo,0.0,dtdd,dtdh,.true.)
		 if (tt.lt.0) pause '**Bad diffracted search'
		 dlim(ix) = dtdd
		 rlim(ix) = rlo
	      endif
	      if (p .lt. dlim(ix)) then
	         mphs(2:) = 'c' // mphs(1:1)
		 rlo = 0.0
	      else
	         rlo = 25.0
	      endif
	      rhi = rlim(ix)
	      do i=1,20
	         rng = 0.5*(rlo+rhi)
		 tt = tpttxx(mphs,rng,0.0,dtdd,dtdh,.false.)
		 dif = abs(dtdd-p)
		 if (dif .lt. tol) exit
		 if (dtdd .lt. p) then
		    if (mphs(2:2) .eq. 'c') then
		       rlo = rng
		    else
		       rhi = rng
		    endif
	         else
		    if (mphs(2:2) .eq. 'c') then
		       rhi = rng
		    else
		       rlo = rng
		    endif
	         endif
	      enddo
	      call dazell(stla,stlo,rng,90.0,evla,evlo)
	      call setfhv('evla',evla,nerr)
	      call setfhv('evlo',evlo,nerr)
	      call setfhv('evdp',0.0,nerr)
	      call setfhv('gcarc',rng,nerr)
	      call setkhv('ka',mphs,nerr)
	      call wsach(inline(1:index(inline,' ')),nerr)
	      if (nerr .ne. 0) then
	         write(0,'(3a)') '**',inline(1:ix),
     &              ':  Trouble writing header.'
              endif
           else
C             Get travel time of main phase.
	      ttP = tpttxx(mphs,d,evdp,dtdd,dtdh,.false.)
	      if (okm) dtdd = 180/(pi*Re)*dtdd

C             Print info and get next file if that's all that's wanted.
	      ix = lenb(inline)
	      if (verbose) then
	         write(*,*) inline(1:ix),' ',mphs(1:index(mphs,' ')-1),
     &              ' ',dtdd
	      endif

C             Set value into file header if required.
              if (hdr .ne. ' ') then
	         call setfhv(hdr,dtdd,nerr)
	         if (nerr .ne. 0) stop '**Bad header variable name used.'
	         call wsach(inline,nerr)
	         if (nerr .ne. 0) then
	            write(0,'(3a)') '**',inline(1:ix),
     &                 ':  Trouble writing header.'
                 endif
	      endif
	   endif

        go to 10

999     continue
	end

	function tpttxx(id,delta,depth,dtddel,dtddep,exact)
C       tpttxx -- Version of tptt that handles multiple arrivals and
C          returns the first one that matches.  This is hopefully the
C          phase of interest.
	parameter (nmax=30)
	real tt(nmax),dtdd(nmax),dtdh(nmax)
	character idphs(nmax)*8, id*(*)
	logical exact

        n = mtptt('all',delta,depth,nmax,idphs,tt,dtdd,dtdh)
	do 10 i=1,min(n,nmax)
	   if (id .eq. idphs(i)) go to 12
	   if (.not. exact) then
	      if (id .eq. 'P' .and. idphs(i) .eq. 'Pg') go to 12
	      if (id .eq. 'S' .and. idphs(i) .eq. 'Sg') go to 12
	      if (id .eq. 'P' .and. idphs(i) .eq. 'Pb') go to 12
	      if (id .eq. 'S' .and. idphs(i) .eq. 'Sb') go to 12
	      if (id .eq. 'P' .and. idphs(i) .eq. 'Pn') go to 12
	      if (id .eq. 'S' .and. idphs(i) .eq. 'Sn') go to 12
	      if (id .eq. 'P' .and. idphs(i) .eq. 'Pdiff') go to 12
	      if (id .eq. 'S' .and. idphs(i) .eq. 'Sdiff') go to 12
	   endif
10      continue
C	write(0,*) 'TPTTXX:  No ',id,' arrival at ',delta
        tpttxx = -1.0
        return

12      continue
        dtddel = dtdd(i)
	dtddep = dtdh(i)
	tpttxx = tt(i)
	end

	integer function lenb(string)
	character*(*) string

	do 10 i=len(string),2,-1
	   if (string(i:i) .ne. ' ') go to 12
10      continue
12      continue
	lenb = i
	end
