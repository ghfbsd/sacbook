C     RSACH, WSACH -- Read and write SAC file headers.  This subroutine should
C        exist in SAC, but in it isn't supplied in a user-callable form.
C
C     This subroutine relies on SAC using a common block called 'cmhdr'.
C     We check that this is the case by calling newhdr and looking at its
C     effect on 'cmhdr'.  If it isn't as predicted, we assume things are
C     messed up and abort.
C
C     RSACH -- Read SAC header.
C
C        Called via:
C           call rsach(file, nerr)
C
C        With:
C           file - Character string giving SAC file name
C           nerr - Integer value (returned).  If nonzero, error reading file 
C              header.

      subroutine rsach(file,nerr)
      character file*(*)
      integer nerr

      logical saciou

      include 'rsach.sachdr.com'

      nvhdr = 0
      call newhdr
      if (nvhdr .eq. 0) then
         write(0,*) '**SAC header format changed--giving up.'
	 nerr = -1
	 return
      endif

C     We assume that if the header got set to something, it is OK.
      if (.not. saciou(iunit)) then
	 nerr = -2
	 write(0,*) '**Can''t find unit to read SAC file header!'
	 return
      endif

C     Read file and get header contents.
      open(iunit,file=file,status='old',access='direct',
     +   recl=lencom*4,iostat=nerr)
      if (nerr .eq. 0) read(iunit,rec=1,iostat=nerr) sachdr
      if (nerr .ne. 0) then
	 write(0,*) '**Can''t read SAC file ',file(1:index(file,' '))
	 return
      endif
      close(iunit)

C     The header is actually split among two common areas.  Move the
C        character information from one to the other.
      do 10 i=1,lenchr
	 khdr(i) = kcom(i)
10    continue
      end

      subroutine wsach(file,nerr)
C     WSACH -- Write SAC header from memory.
C
C        Called via:
C           call wsach(file, nerr)
C
C        With:
C           file - Character string giving SAC file name
C           nerr - Integer value (returned).  If nonzero, error writing file 
C              header.
      character file*(*)
      integer nerr

      logical saciou

      include 'rsach.sachdr.com'
      character record*(lencom*4)
      equivalence (sachdr(1),record)

C     Check that header version number is correct.
      if (nvhdr .ne. 6) then
	 nerr = -1
	 write(0,*) '**Version number isn''t right in header.',
     +      '  Initialized?'
	 return
      endif

C     The header is actually split among two common areas.  Move the
C        character information from one to the other.
      do 10 i=1,lenchr
	 kcom(i) = khdr(i)
10    continue

      if (.not. saciou(iunit)) then
	 nerr = -2
	 write(0,*) '**Can''t find unit to read SAC file header!'
	 return
      endif

C     Open file and write header contents.
      open(iunit,file=file,status='old',access='direct',
     +   recl=lencom*4,iostat=nerr)
      if (nerr .eq. 0) write(iunit,rec=1,iostat=nerr) record
      if (nerr .ne. 0) then
	 write(0,*) '**Can''t write SAC file ',file(1:index(file,' '))
	 return
      endif
      close(iunit)
      end

      logical function saciou(iunit)
C     SACIOU  --  Return an unused I/O unit number.  If all in use, return
C        .false.  Otherwise, return .true. and the unit number.
      integer iunit
      logical usable, inuse

C     Look for a FORTRAN logical I/O unit with which to read the header.
      do 10 iu=9,99
	 inquire(unit=iu,exist=usable,opened=inuse)
	 if (usable .and. .not. inuse) go to 19
10    continue
      saciou = .false.
      return

19    continue
      saciou = .true.
      iunit = iu
      end
