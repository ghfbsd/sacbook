      subroutine rsach(file, nerr)
C     RSACH -- Substitute routine for RSACH when your SAC I/O library does not
C              have one.
C
C     Called via:  call rsach(file, nerr)
C
C     Assumes:
C        file - character variable giving file name.
C
C     Returns:
C        nerr - integer variable = 0 if file read OK, not zero otherwise
C
C     Function result:
C        SAC header information read in so that GETxHV and SETxHV routines
C        can be used.
C
C     This implementation relies on the GCC Fortran routine STAT.
C
C     G. Helffrich/U. Bristol
C         7 Dec. 2013

      character file*(*)
      integer nerr
      integer info(13)
      character ktype*8


      call stat(file, info, nerr)
      if (nerr .eq. 0) then
         nwds = info(8)/4 - 158
	 call rsac1(file, data, nlen, beg, dt, 0, nerr)
	 if (nerr .ne. 0) return
	 call getnhv('nvhdr',nvhdr,nerr)
	 nerr = nvhdr-6
	 if (nerr .ne. 0) return
	 call getihv('iftype',ktype,nerr)
	 if (ktype .eq. 'ITIME' .or. ktype .eq. 'IXYZ') then
	    npts = nwds
	 else
	    npts = nwds/2
	 endif
	 call setnhv('npts',npts,nerr)
      endif
      end

      subroutine wsach(file, nerr)
C     WSACH -- Substitute routine for WSACH when your SAC I/O library does not
C              have one.
C
C     Called via:  call wsach(file, nerr)
C
C     Assumes:
C        file - character variable giving file name.
C
C     Returns:
C        nerr - integer variable = 0 if header successfully written OK,
C           not zero otherwise
C
C     Function result:
C        SAC header information written out to file.
C
C     This implementation relies on the GCC Fortran routine GETPID, UNLINK,
C     MALLOC and implementation of Cray pointer extension to Fortran 77.
C     (On g77 and gfortran, this requires compiling with -fcray-pointer)
C
C     G. Helffrich/U. Bristol
C         7 Dec. 2013

      character file*(*)
      integer nerr
      real ypts(*),xpts(*)
      character ktype*8, scr*128
      pointer (pxpts,xpts),(pypts,ypts)


      call getnhv('npts', npts, nerr)
      if (nerr .ne. 0) return
      call getihv('iftype',ktype,nerr)
      if (nerr .ne. 0) return

C     Temporarily preserve updated header data in scratch file
      write(scr,'(a,i8.8,a)') '/tmp/sac',getpid(),'.tmp'
      call setnhv('npts', 0, nerr)
      call wsac0(scr,data,data,nerr)
      if (nerr .ne. 0) return

      pypts = malloc(npts*4)
      if (ktype .eq. 'ITIME' .or. ktype .eq. 'IXYZ') then
	 call rsac1(file, ypts, n, beg, dt, npts, nerr)
	 if (nerr .ne. 0) return
	 call rsac1(scr, data, n, beg, dt, 0, nerr)
	 if (nerr .ne. 0) return
	 call setnhv('npts', npts, nerr)
	 if (nerr .ne. 0) return
	 call wsac0(file, ypts, ypts, nerr)
      else
         pxpts = malloc(npts*4)
	 call rsac2(file, ypts, n, xpts, npts, nerr)
	 if (nerr .ne. 0) return
	 call rsac2(scr, data, n, data, 0, nerr)
	 if (nerr .ne. 0) return
	 call setnhv('npts', npts, nerr)
	 if (nerr .ne. 0) return
	 call wsac0(file, xpts, ypts, nerr)
      endif
      call unlink(scr)
      end
