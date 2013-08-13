C     Set O marker of a file to make the arrival time of a names phase occur
C     at a specified time in the SAC file.
C
C     Command line parameters:
C        -model xxxx - Use xxxx as the travel time model.
C
C     Input is through standard input:
C        file phase time
C     where
C        file - is the file name
C        phase - is the phase name
C        time - is the time in the file for the phase arrival alignment.
C
C     G. Helffrich/U. Bristol 22 Feb. 2002

      parameter (maxph=32, ntok=3)
      character line*256, token(ntok)*128, sacfile*128
      real tt(maxph), dtdd(maxph), dtdh(maxph)
      character idphs(maxph)*8
      equivalence (sacfile,token(1))

C     Parse command line parameters.
      iskip = 0
      do 5 i=1,iargc()
	 if (i .le. iskip) go to 5
	 call getarg(i,line)
	 if (line .eq. '-model') then
	    call getarg(i+1,line)
	    call tpmod(line)
	    iskip = i+1
	 else
	    write(0,*) '**Invalid parameter: ',line(1:lennb(line))
	 endif
5     continue
	    
1000  continue
	 read(5,'(a)',iostat=ios) line
	 if(ios .ne. 0) stop
	 if (line .eq. ' ') go to 1000
	 call tokens(line,ntok,n,token)
	 if (n .lt. 3) then
	    write(0,*) '**SACOMARKER:  Invalid input.'
	    go to 1000
	 endif

	 call rsach(sacfile,nerr)
	 if (nerr .ne. 0) then
	    write(0,*) ' **SACOMARKER:  Unable to read SAC file ',
     +         sacfile(1:lennb(sacfile))
	    stop
	 endif

	 call getfhv('B',fbeg,nerr)
	 if (nerr .ne. 0) go to 999
	 call getfhv('E',fend,nerr)
	 if (nerr .ne. 0) go to 999

	 call getfhv('EVLA',elat,nerr)
	 if (nerr .ne. 0) go to 999
	 call getfhv('EVLO',elon,nerr)
	 if (nerr .ne. 0) go to 999
	 call getfhv('EVDP',edep,nerr)
	 if (nerr .ne. 0) go to 999

	 call getfhv('STLA',slat,nerr)
	 if (nerr .ne. 0) go to 999
	 call getfhv('STLO',slon,nerr)
	 if (nerr .ne. 0) go to 999

	 read(token(3),*,iostat=ios) secs
	 if (ios .ne. 0) then
	    write(0,*) '**SACOMARKER:  Invalid alignment time',
     +         ' given for ',sacfile(1:lennb(sacfile)),'.'
            go to 1000
	 endif

	 call gcd(elat,elon,slat,slon,del,delkm,az,baz)
	 nphs = mtptt('all',del,edep,maxph,idphs,tt,dtdd,dtdh)
	 do i=1,min(maxph,nphs)
	    if (idphs(i) .eq. token(2)) then
	       call setfhv('O',secs-tt(i),nerr)
	       call wsach(sacfile,nerr)
	       if (nerr .ne. 0) then
	          write(0,*) '**SACOMARKER:  Can''t write file header',
     +               ' for ',sacfile(1:lennb(sacfile)),'.'
	       endif
	       go to 1000
	    endif
	 enddo
	       
	 write(0,*) '**SACOMARKER:  No ',token(2)(1:lennb(token(2))),
     +      ' in ',sacfile(1:lennb(sacfile))
         go to 1000

999      continue
	 write(0,*) '**No event info (or it is incomplete) for ',
     +      sacfile(1:lennb(sacfile)),', skipping.'
      go to 1000
      end

      integer function lennb(string)
C     lennb  --  Return non-blank string length.

      character string*(*)

      do 10 i=1,len(string)
	 if (string(i:i) .eq. ' ') go to 15
10    continue
15    continue
      lennb = i-1
      end
