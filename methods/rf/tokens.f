      subroutine tokens(line,nmax,n,token)
C     tokens -- Split input line into blank-separated tokens and return
C               each.
C
C     George Helffrich/U. Bristol
C        Copyright 1992-2013.
      character line*(*), token(nmax)*(*)

C     Find first nonblank
      n = 0
      k = 1
5     continue
	 do 10 i=k,len(line)
	    if (line(i:i) .ne. ' ') go to 20
10       continue
	 return

20       continue
	 do 25 k=i,len(line)
	    if (line(k:k) .eq. ' ') go to 30
25       continue
         k = len(line)+1
30       continue
	 if (n .lt. nmax) then
	    n = n+1
	    token(n) = line(i:k-1)
	 endif
	 i = k
      if (i .lt. len(line)) go to 5
      end

      logical function getnbs(str,i,nbs)
C     getnbs -- Get ith nonblank string from string str, return it in nbs.
C               If found, return true.  If not there, return false.

      character str*(*), nbs*(*)

      ix = 1
      do 10 is=1,i
C        Find next nonblank.
	 do 15 j=ix,len(str)
	    if (str(j:j) .ne. ' ') go to 20
15       continue
         getnbs = .false.
	 return

20       continue
	 ibeg = ix
	 icnt = index(str(ix:),' ')
	 if (icnt .eq. 0) icnt = len(str)-ix+1
	 ix = ix + icnt
10    continue
      nbs = str(ibeg:ibeg+icnt-1)
      getnbs = .true.
      end
