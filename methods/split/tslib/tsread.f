c
c   routine TSREAD reads ss from an unformatted file open on unit iunit.
c   The number of points npts is taken from the file.   
c   If an error is encountered, ierr is returned as 1 (normal return 0)
c   
      subroutine tsread(iunit, npts, ss, ierr)
c
      real ss(*)
      ierr = 0
      read(unit=iunit, err=1000) npts, (tdum, ss(i),i = 1, npts)
c
      return 
 1000 ierr = 1
      return 
      end
