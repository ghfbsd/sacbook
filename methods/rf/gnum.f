      logical function gnum(string,value,errmsg)
C     gnum -- Get a number from a string, and generate error message
C             if bad.
      character string*(*),errmsg*(*)

      ios = -1
      if (string .ne. ' ') read(string,*,iostat=ios) value
      if (ios .ne. 0) then
	 if (len(errmsg).gt.0) write(0,*) '**Invalid ',errmsg,'.'
	 gnum = .true.
      else
	 gnum = .false.
      endif
      end
