      subroutine addtim(dyr,ddy,dhh,dmm,dsc,syr,sdy,shh,smm,ssc,secs)
C  ADDTIM  --  Add a time, in seconds, to a time expressed as year, jday,
C              hour, minute, seconds and make it work.
C
C  Assumes:
C     syr, sdy - Beginning year, julian day
C     shh, smm, ssc - Beginning hour, minute, second
C     secs - seconds to add
C  Returns:
C     dyr, ddy - Resulting year, julian day
C     dhh, dmm, dsc - Resulting hour, minute, second
C
C  By G. Helffrich/DTM based on code written by P. Shearer.

      integer  dyr,ddy,dhh,dmm, syr,sdy,shh,smm
      real dsc, ssc, secs

      dyr = syr
      ddy = sdy
      dhh = shh
      dmm = smm
      dsc = ssc+secs
      call CLEANTIME(dyr,ddy,dhh,dmm,dsc)
      end

      subroutine CLEANTIME(yr,jdy,hr,mn,sc)
C   CLEANTIME  --  cleans up oddball yr,jdy,hr,mn,sc times

      integer yr,jdy,hr,mn,dmn,dhr,ddy,dyr

      dmn=int(sc/60.)
      sc=sc-60.*float(dmn)
      if (sc.lt.0.) then
	 sc=sc+60.
	 dmn=dmn-1
      end if
      mn=mn+dmn
c
      dhr=mn/60
      mn=mn-60*dhr
      if (mn.lt.0) then
	 mn=mn+60
	 dhr=dhr-1
      end if
      hr=hr+dhr
c
      ddy=hr/24
      hr=hr-24*ddy
      if (hr.lt.0) then
	 hr=hr+24
	 ddy=ddy-1
      end if
      jdy = jdy+ddy
c
      if (mod(yr,4) .eq. 0 .and. mod(yr,400) .ne. 0) then
	 idyr = 366
      else
	 idyr = 365
      endif
      dyr=jdy/idyr
      jdy = jdy - idyr*dyr
      if (jdy .lt. 0) then
	 jdy = jdy + idyr
	 dyr = dyr-1
      endif
      yr = yr+dyr
      end
