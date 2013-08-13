#! /bin/sh
# Shell script to build SAC macro to execute that will label a plot with
# travel time information.  Input is from standard input, line-by-line:
#   1 - event depth
#   2 - event distance (degrees)
#   3 - O value from SAC file header
#   4 - start and end time for plot (seconds relative to file zero)
#   5 - Event ID string
#   6 - Phase id (or "none")
#   7 - Time window around phase (seconds)
#   8 - Name of travel time tables (or default)
#   9 - Explicit phase names to mark or "all"
#
# by G. Helffrich/U. Bristol
#   12 Oct. 2006
# Copyright (c) 2013 by G. Helffrich.
# All rights reserved.
#
# Redistribution and use in source form, with or without
# modification, is permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer:
# THERE IS NO WARRANTY EXPRESSED OR IMPLIED, NOR LIABILITY FOR DAMAGES ASSUMED,
# IN PROVIDING THIS SOFTWARE.

file=./ttsac.pcf ttfile=/tmp/ttimes.$$
# Default travel time calculator:
#   1 - Crotwell et al. taup; 0 - Buland and Kennett ttimes
defttcalc=0

read depth; read dist; read otime; read seism_start seism_end
read infoline; read phase; read window; read tables; read phases
info=`echo $infoline | awk -F_ '{printf("Event:%s  _Sta:%s  _Dist:%.1f  \
_Az:%.1f  _Baz:%.1f  _%.3fN  _%.3fE  %.1fkm\n",$1,$2,$3,$4,$5,$6,$7,$8)}'`
for arg in ${tables} ; do
   case $arg in
      taup) taup=1 ;;
      ttimes) taup=0 ;;
      default) [ -f .ttsac ] && . .ttsac
         [ "_${tables}" = "_default" ] && taup=1 tables=ak135 ;;
      *) tables=$arg
         echo "taup=${taup} tables=${tables}" > .ttsac ;;
   esac
done
[ "_${tables}" = "_default" ] && taup=${defttcalc} tables=ak135

# echo "depth $depth" >> /dev/tty
# echo "dist $dist" >> /dev/tty
# echo "otime $otime" >> /dev/tty
# echo "seism_start $seism_start seism_end $seism_end" >> /dev/tty
# echo "infoline $infoline" >> /dev/tty
# echo "phase $phase" >> /dev/tty
# echo "window $window" >> /dev/tty
# echo "tables $tables" >> /dev/tty
case $taup in
0) which ttimes >/dev/null 2>&1 || taup=-1
   ;;
1) which taup_time >/dev/null 2>&1 || taup=-1
   ;;
*) ;;
esac
if [ $taup -lt 0 ] ; then
   echo "**No travel time tables accessible" > /dev/tty
   > $ttfile
elif [ $taup -eq 0 ] ; then
   eval ttimes -model $tables << EOF $pipe > $ttfile
all

$depth
$dist
EOF
else
   taup_time -mod ${tables} -deg ${dist} -h ${depth} \
      -ph ttall,PKJKP,SKKS,SKKKS,SKKKKS |
   awk 'BEGIN{olddel=-1}
      / = /{del=$1; if (olddel != del) {n=0; olddel=del}; n+=1;
         if ($1==$6) sl=$5; else sl=-1*$5
         if (n==1) printf"%7.2f   1  %-10s %7.2f   0  00.00 %10.4f %10.2e %10.2e\n",$1,$3,$4,sl,0.0,0.0
         else printf"        %3d  %-10s %7.2f   0  00.00 %10.4f %10.2e %10.2e\n",n,$3,$4,sl,0.0,0.0
   }' |
   eval "cat $pipe > $ttfile"
fi
seism_origin_diff=$otime
cat /dev/null > $file
#Forces software text
#echo '[qsf3hlvbss]'           >> $file
echo '[f3hcvbss]'             >> $file
echo 't 0.50 0.94'            >> $file
echo $info                    >> $file
echo '[f3stvb]'               >> $file
echo 't 0.15 0.08'            >> $file
echo "$tables tables"         >> $file
echo '[f3stvb]'               >> $file
echo $seism_start $seism_end $seism_origin_diff "$phase" $window "$phases" |
cat - $ttfile |
awk '
#	This awkscript maps the start and end times (sec_beg and sec_end) in 
#	seconds on the seismogram to plotting device units. Because the travel 
#	times are in seconds relative to the origin time, it also must know 
#	the difference (sec_diff) between the origin time and the start time 
#	of the seismogram.
#	The magic device settings (mach_[x,y][max,min]) are for the SAC 
#	devices "sun" and "sgf".
#	Methodology is to run through the ttimes output and save all phase
#	names and arrival times.  Then those that sit within the time window
#	are output.  Comments are placed in the .pcf file for all phases for
#	error checking purposes.
BEGIN { mach_xmin = 0.1; mach_xmax = 0.9; mach_ymin = 0.1; mach_ymax = 0.87
        mach_offset = 0.024
	np = 0 }
NR == 1 { sec_beg=0.0+$1; sec_end=0.0+$2; sec_diff=0.0+$3
	  phase_match=$4; time_window=0.0+$5
	  phlist=$6; for (i=7;i<=NF;i++) phlist=phlist " " $(i)
	  phase_min_time=864000.0; phase_max_time=0.0; n_matched = 0
#	  If phase name has no branch suffix, let it match any branch.
	  if (phase_match !~ /^.*((ab)|(bc)|(df)|(ac)|[gnb])$/)
	     phase_match = phase_match "((ab)|(bc)|(df)|(ac)|[gnb])?$"
	  if (phlist=="all") {nphl=1; phl[1]=".*"} else nphl=split(phlist,phl)
}
($0 ~ /.*\..*\..*\..*\..*\./) || (NR == 1) {
   if (NR==1) {
      phase="Origin"; tt=0.0
      printf("* Matching %s\n",phase_match)
   } else
      if (NF==9) {phase=$3; tt=0.0+$4} else {phase=$2 ; tt=0.0+$3}
   pmatch=0; for (i=1;i<=nphl && !pmatch;i++) if (phase ~ phl[i]) pmatch=1
   if (!pmatch) next
   np += 1; pid[np] = phase; ptt[np] = tt; ph[phase] += 1
   if ( phase ~ phase_match ) {
      printf("* Matched %s\n",phase)
      n_matched += 1
      if ((tt+sec_diff >= sec_beg) && (tt+sec_diff <= sec_end)) {
#        Only consider as arrival if it arrives in the time window of interest.
	 if (tt < phase_min_time) phase_min_time = tt
	 if (tt > phase_max_time) phase_max_time = tt
      }
   }
}
END {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	If defining a window around some phase arrival, set that bracket.
#	Otherwise use the explicit bounds given.  Run through the collected
#	phases and emit only those that are within the time bounds.
   if ( phase_match != "none") {
      if (n_matched > 0) {
#        If phase arrivals are separated by more than twice the seismogram
#           length, cut off the one that is last.  This keeps things like
#           major arc SKKS arrivals from unreasonably lengthening the display.
#        Do not make window wider than seismogram, however.  Looks rather silly.
	 if (1.5*(sec_end-sec_beg) < (phase_max_time-phase_min_time)) {
	    phase_max_time = phase_min_time
	 }
	 twlo = phase_min_time - time_window + sec_diff
	 twhi = phase_max_time + time_window + sec_diff
	 if (twlo >= sec_beg) sec_beg = twlo
	 if (twhi <= sec_end) sec_end = twhi
      }
      printf "xlim %.0f %.0f\nsetbb ttsacxlim \"%.0f %.0f\"\n",
	 sec_beg,sec_end,sec_beg,sec_end > "./ttsac.xlim" 
   }
   for (n=1; n<=np; n++) {
      printf("* %s %f\n",pid[n],ptt[n])
      mach_time = mach_xmin + ( (mach_xmax - mach_xmin) * ( (ptt[n] + sec_diff - sec_beg) / (sec_end - sec_beg) ) )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Check to see how close two arrivals may be -- if they are closer than 
#       0.02 units units in the x direction, then the label will be "offset" 
#       by 0.02 units in the y direction.
      if ( (mach_time - prev_mach_time) < 0.02) 
	 offset = prev_offset + mach_offset
      else 
	 offset = 0.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Check to see that the arrival is in the window covered by the seismogram
#	-- if not, then do not include in the overlay.
      if ( (mach_time <= mach_xmax) && (mach_time >= mach_xmin) )
	 printf("t %.3f %.3f\n%s\no %.3f %.3f\n[l2]\nl %.3f %.3f\n[l1]\n",
	 mach_time,(mach_ymax+offset),pid[n],mach_time,(mach_ymax+offset),mach_time,(mach_ymin+0.05))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      prev_mach_time=mach_time
      if (offset >= 0.04) prev_offset = -0.02; else prev_offset = offset
   }
   print  "Q"
}
' >> $file
echo 'Q' >> $file
/bin/rm -f $ttfile
exit
