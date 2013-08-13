#!/bin/sh
#  Helper shell script to splitstack macro.  Has two modes of operation:
#  1) annotation of joint uncertainty plot with file name, S/N and back-az
#     args:  annotate phi dt
#  2) jackknife of sequence of input files to get a jackknife estimate from
#     the individual stacks.
#     args:  jackknife filename outputfile
#
#  G. Helffrich/25 March 2005
#
# Copyright (c) 2013 by G. Helffrich.
# All rights reserved.
#
# Redistribution and use in source form, with or without
# modification, is permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer:
# THERE IS NO WARRANTY EXPRESSED OR IMPLIED, NOR LIABILITY FOR DAMAGES ASSUMED,
# IN PROVIDING THIS SOFTWARE.
case $1 in
  ann*)
     shift; phi=${1:-0.0} dt=${2:-0.0}
     awk 'BEGIN{ox=0.125;oy=0.8;rad=0.1;degrad=3.14159/180.
        fac=(0.75-0.25)/(0.9-0.1)*(0.80-0.05)/(0.74-0.05)
        printf "[vbhc]\nt %f %f\nBack-azimuths\n",ox,oy+1.3*rad
     }
     {ang=(90.0-$3)*degrad; printf "o %f %f\nl %f %f\n",ox,oy,ox+rad*cos(ang)*fac,oy+rad*sin(ang)}
     END{
        ang=90.0-('"${phi}"'); dt='"${dt}"' ; rad*=1.1
        if (ang > 90.0) c="r"; else c="l" ; ang*=degrad
        printf "[h%svc]\nt %f %f\nf\n",c,ox+rad*cos(ang)*fac,oy+rad*sin(ang)
        printf "[ln2]\no %f %f\nl %f %f\n",ox-rad*cos(ang)*fac,oy-rad*sin(ang),ox+rad*cos(ang)*fac,oy+rad*sin(ang)
        ang+=3.14159/2;
        printf "o %f %f\nl %f %f\n",ox-rad*cos(ang)*fac,oy-rad*sin(ang),ox+rad*cos(ang)*fac,oy+rad*sin(ang)
     }'
     ;;
  jack*)
     temp=/tmp/temp$$ file=$2 ofile=$3 n=0
     cat $file | while read var rest; do
	if expr $var : '[^*#]' >/dev/null ; then
	   fgrep -v "${var}" $file > $temp; n=`expr $n + 1`
	   splstacksac $temp ${temp}.j${n}
	   echo  ${temp}.j${n} >> $temp.jack
	fi
     done
     splstacksac -jack $temp.jack $ofile
     /bin/rm -f $temp.j*
     ;;
  *) echo "$0:  Invalid operation mode." ; exit 1
     ;;
esac
