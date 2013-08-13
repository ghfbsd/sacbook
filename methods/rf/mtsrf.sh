#!/bin/sh
svel=$4 model=${5:-ak135}
phs=`echo $3 | awk '{phs=$1    # Strip branch suffix from phase name if any
   if (phs ~ /^.*((ab)|(bc)|(df)|(ac))$/) print substr(phs,1,length(phs)-2)
   else if (phs ~ /^.*[gnb]$/) print substr(phs,1,length(phs)-1)
   else print phs
}'`
cat << EOF |
${phs}

$1
$2
EOF
ttimes -model $model |
awk 'BEGIN{phs="'"${3}"'"; pi=4.0*atan2(1,1); Re=6371; vel='"${svel}"';got=0
   if (phs !~ /^.*((ab)|(bc)|(df)|(ac)|[gnb])$/) {
      phase_match = phs "((ab)|(bc)|(df)|(ac)|[gnb])?$"
   } else {
      phase_match = phs
   }
}
/.*\..*\..*\..*\..*\./{
   if (NF==9) {phase=$3; p=$7+0.0} else {phase=$2; p=$6+0.0}
   if (phase ~ phase_match) {
      p*=180/pi; arg=p*vel/Re; print atan2(arg,sqrt(1-arg*arg))*180/pi
      got=1; exit
   }
}
END{if(!got) {print "**No",phs,"for event found!" > "/dev/tty"; print 0.0}}'
