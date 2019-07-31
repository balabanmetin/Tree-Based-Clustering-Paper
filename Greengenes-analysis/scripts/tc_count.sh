#!/bin/bash

nsn=`cat $1 | awk '{print $2}' | sort -n | tail -n 1`
#echo $nsn
sn=`cat $1 | awk '{print $2}' | sort -n | uniq -c | head -n 1 | awk '{print $1}'`

#echo $sn
printf "%d\t%d\t%d\n" $sn $nsn $(($sn + nsn))
