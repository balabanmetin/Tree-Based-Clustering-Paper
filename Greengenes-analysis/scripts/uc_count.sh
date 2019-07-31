#!/bin/bash

tot=`cat $1 | wc -l`
#echo $nsn
sn=`cut -f2- $1 | awk '{print NF}' | sort -n | uniq -c | head -n 1 | awk '{print $1}'`
mx=`cut -f2- $1 | awk '{print NF}' | sort -n | tail -n 1`

#mx=`cut -f2- $1 | awk '{print NF}' | sort -n | tail -n

#echo $sn
printf "%d\t%d\t%d\t%d\n" $sn $((tot - sn)) $tot $mx
