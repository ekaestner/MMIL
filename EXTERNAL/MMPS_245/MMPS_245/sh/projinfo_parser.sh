#!/bin/bash
Usage(){
  echo "#########################################################################"
  echo "Parses MMIL Project info"
  echo "Please note support for double quotes in current version is in beta stage"
  echo "$0 filename_master_projinfo key ProjID worker"
  echo "#########################################################################"
  exit 0
}

[ $# -ne 4 ] && Usage

Projinfo=$1
key=$2
ProjID=$3
worker=$4

[ ! -f $Projinfo ] && echo "Error: master projinfo file doesn't exist" && exit -1

pos=0
worker_pos=0
projid_pos=0
key_pos=0
for keys in `head -n 1 $Projinfo|sed -e "s/,/ /g" -e "s/\"//g"`
do
  ((pos=pos+1))
  keys=`echo $keys|sed -e 's/"//g'`
  case $keys in
    account)
      worker_pos=$pos
      ;;
    ProjID)
      projid_pos=$pos
      ;;
    $key)
      key_pos=$pos
      ;;
  esac
  [ $worker_pos -gt 0 -a $projid_pos -gt 0 -a $key_pos -gt 0 ] && break
done
#[ $worker_pos -eq 0 -o $projid_pos -eq 0 -o $key_pos -eq 0 ] && echo "Error: No record found" && exit -1
[ $worker_pos -eq 0 -o $projid_pos -eq 0 -o $key_pos -eq 0 ] && exit -1
pattern_pre=""
if [ $projid_pos -lt $worker_pos ]
then
  #[ $projid_pos -gt 1 ] && ((pattern_pre_space=projid_pos-1)) && pattern_pre="^\([^,]*,\)\{$((projid_pos-1))\}"
  [ $projid_pos -gt 1 ] && pattern_pre="^\([^,]*,\)\{$((projid_pos-1))\}\"\{0,1\}${ProjID}\"\{0,1\}," || pattern_pre="^\"\{0,1\}${ProjID}\"\{0,1\},"
  [ $((worker_pos-projid_pos)) -gt 1 ] && pattern_pre="$pattern_pre\([^,]*,\)\{$((worker_pos-projid_pos-1))\}"
  pattern="${pattern_pre}\"\{0,1\}${worker}\"\{0,1\},.*"
else
  [ $worker_pos -gt 1 ] && pattern_pre="^\([^,]*,\)\{$((worker_pos-1))\}\"\{0,1\}${worker}\"\{0,1\}," || pattern_pre="^\"\{0,1\}${worker}\"\{0,1\},"
  [ $((projid_pos-worker_pos)) -gt 1 ] && pattern_pre="$pattern_pre\([^,]*,\)\{$((projid_pos-worker_pos-1))\}"
  pattern="${pattern_pre}\"\{0,1\}${ProjID}\"\{0,1\},.*"
fi
record=`grep $pattern $Projinfo`
#[ ${#record} -eq 0 ] && echo "Error: No record found" && exit -1
[ ${#record} -eq 0 ] && exit -1
case `echo $record|wc -l` in
  0)
    #echo "Error: No record found"
    exit -1
    ;;
  1)
    echo $record|cut -d, -f$key_pos|sed -e 's/"//g'
    #eval setenv result `echo $record|cut -d, -f$key_pos`
    #result=`echo $record|cut -d, -f$key_pos`
    #eval result="$result"
    ;;
  *)
    #echo "Error: Multiple records found"
    exit -1
    ;;
esac
