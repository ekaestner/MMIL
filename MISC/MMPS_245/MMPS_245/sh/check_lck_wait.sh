#!/bin/bash
lckfile=$1
delay=$2
while [ -f $lckfile ]
do
  sleep $2
done
