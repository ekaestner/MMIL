#!/bin/bash

#
# tractoview.sh
#   wrapper for tractoview_bin
#
# created:  04/16/12 by Trevor Cooper
# last mod: 05/17/12 by Don Hagler
# 

appname=`basename $0 | sed s,\.sh$,,`
#echo "  appname=$appname"

dirname=`dirname $0`
#echo "  dirname=$dirname"

tmp="${dirname#?}"
#echo "  tmp=$tmp"

if [ "${dirname%$tmp}" != "/" ]; then
  dirname=$PWD/$dirname
#  echo "  dirname=$dirname"
fi

#echo "  LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
if [ -n "${LD_LIBRARY_PATH:+1}" ] ; then
  LD_LIBRARY_PATH_SAV=$LD_LIBRARY_PATH
  LD_LIBRARY_PATH="${MMPS_LIB}/qt-3.3.3:$LD_LIBRARY_PATH"
else
  LD_LIBRARY_PATH_SAV=""
  LD_LIBRARY_PATH="${MMPS_LIB}/qt-3.3.3"
fi
#cat ${MMPS_LIB}/qt-3.3.3/README
export LD_LIBRARY_PATH
#echo "  LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
#echo "  LDD Result..."
ldd ${dirname}/${appname}_bin > /dev/null

#echo "  Executing... $dirname/$appname $@"
${dirname}/${appname}_bin "$@"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH_SAV;
export LD_LIBRARY_PATH
#echo "  LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
