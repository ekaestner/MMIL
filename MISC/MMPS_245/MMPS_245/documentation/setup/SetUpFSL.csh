#!/bin/tcsh -f

################################################################################
#   Default MMIL SetUpFSL.csh file
#   Created 2010-02-22 tcooper
#   Master Copy Stored at /usr/pubsw/scripts/SetUpFSL.csh
#
#   Currently supported versions of FSL on workstations and all clusters are
#
#   fsl-3.3.3_64  fsl-3.3.8_32RH9  fsl-4.0.0_32  fsl-4.0.0_64  fsl-4.0.4_RH4_64
#
################################################################################

set cpu=`uname -m`
set release=`cat /etc/redhat-release`

################################################################################
### Default pubsw location
###
### NOTE: This definition places the 'determination' of 'where' /usr/pubsw IS
###       into the hands of the machine maintainter. /usr/pubsw can be linked
###       the the shared storage location (/md7/1/pubsw) or locally.
###
setenv PUBSW /usr/pubsw


################################################################################
### Remove previous FSL Setup
if ( `env | grep -c FSLVER` == 1 ) then
#  echo "Removing previous FSL Setup..."
  unsetenv FSLVER
  unsetenv FSLDIR
  set path=`echo $path | tr ' ' '\n' | awk '$0 !~ "/fsl/"' | paste -sd' '`
endif

################################################################################
### Supply the desired FSL version string as an input. Use 'default' FSL
### version of 3.3.8_32RH9 if not supplied.
setenv FSLVER 3.3.8_32RH9

################################################################################
### Parse input
if ($#argv == 0 ) then
#  echo "No arguments supplied... using default FSL Version (fsl-3.3.8_32RH9)"
else
  if (! -e ${PUBSW}/packages/fsl/fsl-$argv[1]) then
#    echo "Could not find FSL version specified... using default FSL Version (fsl-3.3.8_32RH9)"
  else
    setenv FSLVER $argv[1]
  endif
endif

################################################################################
### Set FSL Environment
setenv FSLDIR ${PUBSW}/packages/fsl/fsl-${FSLVER}

set addpathlist=(${FSLDIR}/bin)
foreach dir ( $addpathlist )
  if ("$path" !~ *$dir*) set path = ($dir $path)
end
unset addpathlist dir

setenv TCLLIBPATH ${FSLDIR}/lib
source ${FSLDIR}/etc/fslconf/fsl.csh

################################################################################
### Clean Path
set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`
