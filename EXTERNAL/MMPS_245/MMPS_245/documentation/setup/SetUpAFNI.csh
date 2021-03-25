#!/bin/tcsh -f

##############################################################################################################
#   Default MMIL SetUpAFNI_new.csh file
#   Created 2010-02-05 tcooper
#   Master Copy Stored at /usr/pubsw/scripts/SetUpAFNI.csh
#
#   AFNI_2007_03_06_0841  AFNI_2007_05_29_1644  AFNI_2008_02_01_1144
#   AFNI_2010_10_19_1028  AFNI_2011_12_21_1014
#
#   AFNI_LATEST
#
#   AFNI_LATEST will always be the most current version of AFNI downloaded from the website.
#   Individual binaries MAY CHANGE within LATEST even though the AFNI_VERSION_LABEL compiled into the code may
#   match the latest 'tagged' version available.
#
#   The AFNI_VERSION_LABEL for the installed version of AFNI in AFNI_LATEST can be found in the file...
#
#      /usr/pubsw/packages/afni/AFNI_LATEST/AFNI_VERSION_LABEL
#
#   If the AFNI environment variable AFNI_VERSION_CHECK is NOT set to NO then AFNI will check the website
#   for the latest released version of AFNI and compare the AFNI_VERSION_LABEL of the running copy of AFNI
#   with that available on the website.
#
#   You can (should?) manually check and record the running/available version of AFNI during automated
#   processing to ensure traceability. This can easily be done with the 'afni_vcheck' binary. For
#   example...
#
#       [~]$ afni_vcheck
#       This program was compiled with the following settings:
#         Version ID   = AFNI_2011_12_21_1014
#       ++ now fetching http://afni.nimh.nih.gov/pub/dist/AFNI.version
#       Latest version listed at AFNI web site:
#         Version ID   = AFNI_2011_12_21_1014
#       
#   Finally, the last version of AFNI to be run by a user (account) should be listed in...
#
#      ~/.afni.vctime
#
##############################################################################################################

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
### Remove previous AFNI Setup
if ( `printenv | grep -c AFNI_VER` >= 1 ) then
  #echo "Removing previous AFNI Setup..."
  unsetenv AFNI_VER
  unsetenv AFNI_HOME
  unsetenv AFNI_NOSPLASH
  unsetenv AFNI_VERSION_CHECK
  unsetenv AFNI_PLUGINPATH

  set path=`echo $path | tr ' ' '\n' | awk '$0 !~ "/afni/"' | paste -sd' '`
endif

################################################################################
### Parse AFNI Version
if ($#argv == 0 ) then
  #echo "No arguments supplied... using LATEST AFNI Version"
  setenv AFNI_VER LATEST
else
  setenv AFNI_VER ${argv[1]}
endif

################################################################################
### Set AFNI path based on CPU architecture
if ($cpu == "x86_64") then
  if ("${AFNI_VER}" == "LATEST") then
    setenv AFNIDIR ${PUBSW}/packages/afni/AFNI_${AFNI_VER}/linux_xorg7_64
  else
    setenv AFNIDIR ${PUBSW}/packages/afni/AFNI_${AFNI_VER}/linux_gcc33_64
  endif
else  
  if ("${AFNI_VER}" == "LATEST") then
    setenv AFNIDIR ${PUBSW}/packages/afni/AFNI_${AFNI_VER}/linux_xorg7
  else
    setenv AFNIDIR ${PUBSW}/packages/afni/AFNI_${AFNI_VER}/linux_gcc32
  endif
endif

set    addpathlist=(${AFNIDIR})
foreach dir ( $addpathlist )
  if ("$path" !~ *$dir*) set path = ($dir $path)
end
unset addpathlist dir

################################################################################
### Setup AFNI Environment
setenv AFNI_NOSPLASH YES
setenv AFNI_VERSION_CHECK YES
setenv AFNI_PLUGINPATH ${AFNIDIR}

################################################################################
### Clean Path
set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`
