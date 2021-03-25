#!/bin/tcsh -f

################################################################################
#   Default MMIL SetUpDCMTK.csh file
#   Created 2011-01-07 tcooper
#   Master Copy Stored at /usr/pubsw/scripts/SetUpDCMTK.csh
#
#   Currently supported versions of DCMTK on workstations and all clusters are...
#
#   3.5.3    3.5.4    3.6.0
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
### Remove previous DCMTK Setup
if ( `env | grep -c DCMTK_VER` == 1 ) then
i# echo "Removing previous DCMTK Setup..."
  unsetenv DCMTK_VER
  unsetenv DCMTK_HOME
  set path=`echo $path | tr ' ' '\n' | awk '$0 !~ "/DCMTK/"' | paste -sd' '`
endif

################################################################################
### Supply the desired DCMTK version string as an input. Use 'default'
### DCMTK version of 360 if not supplied.
setenv DCMTK_VER 360


################################################################################
### Parse input
if ($#argv == 0 ) then
# echo "No arguments supplied... using default DCMTK Version (RH4-x86_64-R${DCMTK_VER})"
else
    setenv DCMTK_HOME ${PUBSW}/packages/dcmtk/$argv[1]
  endif
  if (! -e ${DCMTK_HOME}) then
    # echo "Could not find DCMTK version specified... using default DCMTK Version (${DCMTK_VER})"
    unsetenv DCMTK_HOME
  else
    setenv DCMTK_VER $argv[1]
  endif
endif

setenv DCMTK_HOME ${PUBSW}/packages/dcmtk/${DCMTK_VER}
set addpathlist=($DCMTK_HOME/bin)

################################################################################
### Set DCMTK MAN Pages
switch ($DCMTK_VER)
  case [3].[6].[0-9]:
    setenv MANPATH "${DCMTK_HOME}/share/man:$MANPATH"
    breaksw

  case [3].[5].[0-9]:
    setenv MANPATH "${DCMTK_HOME}/man:$MANPATH"
    breaksw
endsw

foreach dir ( $addpathlist )
  if ("$path" !~ *$dir*) set path = ($dir $path)
end
unset addpathlist dir

################################################################################
### Clean Path
set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`
