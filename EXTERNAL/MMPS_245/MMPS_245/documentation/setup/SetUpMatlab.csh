#!/bin/tcsh -f

################################################################################
#   Default MMIL SetUpMatlab.csh file
#   Created 2011-03-30 tcooper
#   Master Copy Stored at /usr/pubsw/scripts/SetUpMatlab.csh
#
#   Currently supported versions of MATLAB  on workstations and all clusters are...
#
#   R2006a_32  R2006b_32 *R2007a  R2008a  R2009a R2010a R2011a R2012a
#   R2006a_64  R2006b_64  R2007b  R2008b *R2009b R2010b R2011b R2012b
#
#   * - installed locally on clusters (others are linked to pubsw NFS share)
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
### Remove previous MATLAB Setup
if ( `env | grep -c MATLABVER` == 1 ) then
#  echo "Removing previous MATLAB Setup..."
  unsetenv MATLABVER
  unsetenv MATLABDIR
  set path=`echo $path | tr ' ' '\n' | awk '$0 !~ "/matlab/"' | paste -sd' '`
endif

################################################################################
### Supply the desired MATLAB version string as an input. Use 'default' MATLAB
### version of R2007a if not supplied.
setenv MATLABVER R2007a

################################################################################
### Parse input
if ($#argv == 0 ) then
#  echo "No arguments supplied... using default MATLAB Version (R2007a)"
else
  if (! -e ${PUBSW}/packages/matlab/$argv[1]) then
#    echo "Could not find MATLAB version specified... using default MATLAB Version (R2007a)"
  else
    setenv MATLABVER $argv[1]
  endif
endif

################################################################################
### Set MATLAB Environment
setenv MATLABDIR ${PUBSW}/packages/matlab/${MATLABVER}

set addpathlist=(${MATLABDIR}/bin)
foreach dir ( $addpathlist )
  if ("$path" !~ *$dir*) set path = ($dir $path)
end
unset addpathlist dir

if( -e /etc/rocks-release ) then
  alias matlab 'matlab -singleCompThread'
endif
alias mat 'matlab -nosplash -nojvm \!*'
alias mat32 'matlab -nosplash -nojvm -glnx86 \!*'
alias matlab32 'matlab -glnx86 \!*'

################################################################################
### Clean Path
set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`
