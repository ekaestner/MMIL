#!/bin/tcsh -f

############################################################################
#   Default MMIL SetUpMNE.csh file
#   Created:  03/19/12 by Trevor Cooper
#   Last Mod: 11/07/13 by Don Hagler
#
#   Supported versions of MNE are...
#
#   2.1 & 2.7.0
#
#   NOTE: MUST be sourced AFTER SetUpMatlab.csh in order for MNE/Matlab
#         to be properly integrated.
############################################################################

### Setup the MEG version to run.
set cpu=`uname -m`
set release=`cat /etc/redhat-release`

### Add MNE path
if ($#argv == 0 ) then
#  echo "No arguments supplied... using default MNE Version"
  setenv MNE_VER 2.7.0
else
  unsetenv MNE_VER
  unsetenv MNE_ROOT
  setenv MNE_VER $argv[1]
endif

############################################################################
### Configuration method varies depending on MNE version
switch ($MNE_VER)
  case 2.1:
    setenv MNE_ROOT ${PUBSW}/packages/mnedir/mne-${MNE_VER}
    set addpathlist=($MNE_ROOT/bin/mne)

    foreach dir ( $addpathlist )
      if ("$path" !~ *$dir*) set path = ($dir $path)
    end
    unset addpathlist dir
    breaksw

  case 2.7.0:
    setenv MNE_ROOT ${PUBSW}/packages/mnedir/mne-${MNE_VER}
    setenv MATLAB_ROOT ${MATLABDIR}
    source $MNE_ROOT/bin/mne_setup
    breaksw
endsw
