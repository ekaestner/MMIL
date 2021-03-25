#!/bin/tcsh -f

############################################################################
#   Default MMIL SetUpFreeSurfer.csh file
#   Created 2010-02-08 tcooper
#   Master Copy Stored at /usr/pubsw/scripts/SetUpFreeSurfer.csh
############################################################################

### Setup the FREESURFER version to run.
set cpu=`uname -m`
set release=`cat /etc/redhat-release`

setenv TCL_LIBRARY /usr/lib
setenv TK_LIBRARY /usr/lib
setenv LANG C
setenv NO_FSFAST 1
setenv NO_FSL 1
setenv FS_FREESURFERENV_NO_OUTPUT 1

if ($#argv == 0 ) then
#  echo "No arguments supplied... using default FreeSurfer Version"
  setenv FREESURFER_VER 450
else
  unsetenv FREESURFER_VER
  unsetenv FREESURFER_HOME
  setenv FREESURFER_VER $argv[1]
endif

if ($cpu == "x86_64") then
  setenv FREESURFER_HOME ${PUBSW}/packages/freesurfer/RH4-x86_64-R${FREESURFER_VER}
else  
  setenv FREESURFER_HOME ${PUBSW}/packages/freesurfer/RH4-x86_32-R${FREESURFER_VER}
endif
set addpathlist=($FREESURFER_HOME/bin)


############################################################################
### Source the appropriate FreeSurferEnv*.csh file based on $FREESURFER_VER
switch ($FREESURFER_VER)
  case [5][0-9][0-9]:
    source ${PUBSH}/bin/FreeSurferEnv_5.csh
    breaksw

  case [4][0-9][0-9]:
    source ${PUBSH}/bin/FreeSurferEnv_400.csh
    breaksw

  case [3][0-9][0-9]:
    source ${PUBSH}/bin/FreeSurferEnv_3.csh
    breaksw
endsw

foreach dir ( $addpathlist )
  if ("$path" !~ *$dir*) set path = ($dir $path)
end
unset addpathlist dir


############################################################################
### Default FreeSurfer aliases
alias sete   'setenv SUBJECTS_DIR ${PWD}'
alias tkwm   'tkmedit \!:1 wm.mgz lh.white -aux T1.mgz -aux-surface rh.white'

switch ($FREESURFER_VER)
  case [4-5][0-9][0-9]:
    alias tkwma  'tkmedit \!:1 wm.mgz lh.white -aux T1.mgz -aux-surface rh.white -segmentation aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt'
    breaksw

  case [3][0-9][0-9]:
    alias tkwma 'tkmedit \!:1 wm.mgz lh.white -aux T1.mgz -aux-surface rh.white -segmentation mri/aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt'
    breaksw
endsw
