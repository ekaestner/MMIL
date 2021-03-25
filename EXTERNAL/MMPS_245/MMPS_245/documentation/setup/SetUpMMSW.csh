#!/bin/tcsh -f

################################################################################
#   Default MMIL SetUpMMSW.csh file
#   Created 2010-02-22 tcooper
#   Master Copy Stored at /usr/pubsw/scripts/SetUpMMSW.csh
#
################################################################################

################################################################################
### Default pubsw location
###
### NOTE: This definition places the 'determination' of 'where' /usr/pubsw IS
###       into the hands of the machine maintainter. /usr/pubsw can be linked
###       the the shared storage location (/md7/1/pubsw) or locally.
###
setenv PUBSW /usr/pubsw


################################################################################
### common aliases
alias fixvid '/bin/chmod a+rw /dev/nvidia*'
alias countfiles 'du -a \!:1 | cut -d/ -f2 | sort | uniq -c | sort -nr'


################################################################################
### freesurfer aliases
alias sete  'setenv SUBJECTS_DIR ${PWD}'
alias tkwma 'tkmedit \!:1 wm.mgz lh.white -aux T1.mgz -aux-surface rh.white -segmentation mri/aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt'
alias tkwm  'tkmedit \!:1 wm.mgz lh.white -aux T1.mgz -aux-surface rh.white'


################################################################################
### Setup standard software environment.
### Override by setting $PUBSH to point somewhere else or specify path to specific setup file.
source ${PUBSW}/bin/SetUpMMPS.csh 100
echo "MMPS Version:       ${MMPSVER}"
source ${PUBSW}/bin/SetUpFreeSurfer_new.csh 402
echo "FreeSurfer Version: ${FREESURFER_VER}"
source ${PUBSW}/bin/SetUpMatlab_new.csh R2009b
echo "MATLAB Version:     ${MATLABVER}"
source ${PUBSW}/bin/SetUpAFNI_new.csh
source ${PUBSW}/bin/SetUpFSL_new.csh 4.0.0_32
echo "FSL Version:        ${FSLVER}"


################################################################################
### Clean path
set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`
