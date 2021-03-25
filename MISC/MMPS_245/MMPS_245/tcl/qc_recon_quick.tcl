#include <tcl.h>

##############################################################################

# script for quality control (QC) of surface reconstruction

# Created:  01/22/2016 by Don Hagler
# Last Mod: 10/25/2016 by Don Hagler

# Scripting reference available at:
# https://surfer.nmr.mgh.harvard.edu/fswiki/TkMeditGuide/TkMeditReference/TkMeditScripting

##############################################################################

# SetOrientation orientation 
# orientation:
# 0     coronal
# 1     horizontal
# 2     sagittal
SetOrientation 2
SetSlice 110

##
SetVolumeMinMax 0 0 255
SetVolumeMinMax 1 0 255

##  volume brightness contrast 0=main 1=aux brainmask +T1
SetVolumeBrightnessContrast 0 .35 12
SetVolumeBrightnessContrast 1 .35 12

# set the slice number
SetDisplayFlag 5 0
SetDisplayFlag 15 0
SetDisplayFlag 16 0
SetZoomLevel 1.5

RedrawScreen

