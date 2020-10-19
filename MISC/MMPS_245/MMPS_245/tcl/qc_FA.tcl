#include <tcl.h>

##############################################################################

# script for quality control (QC) of DTI FA reg to T1

# Created:  01/22/2016 by Don Hagler
# Last Mod: 02/01/2017 by Don Hagler

# Scripting reference available at:
# https://surfer.nmr.mgh.harvard.edu/fswiki/TkMeditGuide/TkMeditReference/TkMeditScripting

##############################################################################

##
SetVolumeMinMax 0 0 1
SetVolumeMinMax 1 0 255

##  volume brightness contrast 0=main 1=aux
SetVolumeBrightnessContrast 0 .35 12
SetVolumeBrightnessContrast 1 .35 12

# SetZoomLevel level 
SetZoomLevel 1

SetDisplayFlag 1 0
SetDisplayFlag 15 0 # turn off fseg
SetDisplayFlag 16 0 # turn off fseg

# scan through brain in each orientation
foreach orientation {0 1 2} {
    SetOrientation $orientation
    for { set slice 50 } { $slice < 190 } { incr slice } {
    SetSlice $slice
    RedrawScreen
    after 11
  }
}

SetDisplayFlag 15 0 # turn off fseg
SetDisplayFlag 16 0 # turn off fseg
SetOrientation 2
# check corpus registration
for { set slice 120 } { $slice < 135 } { incr slice } {
    SetSlice $slice
    SetDisplayFlag 1 1
    RedrawScreen
    after 150
    SetDisplayFlag 1 0
    RedrawScreen
    after 150
}
