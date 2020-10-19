#include <tcl.h>

##############################################################################

# script for quality control (QC) of surface reconstruction

# Created:  01/22/2016 by Don Hagler
# Last Mod: 02/22/2016 by Don Hagler

# Scripting reference available at:
# https://surfer.nmr.mgh.harvard.edu/fswiki/TkMeditGuide/TkMeditReference/TkMeditScripting

##############################################################################

# SetOrientation orientation 
# orientation:
# 0     coronal
# 1     horizontal
# 2     sagittal
SetOrientation 2
SetSlice 120

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
# SetZoomLevel level 
SetZoomLevel 1.5

# scan through brain
SetCursor 0 110 118 194
SetOrientation 0
for { set slice 45 } { $slice < 190 } { incr slice } {
    SetSlice $slice
    RedrawScreen
    after 15
}

SetCursor 0 140 118 194
SetOrientation 0
for { set slice 45 } { $slice < 190 } { incr slice } {
    SetSlice $slice
    RedrawScreen
    after 20
}

SetCursor 0 118 118 85
SetOrientation 2
for { set slice 110 } { $slice < 145 } { incr slice } {
    SetSlice $slice
    RedrawScreen
    after 25
}
for { set slice 110 } { $slice < 145 } { incr slice } {
    SetSlice $slice
    RedrawScreen
    after 35
}

SetCursor 0 118 118 137
SetOrientation 2
for { set slice 55 } { $slice < 185 } { incr slice } {
    SetSlice $slice
    RedrawScreen
    after 20
}

SetCursor 0 128 118 130
SetOrientation 1
SetSlice 95
RedrawScreen
after 500
SetOrientation 0
SetSlice 108
RedrawScreen
after 500
SetOrientation 2
SetSlice 150
RedrawScreen

SetSlice 110
RedrawScreen

