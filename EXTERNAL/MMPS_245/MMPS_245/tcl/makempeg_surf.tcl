#!/usr/local/bin/tclsh

set output "surf"
set savepath [subst "$env(SUBJECTS_DIR)/[GetSubjectName 0]/mpeg"]

SetDisplayFlag 5 0 # toggle off original surface
SetDisplayFlag 3 0 # toggle off cursor
SetVolumeBrightnessContrast 0 .53 12
SetVolumeBrightnessContrast 1 .53 12

SetOrientation 0 # do coronal slices
for {set i 0} {$i < 255} {incr i} {
   SetSlice $i
   set fileno [format %03d $i]
   RedrawScreen
   SaveTIFF $savepath/$output\_COR_$fileno.tif
}

SetOrientation 1 # do horizontal slices
for {set i 0} {$i < 255} {incr i} {
   SetSlice $i
   set fileno [format %03d $i]
   RedrawScreen
   SaveTIFF $savepath/$output\_HOR_$fileno.tif  
}

SetOrientation 2 # do sagittal slices
for {set i 0} {$i < 255} {incr i} {
   SetSlice $i
   set fileno [format %03d $i]
   RedrawScreen
   SaveTIFF $savepath/$output\_SAG_$fileno.tif    
}

QuitMedit
