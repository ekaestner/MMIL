#!/usr/local/bin/tclsh

#set savepath [subst "$env(SUBJECTS_DIR)"] #/[GetSubjectName 0]/mpeg
set subjectname [subst "[GetSubjectDir 0]"]
set  subjid "[string range $subjectname 54 70]"
set scantype "[string range $subjectname 90 93]"
#scan subjectname "/space/md10/3/kuperman/Hawaii/Containers/MRIPROCESSED_%s/MPR1_uw.mgh_" subjid



SetDisplayFlag 5 0 # toggle off original surface

SetDisplayFlag 3 0 # toggle off cursor
SetVolumeBrightnessContrast 0 .40 12
#SetVolumeBrightnessContrast 1 .40 12


#SetOrientation 0 # do coronal slices
#for {set i 115} {$i < 135} {incr i} {
#SetSlice $i
#set fileno [format %03d $i]
#RedrawScreen
#SaveTIFF $savepath/$subjectname\_aseg_COR_$fileno.tif
#}
#
#SetOrientation 1 # do horizontal slices
#for {set i 90} {$i < 105} {incr i} {
#SetSlice $i
#set fileno [format %03d $i]
#RedrawScreen
#SaveTIFF $savepath/$subjectname\_aseg_HOR_$fileno.tif
#}

SetOrientation 1 # do sagittal slices
for {set i 130} {$i < 131} {incr i} {
SetSlice $i
set fileno [format %03d $i]
RedrawScreen
SaveTIFF /home/mmilrec/hawaiineoSSs/$subjid\_SAG_$scantype\_$fileno.tif
}

#
#SetOrientation 0 # do coronal slices
#for {set i 0} {$i < 255} {incr i} {
#   SetSlice $i
#   set fileno [format %03d $i]
#   RedrawScreen
#   SaveTIFF $savepath/$output\_COR_$fileno.tif
#}
#
#SetOrientation 1 # do horizontal slices
#for {set i 0} {$i < 255} {incr i} {
#   SetSlice $i
#   set fileno [format %03d $i]
#   RedrawScreen
#   SaveTIFF $savepath/$output\_HOR_$fileno.tif  
#}
#
#SetOrientation 2 # do sagittal slices
#for {set i 0} {$i < 255} {incr i} {
#   SetSlice $i
#   set fileno [format %03d $i]
#   RedrawScreen
#   SaveTIFF $savepath/$output\_SAG_$fileno.tif
#}

QuitMedit
