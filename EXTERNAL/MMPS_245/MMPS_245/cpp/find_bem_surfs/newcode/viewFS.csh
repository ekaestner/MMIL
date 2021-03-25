#!/bin/csh

if($1 == "") then 
echo "Usage: script followed by Subject_Name" 
exit 1 
endif

# create Tcl script here 

set out = tempor.tcl 

# set subject = "loFA_atl"
set inner_skull = "_inner_skull5.tri"
set outer_skull = "_outer_scalp5.tri"
set outer_scalp = "_outer_skull5.tri"
set width = 2 

echo "#! /usr/pubsw/bin/tixwish" > $out 

echo "set subject $1" >> $out 

echo "set width 2" >> $out  

echo "set inner_skull _inner_skull5.tri" >> $out 
echo "set outer_skull _outer_scalp5.tri" >> $out 
echo "set outer_scalp _outer_skull5.tri" >> $out 

echo "LoadMainSurface 0 ./$1$inner_skull" >> $out 
echo "LoadOriginalSurface 0 ./$1$outer_skull" >> $out 
echo "LoadPialSurface 0 ./$1$outer_scalp" >> $out 

echo "SetSurfaceLineWidth 0 0 $width" >> $out  
echo "SetSurfaceLineWidth 0 1 $width" >> $out 
echo "SetSurfaceLineWidth 0 2 $width" >> $out  

echo "set gOrientation 0" >> $out 
echo "RedrawScreen" >> $out 
echo "SaveTIFF v1.tif" >> $out 

echo "set gOrientation 1" >> $out 
echo "RedrawScreen" >> $out 
echo "SaveTIFF v2.tif" >> $out 

echo "set gOrientation 2" >> $out 
echo "RedrawScreen" >> $out 
echo "SaveTIFF v3.tif" >> $out 

echo "exit" >> $out 

# end of Tcl script 

# run and then remove Tcl script here 

set extension = "_viewFS.pdf" 
set title = "Subject # $1"  

tkmedit -f $1.mgh -tcl $out

rm $out 

# create tiled image in pdf 

set extension = "_viewFS.pdf" 
set title = "Subject # $1"  

montage -tile 1x3 -geometry 512x512+30+30 -borderwidth 50 -title "$title"  \
v1.tif v2.tif v3.tif a.tif

/usr/bin/convert a.tif $1$extension

rm *tif 

acroread $1$extension  
