#!/bin/csh

if($1 == "") then 
echo "Usage: script followed by Subject_Name" 
exit 1 
endif

set extension = "_view3D.pdf" 
set title = "Subject # $1" 

set dir = /space/monkeys/1/home/igor/TRIANGLE/NEW_TRI_FILES

$dir/a.out $1 1 1 dummy & 
sleep 1 
import -window dummy v11.tif 
killall a.out 

$dir/a.out $1 1 2 dummy & 
sleep 1
import -window dummy v12.tif 
killall a.out 

$dir/a.out $1 1 3 dummy & 
sleep 1 
import -window dummy v13.tif 
killall a.out 

$dir/a.out $1 2 1 dummy & 
sleep 1 
import -window dummy v21.tif 
killall a.out 

$dir/a.out $1 2 2 dummy & 
sleep 1 
import -window dummy v22.tif 
killall a.out 

$dir/a.out $1 2 3 dummy & 
sleep 1 
import -window dummy v23.tif 
killall a.out 

$dir/a.out $1 3 1 dummy & 
sleep 1 
import -window dummy v31.tif 
killall a.out 

$dir/a.out $1 3 2 dummy & 
sleep 1 
import -window dummy v32.tif 
killall a.out 

$dir/a.out $1 3 3 dummy & 
sleep 1 
import -window dummy v33.tif 
killall a.out 


montage -tile 3x3 -geometry 512x512+30+30 -borderwidth 50 -title "$title"  \
v11.tif v12.tif v13.tif v21.tif v22.tif v23.tif  v31.tif v32.tif v33.tif a.tif

/usr/bin/convert a.tif $1$extension 

rm *tif

acroread $1$extension 


