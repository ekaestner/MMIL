#!/bin/csh

if($1 == "") then 
echo "Usage: script followed by Subject_Name" 
exit 1 
endif

set mainCodeDirectory = /space/monkeys/1/home/igor/TRIANGLE 

set mriDirectory = /space/monkeys/1/home/igor/PROJECTS 
set mriFile = $1.mgh 

echo $mriDirectory/$mriFile > dummy 

set triInputDirectory = /home/igor/TRIANGLE/NEW_TRI_FILES/ 

echo $triInputDirectory >> dummy 

set meshNumber = 5 

echo $meshNumber >> dummy 

# INNER SKULL 

set nMove = 300 
set nRest = 30 
set nRelax = 3 
set coefTangential = 0.03 
set coefNormal = 0.03 
set coefMRI = 0.3 
set directionMRI = 1 
set coefRepelling = 0.
set intensityThreshould = 135. 

echo $nMove >> dummy 
echo $nRest >> dummy
echo $nRelax >> dummy
echo $coefTangential >> dummy 
echo $coefNormal >> dummy 
echo $coefMRI >> dummy 
echo $directionMRI >> dummy 
echo $coefRepelling >> dummy 
echo $intensityThreshould >> dummy 


# OUTER SKULL 

set nMove = 100 
set nRest = 30 
set nRelax = 3 
set coefTangential = 0.05 
set coefNormal = 0.1 
set coefMRI = 0.2
set directionMRI = -1 
set coefRepelling = 0.5
set intensityThreshould = 90.

echo $nMove >> dummy 
echo $nRest >> dummy
echo $nRelax >> dummy
echo $coefTangential >> dummy 
echo $coefNormal >> dummy 
echo $coefMRI >> dummy 
echo $directionMRI >> dummy 
echo $coefRepelling >> dummy 
echo $intensityThreshould >> dummy 

# OUTER SCALP  

set nMove = 300 
set nRest = 30 
set nRelax = 1 
set coefTangential = 0.01 
set coefNormal = 0.01
set coefMRI = 0.3
set directionMRI = 1 
set coefRepelling = 0.5
set intensityThreshould = 45.

echo $nMove >> dummy 
echo $nRest >> dummy
echo $nRelax >> dummy
echo $coefTangential >> dummy 
echo $coefNormal >> dummy 
echo $coefMRI >> dummy 
echo $directionMRI >> dummy 
echo $coefRepelling >> dummy 
echo $intensityThreshould >> dummy 

$mainCodeDirectory/a.out < dummy 

mv dummy $1.par

