#!/bin/ksh
# Give as an agrument a file containing a space-separated list of subject IDs

#touchFilename=conform.touch                  # Oddly enough, this one is for affine registration to baseline
#touchFilename=longitudinal_nonlinreg.touch
touchFilename=longitudinal_nonlinregROI_left_ROIs.touch
#touchFilename=longitudinal_nonlinregROI_right_ROIs.touch

parentPath=/space/md2/2/data/MMILDB/ADNI/Containers

echo Using $1
subjList=`cat $1`
#subjList=`cat ~adniproj/MetaData/ADNI/subjList.txt`



for subj in $subjList; do
    skip=1
    dirList=`ls -d $parentPath/FREESURFERRECON_${subj}*`
    for dir in $dirList; do
	if [ $skip = 1 ]; then
	    skip=0
	    #echo skipping baseline
	    continue
	fi
	touchFile=$dir/touch/$touchFilename
        if [ ! -f $touchFile ]; then
	    echo Could NOT find ${touchFile}
	    tail -1 $dir/scripts/fs_recon.log
#        if [ -f $touchFile ]; then
#	    echo Found ${touchFile}
#	    ls -l $touchFile
#	else
#	    echo Could NOT find ${touchFile}
        fi
    done
done
