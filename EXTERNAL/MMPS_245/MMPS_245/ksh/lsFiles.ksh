#!/bin/ksh
# Give as an agrument a file containing a space-separated list of subject IDs

#touchFile=touch/conform.touch                  # Oddly enough, this one is for affine registration to baseline
#touchFile=touch/longitudinal_nonlinreg.touch
#touchFile=touch/longitudinal_nonlinregROI_left_ROIs.touch
touchFile=touch/longitudinal_nonlinregROI_right_ROIs.touch
#touchFile=touch/longitudinal_nonlinreg.touch
#fullBrainTouchFile=touch/longitudinal_nonlinreg.touch

#imageFile=mri/nonlinregROI_to_baseline/rawavg_reg2baseline_TP2_NonLinReg_5_17_18.mgz    # LEFT
imageFile=mri/nonlinregROI_to_baseline/rawavg_reg2baseline_TP2_NonLinReg_44_53_54.mgz   # RIGHT
#imageFile=mri/nonlinreg_to_baseline/dv.mgz

#imageAltFile=mri/nonlinregROI_to_baseline_NO_FLOW/rawavg_reg2baseline_TP2_NonLinReg_5_17_18.mgz
#imageAltFile=mri/nonlinregROI_to_baseline_NO_FLOW/rawavg_reg2baseline_TP2_NonLinReg_44_53_54.mgz

parentDir=/space/md2/2/data/MMILDB/ADNI/Containers

#subjDirList=`ls -d /space/md2/2/data/MMILDB/ADNI/Containers/FREESURFERRECON*`
subjDirList=`ls $parentDir`

#echo $subjDirList

for subjDir in $subjDirList; do
    
    tFilePath=$parentDir/$subjDir/$touchFile
    iFilePath=$parentDir/$subjDir/$imageFile
    #iAltFilePath=$parentDir/$subjDir/$imageAltFile
    #fullBrainTouchFilePath=$parentDir/$subjDir/$fullBrainTouchFile

    #echo $tFilePath
    #echo $iFilePath
    
    if [ -f $tFilePath ]; then
#    if [ -f $fullBrainTouchFilePath ]; then
	if [ ! -f $iFilePath ]; then
#	if [ -f $iFilePath ]; then
	    #if [ ! $iAltFilePath ]; then
		#echo $subjDir
	    #fi
#	    echo Found "              "$tFilePath
#	    echo Could not find found $iFilePath
	    
#	    echo Found $tFilePath
#	    echo found $iFilePath
	    echo Deleting $tFilePath
	    #ls -l $tFilePath
	    rm $tFilePath
	fi
    fi
done
