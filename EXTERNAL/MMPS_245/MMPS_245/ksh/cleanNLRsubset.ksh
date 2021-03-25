#!/bin/ksh
# Give as an agrument a file containing a space-separated list of subject IDs

#touchFilename=conform.touch                  # Oddly enough, this one is for affine registration to baseline
touchFilenameNLR=longitudinal_nonlinreg.touch
touchFilenameNLRFlow=longitudinal_nonlinregFlow.touch
touchFilenameNLRROIL=longitudinal_nonlinregROI_left_ROIs.touch
touchFilenameNLRROIR=longitudinal_nonlinregROI_right_ROIs.touch

segstatsFilenameNLR=segstats_longitudinal.mat
segstatsFilenameNLRE=segstats_longitudinal_erode2.0.mat
segstatsFilenameNLRFlow=segstats_longitudinalFlow.mat
segstatsFilenameNLRFlowE=segstats_longitudinalFlow_erode2.0.mat
segstatsFilenameNLRROI=segstats_longitudinalROI.mat
segstatsFilenameNLRROIE=segstats_longitudinalROI_erode2.0.mat

parentPath=/space/md2/2/data/MMILDB/ADNI/Containers

echo Using $1
subjList=`cat $1`
#subjList=`cat ~adniproj/MetaData/ADNI/subjList.txt`

for subj in $subjList; do
    
    ignoreSubject=0
    
    skip=1
    dirList=`ls -d $parentPath/FREESURFERRECON_${subj}*`
    for dir in $dirList; do
	if [ $skip = 1 ]; then
	    skip=0
	    #echo skipping baseline
	    continue
	fi
	
	mkdir -p $dir/mri/old
	
	nlrdir=$dir/mri/nonlinreg_to_baseline
	nlrdirF=$dir/mri/nonlinregFlow_to_baseline
	nlrdirR=$dir/mri/nonlinregROI_to_baseline
	
	if [ -d $nlrdir ]; then
	    mv $nlrdir $dir/mri/old/nonlinreg_to_baseline
	fi
	if [ -d $nlrdirF ]; then
	    mv $nlrdirF $dir/mri/old/nonlinregFlow_to_baseline
	fi
	if [ -d $nlrdirR ]; then
	    mv $nlrdirR $dir/mri/old/nonlinregROI_to_baseline
	fi
	
	touchFileNLR=$dir/touch/$touchFilenameNLR
	touchFileNLRFlow=$dir/touch/$touchFilenameNLRFlow
	touchFileNLRROIL=$dir/touch/$touchFilenameNLRROIL
	touchFileNLRROIR=$dir/touch/$touchFilenameNLRROIR
	
        if [ -f $touchFileNLR ]; then
	    rm $touchFileNLR
	fi
        if [ -f $touchFileNLRFlow ]; then
	    rm $touchFileNLRFlow
	fi
        if [ -f $touchFileNLRROIL ]; then
	    rm $touchFileNLRROIL
	fi
        if [ -f $touchFileNLRROIR ]; then
	    rm $touchFileNLRROIR
	fi
	
	segstatsFileNLR=$dir/stats/$segstatsFilenameNLR
	segstatsFileNLRE=$dir/stats/$segstatsFilenameNLRE
	segstatsFileNLRFlow=$dir/stats/$segstatsFilenameNLRFlow
	segstatsFileNLRFlowE=$dir/stats/$segstatsFilenameNLRFlowE
	segstatsFileNLRROI=$dir/stats/$segstatsFilenameNLRROI
	segstatsFileNLRROIE=$dir/stats/$segstatsFilenameNLRROIE

        if [ -f $segstatsFileNLR ]; then
	    rm $segstatsFileNLR
	fi
        if [ -f $segstatsFileNLRE ]; then
	    rm $segstatsFileNLRE
	fi
        if [ -f $segstatsFileNLRFlow ]; then
	    rm $segstatsFileNLRFlow
	fi
        if [ -f $segstatsFileNLRFlowE ]; then
	    rm $segstatsFileNLRFlowE
	fi
        if [ -f $segstatsFileNLRROI ]; then
	    rm $segstatsFileNLRROI
	fi
        if [ -f $segstatsFileNLRROIE ]; then
	    rm $segstatsFileNLRROIE
	fi
	
    done
    
done
