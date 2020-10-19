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
	
	nlrdir=$dir/mri/nonlinreg_to_baseline
	nlrdirF=$dir/mri/nonlinregFlow_to_baseline
	nlrdirR=$dir/mri/nonlinregROI_to_baseline
	
	if [ ! -d $nlrdir ]; then
	    ignoreSubject=1
#	    echo $nlrdir does not exist
	fi
	if [ ! -d $nlrdirF ]; then
	    ignoreSubject=1
#	    echo $nlrdirF does not exist
	fi
	if [ ! -d $nlrdirR ]; then
	    ignoreSubject=1
#	    echo $nlrdirR does not exist
	fi
	
	touchFileNLR=$dir/touch/$touchFilenameNLR
	touchFileNLRFlow=$dir/touch/$touchFilenameNLRFlow
	touchFileNLRROIL=$dir/touch/$touchFilenameNLRROIL
	touchFileNLRROIR=$dir/touch/$touchFilenameNLRROIR
	
        if [ ! -f $touchFileNLR ]; then
	    ignoreSubject=1
#	    echo Could NOT find ${touchFileNLR}
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $touchFileNLRFlow ]; then
	    ignoreSubject=1
#	    echo Could NOT find ${touchFileNLRFlow}
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $touchFileNLRROIL ]; then
	    ignoreSubject=1
#	    echo Could NOT find ${touchFileNLRROIL}
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $touchFileNLRROIR ]; then
	    ignoreSubject=1
#	    echo Could NOT find ${touchFileNLRROIR}
#	    tail -1 $dir/scripts/fs_recon.log
	fi
#        if [ -f $touchFile ]; then
#	    echo Found ${touchFile}
#	    ls -l $touchFile
#	else
#	    echo Could NOT find ${touchFile}
#        fi
	
	segstatsFileNLR=$dir/stats/$segstatsFilenameNLR
	segstatsFileNLRE=$dir/stats/$segstatsFilenameNLRE
	segstatsFileNLRFlow=$dir/stats/$segstatsFilenameNLRFlow
	segstatsFileNLRFlowE=$dir/stats/$segstatsFilenameNLRFlowE
	segstatsFileNLRROI=$dir/stats/$segstatsFilenameNLRROI
	segstatsFileNLRROIE=$dir/stats/$segstatsFilenameNLRROIE

        if [ ! -f $segstatsFileNLR ]; then
	    ignoreSubject=1
#	    echo Could NOT find $segstatsFileNLR
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $segstatsFileNLRE ]; then
	    ignoreSubject=1
#	    echo Could NOT find $segstatsFileNLRE
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $segstatsFileNLRFlow ]; then
	    ignoreSubject=1
#	    echo Could NOT find $segstatsFileNLRFlow
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $segstatsFileNLRFlowE ]; then
	    ignoreSubject=1
#	    echo Could NOT find $segstatsFileNLRFlowE
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $segstatsFileNLRROI ]; then
	    ignoreSubject=1
#	    echo Could NOT find $segstatsFileNLRROI
#	    tail -1 $dir/scripts/fs_recon.log
	fi
        if [ ! -f $segstatsFileNLRROIE ]; then
	    ignoreSubject=1
#	    echo Could NOT find $segstatsFileNLRROIE
#	    tail -1 $dir/scripts/fs_recon.log
	fi

    done

    if [ "$ignoreSubject" == "1" ]; then
	echo $ignoreSubject
	echo Ignore $subj
    fi

done
