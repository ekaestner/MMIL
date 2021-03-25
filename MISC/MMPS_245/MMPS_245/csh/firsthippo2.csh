#!/bin/tcsh
####################################################################
# FSL first segmentation for hippocampal volume
#
# usage: firsthippo.csh <fname1> <fname2> <fname3> <outdir> <outstem>
#
# Required Arguments:
#   fname1: full path file fname of T1 weighted image
#   fname2: full path to aseg.mgz or wmparc.mgz
#   fname3: full path to nu.mgz
# Optional Arguments:
#   outdir: output directory
#    {default = pwd/firsthippo}
#   outstem: output file stem
#    {default = subj}
#
# created: 11/01/11 by Lars T. Westlye  l.t.westlye@psykologi.uio.no
#      mod 11/02/11 by Don Hagler
#      mod 11/08/11 by Lars T. Westlye (added calculation of eTIV and TBV)
#
####################################################################

# parse input
if ($#argv == 0) then
  echo " "
  echo "USAGE: firsthippo.csh <fname1> <fname2 <fname3> <outdir> <outstem>"
  echo " "
  echo "Required Arguments:";
  echo "  fname1: full path file fname of T1 weighted image"
  echo "  fname2: full path to aseg.mgz"
  echo "  fname3: full path to nu.mgz"
  echo " "
  echo "Optional Arguments:";
  echo "  outdir: output directory"
  echo "    {default = pwd/firsthippo}"
  echo "  outstem: output file stem"
  echo "    {default = subj}"
  echo " "
  exit 1;
else
  set fname1=$1
  set fname2=$2
  set fname3=$3
endif

if ($#argv < 4) then
  set outdir=`pwd`"/firsthippo"
else
  set outdir=$4
endif

if($#argv < 5) then
  set outstem="subj"
else
  set outstem=$5
endif

####################################################################

# setup FSL
#if ( `env | grep -c FSLDIR` == 0 ) then
source /usr/pubsw/bin/SetUpFSL.csh 4.1.9_RH5_64
#endif
echo "FSLDIR is ${FSLDIR}" 
setenv FSLOUTPUTTYPE NIFTI_GZ

####################################################################

# set some vars
set struc="L_Hipp R_Hipp"
set flirtopts=""
set template=${FSLDIR}/data/standard/MNI152_T1_1mm_brain

####################################################################

# create output directory
if (! -ed $outdir) then
  mkdir -p $outdir
endif

####################################################################

# convert mgz to nii
mri_convert -i $fname1 -o ${outdir}/${outstem}.nii.gz 

# reorient to MNI space
cd $outdir
fslreorient2std ${outstem} ${outstem}

# flirt data to MNI space
echo "flirt data..."
first_flirt ${outstem} ${outstem}_to_std_sub $flirtopts

# run FIRST
echo "run FIRST..."
foreach struct ($struc)
  run_first -i ${outstem} \
    -t ${outstem}_to_std_sub.mat \
    -n 30 -o ${outstem}-${struct}_first \
    -m ${FSLDIR}/data/first/models_336_bin/intref_thal/${struct}.bmv
  first_boundary_corr \
    -s ${outstem}-${struct}_first \
    -o ${outstem}-${struct}_corr \
    -i ${outstem} -b fast 
end

# merge segmentations (3D -> 4D)
echo "merge segmentations..."
fslmerge -t ${outstem}_all_fast_firstseg ${outstem}-L_Hipp_corr ${outstem}-R_Hipp_corr 
fslmerge -t ${outstem}_all_fast_origsegs ${outstem}-L_Hipp_first ${outstem}-R_Hipp_first 

# boundary correction
echo "boundary correction..."
first_mult_bcorr -i ${outstem} -u ${outstem}_all_fast_origsegs -c ${outstem}_all_fast_firstseg -o ${outstem}_all_fast_firstseg 

# count voxels
echo "counting voxels..."
echo "${outstem} Left-Hippocampus `fslstats ${outstem}_all_fast_firstseg -l 16.5 -u 17.5 -V`"  > ${outstem}.volume
echo "${outstem} Right-Hippocampus `fslstats ${outstem}_all_fast_firstseg -l 52.5 -u 53.5 -V`" >> ${outstem}.volume

# housekeeping
imrm ${outstem}-*_first* ${outstem}-*_corr*
echo "......................"
echo "`cat ${outstem}*volume`"

####################################################################

# get brainmask from freesurfer aseg
echo "generating brainmask from freesurfer..." 
cd ${outdir}
mri_binarize --i $fname2 --min 0.2 --dilate 1 --o ${outstem}_grot1.mgz
mri_mask $fname3 ${outstem}_grot1.mgz ${outstem}_masked_grot2.mgz 
mri_convert -i ${outstem}_masked_grot2.mgz -o ${outstem}_brain.nii.gz
fslreorient2std ${outstem}_brain ${outstem}_brain 
rm ${outstem}*mgz

# run flirt and fast on masked nu
echo "running FLIRT to estimate native -> standard xform..." 
flirt -in ${outstem}_brain -ref $template -omat ${outstem}_to_T_brain.mat
echo "running FAST segmentation (GM, WM. CSF)..."
fast ${outstem}_brain

# calculate volumes
echo "calculate volumes..." 
echo "subj_id,eTIV_FLIRT,FASTvol_noCSF, FASTvolGM, FASTvolWM" > ${outstem}_global_size_FSL.csv
set eTIV=`~ltwestlye/programs/scripts/mat2det ${outstem}_to_T_brain.mat | awk '{ print $2 }'`
set volGM=`fslstats ${outstem}_brain_pve_1 -V -M | awk '{ vol = $2 * $3 ; print vol }'`
set volWM=`fslstats ${outstem}_brain_pve_2 -V -M | awk '{ vol = $2 * $3 ; print vol }'`
set voltissue=`expr ${volGM} + ${volWM}`
echo "${outstem},${eTIV},${voltissue},${volGM},${volWM}" >> ${outstem}_global_size_FSL.csv
echo "${outstem} eTIV $eTIV $eTIV" >> ${outstem}.volume
echo "${outstem} TBV $voltissue $voltissue" >> ${outstem}.volume
echo "${outstem} volGM $volGM $volGM" >> ${outstem}.volume
echo "${outstem} volWM $volWM $volWM" >> ${outstem}.volume
echo "......................"
echo "`cat ${outstem}_global_size_FSL.csv`"
echo "......................"
echo "$outstem finished at `date`"
