#!/bin/tcsh
####################################################################
# FSL first segmentation for hippocampal volume
#
# usage: firsthippo.csh <fname> <outdir> <outstem>
#
# Required Arguments:
#   fname: full path file fname of T1 weighted image
#
# Optional Arguments:
#   outdir: output directory
#    {default = pwd/firsthippo}
#   outstem: output file stem
#    {default = subj}
#
# created: 11/01/11 by Lars T. Wetlye  l.t.westlye@psykologi.uio.no
# last mod 11/02/11 by Don Hagler
#
####################################################################

# parse input
if ($#argv == 0) then
  echo " "
  echo "USAGE: firsthippo.csh <fname> <outdir> <outstem>"
  echo " "
  echo "Required Arguments:";
  echo "  fname: full path file fname of T1 weighted image"
  echo " "
  echo "Optional Arguments:";
  echo "  outdir: output directory"
  echo "    {default = pwd/firsthippo}"
  echo "  outstem: output file stem"
  echo "    {default = subj}"
  echo " "
  exit 1;
else
  set fname=$1
endif

if ($#argv < 2) then
  set outdir=`pwd`"/firsthippo"
else
  set outdir=$2
endif

if($#argv < 3) then
  set outstem="subj"
else
  set outstem=$3
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

####################################################################

# create output directory
if (! -ed $outdir) then
  mkdir -p $outdir
endif

####################################################################

# convert mgz to nii
mri_convert -i $fname -o ${outdir}/raw.nii.gz 

# reorient to MNI space
cd $outdir
fslreorient2std raw ${outstem}

# flirt data to MNI space
echo "flirt data.."
first_flirt ${outstem} ${outstem}_to_std_sub $flirtopts

# run FIRST
echo "run FIRST"
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
echo "$outstem finished at `date`"
echo "`cat ${outstem}*volume`"
