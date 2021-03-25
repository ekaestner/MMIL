#!/bin/tcsh
####################################################################
# FSL processing of DT measures data
#
# usage: fsl_dti_measlist.csh <fa> <md> <rd> <l1> <b0> <outdir> <forceflag>
#
#
# Required Arguments:
#  fa:		full path to fa vol "
#  md:    full path to md vol (MD)"  
#  rd:		full path to rd vol (TD)"
#  l1:		full path to l1 vol (LD)"
#  b0:		full path to b0 vol (b0)"
#  outdir:	full path to output directory"
#  forceflag: overwrite existing output [0|1]"
#	
# Created:  11/18/11 by Lars T. Westlye  l.t.westlye@psykologi.uio.no
# Last Mod: 11/21/13 by Don Hagler
#
####################################################################

# parse input
if ($#argv == 0) then
  echo " "
  echo "USAGE: fsl_dti_measlist.csh <fa> <md> <rd> <l1> <b0> <outdir> <forceflag>"
  echo " "
  echo "Required Arguments:"
  echo "   fa:		full path to fa vol "
  echo "   md:    full path to md vol (MD)"  
  echo "   rd:		full path to rd vol (TD)"
  echo "   l1:		full path to l1 vol (LD)"
  echo "   b0:		full path to b0 vol (b0)"
  echo "   outdir:	full path to output directory"         
  echo "   forceflag: overwrite existing output [0|1]"
  exit 1;
else
  set fa=$1
  set md=$2
  set rd=$3
  set l1=$4
  set b0=$5
  set outdir=$6
  set forceflag=$7
endif

if ($#argv < 7) then
  echo " "
  echo "ERROR: you need to set all args"
  echo " "
endif

####################################################################

# setup FSL
#if ( `env | grep -c FSLDIR` == 0 ) then
source /usr/pubsw/bin/SetUpFSL.csh 4.1.9_RH5_64
#endif
echo "FSLDIR is ${FSLDIR}" 
setenv FSLOUTPUTTYPE NIFTI_GZ

####################################################################

# set subjID 
set s="subj"
 
# check if structure is ok
if (! -ed $outdir) then
 mkdir $outdir
endif

# convert and reorient
cd $outdir

# b0 
if (! -e $b0) then
  echo "ERROR: $b0 does not exist" >> ${outdir}/errorlog.log
else
  if (! -e data_nodif.nii.gz || ! -e nodif_brain.nii.gz || $forceflag == 1) then
    mri_convert $b0 data_nodif.nii.gz
    fslreorient2std data_nodif data_nodif
    bet data_nodif nodif_brain -f .2 -m
  endif
endif

# FA 
if (! -e $fa) then
  echo "ERROR: $fa does not exist" >> ${outdir}/errorlog.log
else
  if (! -e data_FA.nii.gz || $forceflag == 1) then
    mri_convert $fa data_FA.nii.gz
    fslreorient2std data_FA data_FA 
    fslmaths data_FA -mas nodif_brain_mask data_FA
  endif
endif

# MD
if (! -e $md) then
  echo "ERROR: $md does not exist" >> ${outdir}/errorlog.log
else
  if (! -e data_MD.nii.gz || $forceflag == 1) then
    mri_convert $md data_MD.nii.gz
    fslreorient2std data_MD data_MD 
    fslmaths data_MD -mas nodif_brain_mask data_MD
  endif
endif
  
# RD
if (! -e $rd) then
  echo "ERROR: $rd does not exist" >> ${outdir}/errorlog.log
else
  if (! -e data_RD.nii.gz || $forceflag == 1) then
    mri_convert $rd data_RD.nii.gz
    fslreorient2std data_RD data_RD 
    fslmaths data_RD -mas nodif_brain_mask data_RD
  endif
endif

# L1
if (! -e $l1) then
  echo "ERROR: $l1 does not exist" >> ${outdir}/errorlog.log
else
  if (! -e data_L1.nii.gz || $forceflag == 1) then
    mri_convert $l1 data_L1.nii.gz
    fslreorient2std data_L1 data_L1 
    fslmaths data_L1 -mas nodif_brain_mask data_L1
  endif
endif

