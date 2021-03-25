#!/bin/tcsh
####################################################################
# FSL processing of DTI data
#
# usage: fsl_dti_data.csh <data> <outdir> <forceflag>
#
#
# Required Arguments:
#  data:	full path to corrected 4D dataset "
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
  echo "USAGE: fsl_dti_data.csh <data> <outdir> <forceflag>"
  echo " "
  echo "Required Arguments:"
  echo "   data:	full path to corrected 4D dataset (assumes NEX=1)"
  echo "   outdir:	full path to output directory"         
  echo "   forceflag: overwrite existing output [0|1]"
  exit 1;
else
  set data=$1
  set outdir=$2
  set forceflag=$3
endif

if ($#argv < 3) then
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

# data 
if (! -e $data) then
  echo "ERROR: $data does not exist" >> ${outdir}/errorlog.log
else if (! -e data.nii.gz || $forceflag == 1) then
  mri_convert $data data.nii.gz
#  fslreorient2std data data 
endif

