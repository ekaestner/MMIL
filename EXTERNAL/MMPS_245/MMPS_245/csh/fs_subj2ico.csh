#!/bin/csh -f
# fs_subj2ico.csh -- resample subject to ico
#   optionally copy bem surfs, create dip files
#
# Note: must have run setup for FreeSurfer v3.0.x
#
# Created:  08/29/08 by Don Hagler
# Last Mod: 04/14/10 by Don Hagler
#

set fname_identity = $MMPS_PARMS/identity.xfm

## initialize variables ##
set subj = ()
set subjdir = ()
set outsubj = ()
set outsubjdir = ()
set ico = 4
set forceflag = 0
set bemflag = 0
set dipflag = 0
set hemilist = (lh rh)

if($#argv < 1) then
  goto usage_exit;
  exit 1;
endif

## parse and check params ##
goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

###############################################################################

# replace existing talairach file with identity matrix (for resample to avg)
set fname_talairach = $SUBJECTS_DIR/$subj/mri/transforms/talairach.xfm
set fname_talairach_orig = $fname_talairach'.orig'
if (! -e $fname_talairach_orig) then
  if (-e $fname_talairach) then
    echo cp $fname_talairach $fname_talairach_orig
    cp $fname_talairach $fname_talairach_orig
  endif
  echo cp $fname_identity $fname_talairach
  cp $fname_identity $fname_talairach
endif

# create _ico version of subject resampled to spherical average space
if (! -e $outsubjdir/$outsubj) then
  make_average_subject_v3 --subjects $subj --out $outsubj --ico $ico \
    --sd-out $outsubjdir
else if ($forceflag) then
  make_average_subject_v3 --subjects $subj --out $outsubj --ico $ico \
    --sd-out $outsubjdir --force
endif

# copy bem surfaces for resampled subject
if ($bemflag) then
  set outdir = $outsubjdir/$outsubj/bem
  mkdir -p $outdir
  echo \
  cp -rpd \
    $SUBJECTS_DIR/$subj'/bem/*.tri' \
    $outdir
  cp -rpd \
    $SUBJECTS_DIR/$subj/bem/*.tri \
    $outdir
endif

# create dip files for resampled subject
if ($dipflag) then
  setenv SUBJECTS_DIR $outsubjdir
  cd $outdir
  foreach hemi ($hemilist)
    set fname_out = $outdir/$hemi'_white.dip'
    if (! -e $fname_out) then
      fs_surfdip $outsubj $hemi
    endif
  end
endif

exit 0;
###############################################

############--------------##################
parse_args:
set cmdline = ($argv);

set subj = $argv[1]; shift;
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "-subjdir":
      if ( $#argv == 0) goto arg1err;
      set subjdir = $argv[1]; shift;
      breaksw
    case "-outsubj":
      if ( $#argv == 0) goto arg1err;
      set outsubj = $argv[1]; shift;
      breaksw
    case "-outsubjdir":
      if ( $#argv == 0) goto arg1err;
      set outsubjdir = $argv[1]; shift;
      breaksw
    case "-ico":
      if ( $#argv == 0) goto arg1err;
      set ico = $argv[1]; shift;
      breaksw
    case "-dip"
      set dipflag = 1
      breaksw
    case "-bem"
      set bemflag = 1
      breaksw
    case "-force"
      set forceflag = 1
      breaksw
    default:
      echo ERROR: Flag $flag unrecognized. 
      echo $cmdline
      exit 1
      breaksw
  endsw

end

goto parse_args_return;
############--------------##################

############--------------##################
check_params:

  if ($#subjdir != 0) then
    setenv SUBJECTS_DIR $subjdir
  endif
  if ($#outsubjdir == 0) then
    set outsubjdir = $SUBJECTS_DIR;
  endif
  if ($#outsubj == 0) then
    set outsubj = $subj'_ico'$ico;
  endif

goto check_params_return;
############--------------##################

############--------------##################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1
############--------------##################

############--------------##################
usage_exit:
  echo " "
  echo "USAGE: fs_subj2ico <subj> [options]"
  echo " "
  echo "Required Arguments:";
  echo "  subj : subject name"
  echo " "
  echo "Optional Arguments:";
  echo '  -subjdir root_subjects_dir : [$SUBJECTS_DIR]'
  echo "  -outsubj  output_subject_name : [subj_ico$ico]"
  echo '  -outsubjdir output_subjects_dir : [$SUBJECTS_DIR]'
  echo "  -ico iconumber : [$ico]"
  echo "  -bem               copy BEM surfaces from bem subdir"
  echo "  -dip               create dip files with tksurfer"
  echo "  -force             overwrite existing output files"
  echo " "
exit 1;
