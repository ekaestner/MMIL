#!/bin/csh -f

###############################################################
# mmil_prep_MRI.csh -- copy freesurfer recon to analysis directory
#   along with BEM surfaces, optionally make dip/dec files
#
# created:  10/22/09 by Don Hagler
# last mod: 04/05/10 by Don Hagler
###############################################################

# initialize variables ##
set proc_rootdir = ()
set fsurf_rootdir = ()
set out_rootdir = ()
set subj_ID = ()
set sess_ID = ()
set ico = ()
set hemilist = (lh rh)
set forceflag = 0
set graymidflag = 0
set dipdecflag = 0
set bemtype = 'T1'

if($#argv < 3) then
  goto usage_exit;
  exit 1;
endif

## parse and check params ##
goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

###########################################################

set subj = 'FREESURFERRECON_'$subj_ID
set outsubj = 'FREESURFERRECON_'$subj_ID'_ico'$ico
set subjdir = $out_rootdir/$subj_ID
setenv SUBJECTS_DIR $subjdir

# copy native space freesurfer recon
if (! -e $out_rootdir/$subj_ID) then
  mkdir -p $out_rootdir/$subj_ID
endif

if (! -e $out_rootdir/$subj_ID/$subj) then
  echo \
  cp -rpd \
    $fsurf_rootdir/$subj'_'$sess_ID \
    $out_rootdir/$subj_ID/$subj
  cp -rpd \
    $fsurf_rootdir/$subj'_'$sess_ID \
    $out_rootdir/$subj_ID/$subj
endif

if ($ico > 0) then
  fs_subj2ico.csh $subj -subjdir $subjdir -outsubj $outsubj -ico $ico

  # copy bem surfaces for resampled subject
  set outdir = $subjdir/$outsubj/bem
  mkdir -p $outdir
  if (! -e $outdir/make_bem_surfs.csh) then
    echo \
    cp -rpd \
      $proc_rootdir/MRIPROCESSED_$subj_ID'_'$sess_ID/'bem_'$bemtype'/*' \
      $outdir
    cp -rpd \
      $proc_rootdir/MRIPROCESSED_$subj_ID'_'$sess_ID/'bem_'$bemtype/* \
      $outdir
  endif
endif

# copy bem surfaces for native subject
set outdir = $subjdir/$subj/bem
mkdir -p $outdir
if (! -e $outdir/make_bem_surfs.csh) then
  echo \
  cp -rpd \
    $proc_rootdir/MRIPROCESSED_$subj_ID'_'$sess_ID/'bem_'$bemtype'/*' \
    $outdir
  cp -rpd \
    $proc_rootdir/MRIPROCESSED_$subj_ID'_'$sess_ID/'bem_'$bemtype/* \
    $outdir
endif

if ($graymidflag) then
  # create graymid files
  set outdir = $subjdir/$subj/surf
  cd $outdir
  foreach hemi ($hemilist)
    set fname_out = $outdir/$hemi'.graymid'
    if (! -e $fname_out) then
      mris_expand -thickness $hemi.white 0.5 $hemi.graymid
    endif
  end
endif

# create dip/dec files
if ($dipdecflag) then
  # for native subject
  set outdir = $subjdir/$subj/bem
  cd $outdir
  foreach hemi ($hemilist)
    set fname_out = $outdir/$hemi'_white.dip'
    if (! -e $fname_out) then
      fs_surfdip $subj $hemi
    endif
    set fname_out = $outdir/$hemi'_white_7.dec'
    if (! -e $fname_out) then
      fs_surfdec $subj $hemi
    endif
  end

  # for native graymid
  if ($graymidflag) then
    set outdir = $subjdir/$subj/bem
    cd $outdir
    foreach hemi ($hemilist)
      set fname_out = $outdir/$hemi'_graymid.dip'
      if (! -e $fname_out) then
        fs_surfdip $subj $hemi -surf graymid
      endif
      set fname_out = $outdir/$hemi'_graymid_7.dec'
      if (! -e $fname_out) then
        fs_surfdec $subj $hemi -surf graymid
      endif
    end
  endif

  # for ico
  if ($ico > 0) then
    set outdir = $subjdir/$outsubj/bem
    cd $outdir
    foreach hemi ($hemilist)
      set fname_out = $outdir/$hemi'_white.dip'
      if (! -e $fname_out) then
        fs_surfdip $outsubj $hemi -surf white
      endif
    end
  endif
endif

exit 0


############--------------##################
parse_args:
set cmdline = ($argv);
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "-proc_rootdir"
      if ( $#argv == 0) goto arg1err;
      set proc_rootdir = $argv[1]; shift;
      breaksw

    case "-fsurf_rootdir"
      if ( $#argv == 0) goto arg1err;
      set fsurf_rootdir = $argv[1]; shift;
      breaksw

    case "-out_rootdir"
      if ( $#argv == 0) goto arg1err;
      set out_rootdir = $argv[1]; shift;
      breaksw

    case "-subj_ID"
      if ( $#argv == 0) goto arg1err;
      set subj_ID = $argv[1]; shift;
      breaksw

    case "-sess_ID"
      if ( $#argv == 0) goto arg1err;
      set sess_ID = $argv[1]; shift;
      breaksw

    case "-ico"
      if ( $#argv == 0) goto arg1err;
      set ico = $argv[1]; shift;
      breaksw

    case "-bemtype"
      if ( $#argv == 0) goto arg1err;
      set bemtype = $argv[1]; shift;
      breaksw

    case "-force"
      set forceflag = 1
      breaksw

    case "-graymid"
      set graymidflag = 1
      breaksw

    case "-dipdec"
      set dipdecflag = 1
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
  echo "USAGE: mmil_prep_MRI.csh -[options] [values]"
  echo " "
  echo "Required Arguments (with values):";
  echo "  -proc_rootdir <proc_rootdir>  source proc rootdir"
  echo "  -fsurf_rootdir <proc_rootdir> source fsurf rootdir"
  echo "  -out_rootdir <rootdir>        output rootdir
  echo "  -subj_ID <subjID>             subject ID"
  echo "  -sess_ID <sessID>             session ID"
  echo " "
  echo "Optional Arguments (with values):";
  echo "  -ico <ico num>            icosahedral order"
  echo "     if 0, do not resample to ico
  echo "     { default = 0 }
  echo "  -bemtype <bemtype>        source of BEM surfaces"
  echo "     { default = T1 }
  echo " "
  echo "Optional Arguments (flags):";
  echo "  -graymid                  create gray mid surfaces"
  echo "  -dipdec                   create dip and dec files with tksurfer"
  echo "  -force                    overwrite existing output files"
  echo " ";
exit 1;
