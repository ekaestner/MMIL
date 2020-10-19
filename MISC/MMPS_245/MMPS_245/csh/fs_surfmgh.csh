#!/bin/csh -f

##############################################
#
# fs_surfmgh -- view painted stats on a surface with tksurfer (mgh file format)
#
# edit log:
#   bcipolli on 06/13/2007
#   bcipolli on 09/21/2007
#   bcipolli on 10/07/2007
#   dhagler  on 10/08/2007
#   dhagler  on 01/23/2008
#   dhagler  on 05/30/2008
#   dhagler  on 06/03/2008
#   dhagler  on 09/23/2008
#   dhagler  on 05/11/2009   revphase
#   dhagler  on 11/04/2009   padding ndigits
#   dhagler  on 11/05/2009   rotate
#   dhagler  on 08/24/2010   default colscale = 1
#   dhagler  on 07/07/2011   ant view
#   dhagler  on 11/05/2012   mgz
#
##############################################

# TO DO:
# outstem - make default outstem as input filename
# allow list of files; load each using tcl

# initialize variables ##
set mghfile = ()
set subj = ()
set hemi = ()
set surf = inflated
set annot = ''
set label = ''
set patch = occip
set fthresh = 0.0
set fslope = 1.5
set fmid = 1.5
set fadef = 0.7
set colscale = -1
set default_colscale = 1 # heat
set smoothsteps = 0
set scale = 1.0
set view = lat
set cust_rot_x = 0
set cust_rot_y = 0
set cust_rot_z = 0
set trans_x = 0
set trans_y = 0
set trans_z = 0
set cust_rot_z = 0
set flatflag = 0
set polarflag = 0
set eccenflag = 0
set complexflag = 0
set revphaseflag = 0
set angoffset = -1
set angcycles = -1
set saveflag = 0
set tiff_flag = 0
set scalebarflag = 0
set colscalebarflag = 0
set truncflag = 0
set offset = 0.2
set cvfact = 1.5
set outstem = output
set outdir = `pwd`
set indir = "";
set offscreenflag = 0
set rmintermedflag = 0
set silenceflag = 0
set tfirst = 0
set tlast = -1
set tmax = -1
set ndigits = 1
set forceflag = 0

if ( $?SUBJECTS_DIR ) then
  set sdir = "$SUBJECTS_DIR";
endif

if($#argv < 3) then
  goto usage_exit;
  exit 1;
endif

##############################################

## parse and check params ##
goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

# get nframes from mri_info
set nframes = `mri_info $mghfile --nframes`

# figure out the # of digits to pad
if ($nframes > 9999) then
  set ndigits = 5
else if ($nframes > 999) then
  set ndigits = 4
else if ($nframes > 99) then
  set ndigits = 3
else if ($nframes > 9) then
  set ndigits = 2
else
  set ndigits = 1
endif

##############################################

# check dependencies

if ( "`which tksurfer_offscreen | grep /`" != "" ) then
  setenv TKSURFER_PROG "tksurfer_offscreen"
else if ( "`which tksurfer | grep /`" != "" ) then
  setenv TKSURFER_PROG "tksurfer"
endif

# Pick the appropriate tksurfer
if ( ! $offscreenflag )  then
  set tksurfer_local = tksurfer;
else if ( $?TKSURFER_PROG ) then
  set tksurfer_local = $TKSURFER_PROG;
else
  echo "  ** WARNING: using local version of tksurfer to render offscreen; may not work!";
  set tksurfer_local = tksurfer;
endif

# Make sure the one we chose exists
if ( "`which $tksurfer_local | grep /`" == "" ) goto no_tksurfer_err;

###############################################################################

# set environment
if ( $?SUBJECTS_DIR ) then
  set oldsdir = "$SUBJECTS_DIR";
  set oldsilence = "$SUBJECTS_DIR";
endif

setenv SUBJECTS_DIR $sdir;
setenv FS_FREESURFERENV_NO_OUTPUT $silenceflag

# create tcl script
set pidlist = `ps | grep surfmgh`
set pid = $pidlist[1]
set out = tempsurf.tcl.$pid
rm -f $out
echo "set overlayflag 1" >> $out
echo "set surfcolor 1" >> $out
echo "set avgflag 1" >> $out
echo "set scalebarflag $scalebarflag" >> $out
echo "set colscalebarflag $colscalebarflag" >> $out
if ($complexflag == 1) then
  echo "set complexvalflag 1" >> $out
endif
if ($polarflag == 1) then
  if ($hemi == rh) then
    if ($angoffset == -1) then
      echo "set angle_offset 0.0" >> $out
    else
      echo "set angle_offset $angoffset" >> $out
    endif    
    if ($revphaseflag) then
      echo "set revphaseflag 1" >> $out
    else
      echo "set revphaseflag 0" >> $out
    endif
  else
    if ($angoffset == -1) then
      echo "set angle_offset 0.5" >> $out
    else
      echo "set angle_offset $angoffset" >> $out
    endif    
    if ($revphaseflag) then
      echo "set revphaseflag 0" >> $out
    else
      echo "set revphaseflag 1" >> $out
    endif
  endif
  if ($angcycles == -1) then
    echo "set angle_cycles 2.2" >> $out
  else
    echo "set angle_cycles $angcycles" >> $out
  endif
  if ($colscale == -1) then
    set colscale = 0
  endif
else if ($eccenflag == 1) then
  if ($angoffset == -1) then
    echo "set angle_offset 0.83" >> $out
  else
    echo "set angle_offset $angoffset" >> $out
  endif    
  if ($angcycles == -1) then
    echo "set angle_cycles 1.0" >> $out
  else
    echo "set angle_cycles $angcycles" >> $out
  endif
  echo "set revphaseflag $revphaseflag" >> $out
  if ($colscale == -1) then
    set colscale = 0
  endif
else if ($complexflag == 1) then
  if ($angoffset == -1) then
    echo "set angle_offset 0.0" >> $out
  else
    echo "set angle_offset $angoffset" >> $out
  endif
  if ($angcycles == -1) then
    echo "set angle_cycles 1.0" >> $out
  else
    echo "set angle_cycles $angcycles" >> $out
  endif
  echo "set revphaseflag $revphaseflag" >> $out
  if ($colscale == -1) then
    set colscale = 0
  endif
else
  echo "set complexvalflag 0" >> $out
  echo "set angle_offset 0.0" >> $out
  echo "set angle_cycles 1.0" >> $out
  if ($colscale == -1) then
    set colscale = $default_colscale
  endif
endif
echo "set colscale $colscale" >> $out
echo "set fthresh $fthresh" >> $out
echo "set fslope $fslope" >> $out
echo "set fmid   $fmid" >> $out
echo "set truncphaseflag $truncflag" >> $out
echo "set flatzrot 0" >> $out
echo "set flatscale 1.0" >> $out
echo "set smoothsteps $smoothsteps" >> $out
echo "set offset $offset" >> $out
echo "set cvfact $cvfact" >> $out
echo "set fadef $fadef" >> $out

if ($flatflag == 1) then
  echo "setfile patch $hemi.$patch.patch.flat" >> $out
  echo "read_binary_patch" >> $out
endif

echo "read_binary_curv" >> $out

if ($offscreenflag == 1) then
  echo "set renderoffscreen 1" >> $out
endif

echo "open_window" >> $out
if ($flatflag == 0) then
  echo "make_lateral_view" >> $out
  if ($view == med) then
    echo "rotate_brain_y 180" >> $out
  else if ($view == ven) then
    echo "rotate_brain_x 90" >> $out
  else if ($view == pos) then
    echo 'if {$hemi == "rh"} {' >> $out
    echo "  rotate_brain_y 270" >> $out
    echo "} else {" >> $out
    echo "  rotate_brain_y 90" >> $out
    echo "}" >> $out
  else if ($view == ant) then
    echo 'if {$hemi == "rh"} {' >> $out
    echo "  rotate_brain_y 90" >> $out
    echo "} else {" >> $out
    echo "  rotate_brain_y 270" >> $out
    echo "}" >> $out
  else if ($view == dor) then
    echo 'if {$hemi == "rh"} {' >> $out
    echo "  rotate_brain_y 270" >> $out
    echo "  rotate_brain_x 270" >> $out
    echo "} else {" >> $out
    echo "  rotate_brain_y 90" >> $out
    echo "  rotate_brain_x 270" >> $out
    echo "}" >> $out
  endif
  if ($cust_rot_x != 0) then
    echo "rotate_brain_x $cust_rot_x" >> $out
    set view = "cus"
  endif
  if ($cust_rot_y != 0) then
    echo "rotate_brain_y $cust_rot_y" >> $out
    set view = "cus"
  endif
  if ($cust_rot_z != 0) then
    echo "rotate_brain_z $cust_rot_z" >> $out
    set view = "cus"
  endif
else
  echo "restore" >> $out
  if ($cust_rot_z != 0) then
    echo "rotate_brain_z $cust_rot_z" >> $out
  endif
endif
if ($trans_x != 0) then
  echo "translate_brain_x $trans_x" >> $out
endif
if ($trans_y != 0) then
  echo "translate_brain_y $trans_y" >> $out
endif
if ($trans_z != 0) then
  echo "translate_brain_z $trans_z" >> $out
endif
echo "do_lighting_model -1 -1 -1 -1" '$offset' >> $out
echo "scale_brain $scale" >> $out

if ($annot != "") then
  echo "set labelstyle 1" >> $out
  echo "labl_import_annotation $annot" >> $out
endif

if ($label != "") then
  echo "set labelstyle 1" >> $out
  echo "labl_load $label" >> $out
endif

if ($tlast == -1) then
  if ($nframes <= $tmax || $tmax == -1) then
    @ tlast = $nframes - 1
  else
    set tlast = $tmax
  endif
endif

echo "sclv_read_from_volume 0 $mghfile 2" >> $out
if ($complexflag == 1) then
  echo "sclv_read_from_volume 1 $imagfile 2" >> $out
endif

set mustrun = 0
set f = $tfirst
while ($f <= $tlast)
  set fstem = "$outdir/$outstem"
  if ($nframes>1) then
    set fpadded    = `count -dig $ndigits $f $f`
    set fstem = "$fstem-$fpadded"
  endif
  set fstem = "$fstem-$hemi"
  if ($flatflag) then
    set fstem = "$fstem-$patch.patch.flat"
  else
    set fstem = "$fstem-$surf-$view"
  endif
  set rgbname = "$fstem.rgb"
  set tiffname = "$fstem.tif"

  # echo these generation commands only if:
  #  - user is forcing it (forceflag)  OR
  #  - user is doing display only (!saveflag)  OR
  #  - output file does not exist
  #
  if ( ! $saveflag || $forceflag || ( $saveflag && ! -e $rgbname ) || ( $saveflag && $tiff_flag && ! -e $tiffname ) ) then
    echo "sclv_set_current_field 0" >> $out
    echo "sclv_set_current_timepoint $f 0" >> $out
    if ($smoothsteps > 0) then
#      echo "sclv_smooth $smoothsteps 0" >> $out
      echo "smooth_val $smoothsteps" >> $out
    endif

    if ($complexflag == 1) then
      echo "sclv_set_current_field 1" >> $out
      echo "sclv_set_current_timepoint $f 0" >> $out
      if ($smoothsteps > 0) then
        echo "swap_val_val2" >> $out      
        echo "smooth_val $smoothsteps" >> $out
        echo "swap_val_val2" >> $out      
# NOTE: sclv_smooth does not work with tksurfer_offscreen
#        echo "sclv_smooth $smoothsteps 1" >> $out
      endif
    endif

    echo "set fthresh $fthresh" >> $out
    echo "set fslope $fslope" >> $out
    echo "set fmid   $fmid" >> $out

    echo "redraw" >> $out
	
    if ($saveflag) then
      echo "set rgb $rgbname" >> $out
      echo "save_rgb" >> $out
    endif
	
    if ($tiff_flag) then
     echo "save_tiff $tiffname" >> $out
    endif

    set mustrun = 1
  endif

  @ f++
end

if ($saveflag || $offscreenflag || $tiff_flag) then
  echo "exit" >> $out
endif

# Now, we've finished generating the script
#   and we're ready to run; check to see
#   if there's actually any NEED to run!

if ( ($saveflag || $tiff_flag || $offscreenflag) && ! $mustrun ) then
  if (! $silenceflag) echo "All rgb files already exist; exiting.";

else
  if ($silenceflag) then
    set silentarg = " 2&>1 silence.log";
  else
    set silentarg = "";
  endif

## todo: will this work with multi-frame data?
#  $tksurfer_local -$subj $hemi $surf -overlay $mghfile -tcl $out $silentarg
  $tksurfer_local -$subj $hemi $surf -tcl $out $silentarg
endif

#  Do the necessary cleanup
if ($rmintermedflag && ($saveflag || $offscreenflag)) then
  if ( -e $out )        rm -f $out
  if ( -e surfer.log )  rm -f surfer.log
  if ( -e silence.log ) rm -f silence.log
endif

# unset environment
if ( $?oldsdir ) then
  setenv SUBJECTS_DIR $oldsdir;
  setenv FS_FREESURFERENV_NO_OUTPUT $oldsilence
endif

exit 0;
###############################################

############--------------##################
parse_args:
set cmdline = ($argv);

set subj    = $argv[1]; shift;
set mghfile = $argv[1]; shift;
set hemi    = $argv[1]; shift;

while( $#argv != 0 )

  set flag = $argv[1]; shift;

  switch($flag)

    case "-imagfile"
      if ( $#argv == 0) goto arg1err;
      set imagfile = $argv[1]; shift;
      set complexflag = 1;
      breaksw

    case "-polar"
      set polarflag = 1;
      set complexflag = 1;
      breaksw

    case "-eccen"
      set eccenflag = 1;
      set complexflag = 1;
      breaksw

    case "-surf":
      if ( $#argv == 0) goto arg1err;
      set surf = $argv[1]; shift;
      breaksw

    case "-colscale":
  	  if ( $#argv == 0) goto arg1err;
  	  set colscale = $argv[1]; shift;
	  breaksw
	  
    case "-fthresh"
      if ( $#argv == 0) goto arg1err;
      set fthresh = $argv[1]; shift;
      breaksw

    case "-fslope"
      if ( $#argv == 0) goto arg1err;
      set fslope = $argv[1]; shift;
      breaksw

    case "-fmid"
      if ( $#argv == 0) goto arg1err;
      set fmid = $argv[1]; shift;
      breaksw

    case "-fadef"
      if ( $#argv == 0) goto arg1err;
      set fadef = $argv[1]; shift;
      breaksw

    case "-smooth"
      if ( $#argv == 0) goto arg1err;
      set smoothsteps = $argv[1]; shift;
      breaksw

    case "-revphase"
      set revphaseflag = 1;
      breaksw

    case "-angoffset"
      if ( $#argv == 0) goto arg1err;
      set angoffset = $argv[1]; shift;
      breaksw

    case "-angcycles"
      if ( $#argv == 0) goto arg1err;
      set angcycles = $argv[1]; shift;
      breaksw

    case "-offset"
      if ( $#argv == 0) goto arg1err;
      set offset = $argv[1]; shift;
      breaksw

    case "-cvfact"
      if ( $#argv == 0) goto arg1err;
      set cvfact = $argv[1]; shift;
      breaksw

    case "-scale"
      if ( $#argv == 0) goto arg1err;
      set scale = $argv[1]; shift;
      breaksw

    case "-view"
      if ( $#argv == 0) goto arg1err;
      set view = $argv[1]; shift;
      breaksw

    case "-rotx"
      if ( $#argv == 0) goto arg1err;
      set cust_rot_x = $argv[1]; shift;
      breaksw

    case "-roty"
      if ( $#argv == 0) goto arg1err;
      set cust_rot_y = $argv[1]; shift;
      breaksw

    case "-rotz"
      if ( $#argv == 0) goto arg1err;
      set cust_rot_z = $argv[1]; shift;
      breaksw

    case "-transx"
      if ( $#argv == 0) goto arg1err;
      set trans_x = $argv[1]; shift;
      breaksw

    case "-transy"
      if ( $#argv == 0) goto arg1err;
      set trans_y = $argv[1]; shift;
      breaksw

    case "-transz"
      if ( $#argv == 0) goto arg1err;
      set trans_z = $argv[1]; shift;
      breaksw

    case "-trunc"
      set truncflag = 1;
      breaksw

    case "-flat"
      set flatflag = 1;
      breaksw

    case "-patch"
      if ( $#argv == 0) goto arg1err;
      set patch = $argv[1]; shift;
      breaksw
	  
    case "-annot"
      if ( $#argv == 0) goto arg1err;
      set annot = $argv[1]; shift;
      breaksw

    case "-label"
      if ( $#argv == 0) goto arg1err;
      set label = $argv[1]; shift;
      breaksw

    case "-savergb"
      set saveflag = 1
      breaksw

    case "-savetiff"
      set tiff_flag  = 1
      breaksw

    case "-force"
      set forceflag = 1
      breaksw

    case "-offscreen"
      set offscreenflag = 1
      breaksw

    case "-scalebar"
      set scalebarflag = 1
      breaksw

    case "-colscalebar"
      set colscalebarflag = 1
      breaksw

    case "-rmintermed"
      set rmintermedflag = 1;
      breaksw

    case "-silent"
      set silenceflag = 1;
      breaksw

    case "-outstem"
      if ( $#argv == 0) goto arg1err;
      set outstem = $argv[1]; shift;
      breaksw

    case "-outdir"
      if ( $#argv == 0) goto arg1err;
      set outdir = $argv[1]; shift;
      breaksw

    case "-indir"
      if ( $#argv == 0) goto arg1err;
      set indir = $argv[1]; shift;
      breaksw

    case "-frame"
      if ( $#argv == 0) goto arg1err;
      set tfirst = $argv[1]; shift;
      set tlast = $tfirst;
      breaksw

    case "-tfirst"
      if ( $#argv == 0) goto arg1err;
      set tfirst = $argv[1]; shift;
      breaksw

    case "-tlast"
      if ( $#argv == 0) goto arg1err;
      set tlast = $argv[1]; shift;
      breaksw

    case "-tmax"
      if ( $#argv == 0) goto arg1err;
      set tmax = $argv[1]; shift;
      breaksw

    case "-sdir"
      if ( $#argv == 0) goto arg1err;
      set sdir = $argv[1]; shift;
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

  if ("$hemi" != "rh" && "$hemi" != "lh") then
    echo "ERROR: hemi must be rh or lh"
    exit 1;
  endif

  if ($#subj == 0) then
    goto usage_exit;
  endif

  ## check mghfile ##
  if ($#mghfile == 0) then
    goto usage_exit;
  endif
  # allow both hemi-free/extension-free as well as
  # fully qualified filename versions
  if ("`echo $mghfile|grep .mgz`" == "") then
    if ("`echo $mghfile|grep .mgh`" == "") then
      set mghfile=$mghfile-$hemi.mgh
    endif
  endif
  # add the optional 'indir'
  if ($indir != "") then
    set mghfile=$indir/$mghfile
  endif
  # verify the file exists!
  if ( ! -e $mghfile ) then
    echo "ERROR: input file doesn't exist: $mghfile";
    exit 3;
  endif

  if ($complexflag) then
    ## check imagfile ##
    if ($#imagfile == 0) then
      echo "ERROR: imagfile must be supplied if -polar or -eccen are used"
      exit 1;
    endif
    if ("`echo $imagfile|grep .mgh`" == "") then
      set imagfile=$imagfile-$hemi.mgh
    endif
    if ($indir != "") then
      set imagfile=$indir/$imagfile
    endif
    if ( ! -e $imagfile ) then
      echo "ERROR: input file doesn't exist: $imagfile";
      exit 3;
    endif
  endif
    
  if ($view != lat && \
      $view != ven && \
      $view != med && \
      $view != ant && \
      $view != pos && \
      $view != dor && \
      $view != cus) then
    echo "ERROR: view must be lat, ven, med, pos, or dor"
    exit 1;
  endif

  if ( ! $?sdir || $sdir == "" ) then
    echo "ERROR: SUBJECTS_DIR not set and no 'sdir' flag... quitting";
    exit 7
  else if ( ! -e $sdir ) then
    echo "ERROR: SUBJECTS_DIR does not exist: $sdir";
    exit 7
  endif
  
  # colscale must be between 0 and 8
  if ( $colscale < -1 || $colscale > 8 ) then
    echo "ERROR: colscale argument must be between 0 and 8"
	exit 17
  endif

goto check_params_return;
############--------------##################

############--------------##################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1

no_tksurfer_err:
  echo "ERROR: couldn't locate any version of tksurfer; both not in your
  echo "       path and also not in the expected mmildev shared path.";
  echo " ";
  echo "You need to resolve this issue before you can render images.";
  exit 2;

############--------------##################

############--------------##################
usage_exit:
  echo " "
  echo "USAGE: surfmgh <subj> <mghfile> <hemi> [options]"
  echo " "
  echo "Required Arguments:";
  echo "  <subj>     : subject name"
  echo "  <mghfile>  : file name, optionally excluding hemi infix and .mgh extension"
  echo "  <hemi>     : hemisphere"
  echo " "
  echo "Optional Arguments:";
  echo "  [I/O args:]"
  echo "  -imagfile     file name for imaginary component"
  echo "    (only required if -polar or -eccen switches are used)"
  echo "    (real component comes from mghfile parameter)"
  echo "  -indir        input dir              : [$indir]"
  echo "  -outdir       output dir             : [$outdir]"
  echo "  -outstem      output file stem       : [$outstem]"
  echo "  -sdir         SUBJECTS_DIR           : [$sdir]"
  echo " "
  echo "  [Time args:]"
  echo "  -tfirst       first time point       : [$tfirst]"
  echo "  -tlast        last time point        : [$tlast]"
  echo " "
  echo "  [Surface file args:]"
  echo "  -annot        annotation file        : [$annot]"
  echo "  -label        label file             : [$label]"
  echo "  -surf         surface (e.g. inflated): [$surf]"
  echo "  -flat         use flat patch"
  echo "  -patch        patch name (e.g. full) : [$patch]"
  echo " "
  echo "  [Image args:]";
  echo "  -colscale     coloring palette       : [$colscale]"
  echo "     1=heat,  2=cyan-red,    3=blue-green-red"
  echo "     5=gray,  6=blue-red,    7=green-red"
  echo "     0=rgb,   4=gr (2 cond), 8=rygb (cx)  (for complex only)"
  echo "  -fthresh      fthresh                : [$fthresh]"
  echo "  -fmid         fmid                   : [$fmid]"
  echo "  -fslope       fslope                 : [$fslope]"
  echo "  -colscalebar  show colorbar"
  echo "  -scalebar     show scalebar"
  echo "  -cvfact       gyrus/sulcus ratio     : [$cvfact]"
  echo "  -offset       bg surface intensity   : [$offset]"
  echo "  -smooth       smoothsteps            : [$smoothsteps]"
  echo "  -trunc        set negative values to zero"
  echo "  -revphase     reverse phase (for cx)"
  echo "  -angoffset    phase offset (for cx)  : [$angoffset]"
  echo "  -angcycles    phase cycles (for cx)  : [$angcycles]"
  echo "  -fadef        fade factor (for cx)   : [$fadef]"
  echo " "
  echo "  [View args:]"
  echo "  -view         lat,med,ven,ant,pos,dor: [$view]"
  echo "  -scale        zoom factor            : [$scale]"
  echo "  -rotx         x rotation             : [$cust_rot_x]"
  echo "  -roty         y rotation             : [$cust_rot_y]"
  echo "  -rotz         z rotation             : [$cust_rot_z]"
  echo "  -transx       x translation          : [$trans_x]"
  echo "  -transy       y translation          : [$trans_y]"
  echo "  -transz       z translation          : [$trans_z]"
  echo " "
  echo "  [Program args:]";
  echo "  -offscreen    render images offscreen"
  echo "  -savergb      save rendered images to rgb"
  echo "  -savetiff     save rendered images to tiff"
  echo "  -force        force remake of cached images"
  echo "  -rmintermed   remove intermediate files (log files, etc.)"
  echo "  -silent       suppresses tksurfer output";
  echo " ";
exit 1;
