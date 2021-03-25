#!/bin/csh -f
# surfstats -- view painted stats on a surface with tksurfer (w file format)

## initialize variables ##
set wfile = ()
set name = ()
set hemi = ()
set surf = inflated
set patch = occip
set fthresh = 0.0
set fslope = 1.5
set fmid = 1.5
set fadef = 0.7
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
set fieldsignflag = 0
set saveflag = 0
set scalebarflag = 0
set colscalebarflag = 0
set truncflag = 0
set angoffset = -1
set angcycles = 0
set offset = 0.2
set cvfact = 1.5
set outstem = output
set outdir = `pwd`
set indir = `pwd`
set offscreenflag = 0

if($#argv < 3) then
  goto usage_exit;
  exit 1;
endif

## parse and check params ##
goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

###############################################################################

# create tcl script
set pidlist = `ps | grep surfstats`
set pid = $pidlist[1]
set out = tempsurf.tcl.$pid
rm -f $out
echo "set overlayflag 1" >> $out
echo "set surfcolor 1" >> $out
echo "set avgflag 1" >> $out
echo "set scalebarflag $scalebarflag" >> $out
echo "set colscalebarflag $colscalebarflag" >> $out
if ($fieldsignflag == 1) then
  echo "set complexvalflag 0" >> $out
  echo "set autoscaleflag 1" >> $out
  echo "set colscale 9" >> $out
  echo "set angle_offset 0.0" >> $out
  echo "set angle_cycles 1.0" >> $out
else if ($polarflag == 1) then
  set complexflag = 1
  echo "set complexvalflag 1" >> $out
  echo "set colscale 0" >> $out
  if ($hemi == rh) then
    if ($angoffset == -1) then
      echo "set angle_offset 0.0" >> $out
    else
      echo "set angle_offset $angoffset" >> $out
    endif    
    echo "set revphaseflag 0" >> $out
  else
    if ($angoffset == -1) then
      echo "set angle_offset 0.5" >> $out
    else
      echo "set angle_offset $angoffset" >> $out
    endif    
    echo "set revphaseflag 1" >> $out
  endif
  if ($angcycles == 0) then
    echo "set angle_cycles 2.2" >> $out
  else
    echo "set angle_cycles $angcycles" >> $out
  endif
else if ($eccenflag == 1) then
  set complexflag = 1
  echo "set complexvalflag 1" >> $out
  echo "set colscale 0" >> $out
  if ($angoffset == -1) then
    echo "set angle_offset 0.83" >> $out
  else
    echo "set angle_offset $angoffset" >> $out
  endif    
  if ($angcycles == 0) then
    echo "set angle_cycles 1.0" >> $out
  else
    echo "set angle_cycles $angcycles" >> $out
  endif
  echo "set revphaseflag 0" >> $out
else
  echo "set complexvalflag 0" >> $out
  echo "set colscale 6" >> $out
  echo "set angle_offset 0.0" >> $out
  echo "set angle_cycles 1.0" >> $out
endif

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

if ($fieldsignflag == 1) then
  set fsfile = $wfile-$hemi.fs
  set fmfile = $wfile-$hemi.fm
  echo "setfile fs $indir/$fsfile" >> $out
  echo "setfile fm $indir/$fmfile" >> $out
  echo "read_fieldsign" >> $out
  echo "read_fsmask" >> $out
  echo "set complexvalflag 0" >> $out
else if ($complexflag == 1) then
  set realfile = $wfile'_'r-$hemi.w
  set imagfile = $wfile'_'i-$hemi.w
  echo "setfile val $indir/$imagfile" >> $out
  echo "read_binary_values" >> $out
  echo "smooth_val $smoothsteps" >> $out
  echo "shift_values" >> $out
  echo "setfile val $indir/$realfile" >> $out
  echo "read_binary_values" >> $out
  echo "smooth_val $smoothsteps" >> $out
else
  set realfile = $wfile-$hemi.w
  echo "setfile val $indir/$realfile" >> $out
  echo "read_binary_values" >> $out
  echo "smooth_val $smoothsteps" >> $out
endif

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
    echo " rotate_brain_x $cust_rot_x" >> $out
    set view = "cus"
  endif
  if ($cust_rot_y != 0) then
    echo " rotate_brain_y $cust_rot_y" >> $out
    set view = "cus"
  endif
  if ($cust_rot_z != 0) then
    echo " rotate_brain_z $cust_rot_z" >> $out
    set view = "cus"
  endif
  set rgbname = "$outdir/$outstem-$hemi-$surf-$view.rgb"
else
  echo "restore" >> $out
  if ($cust_rot_z != 0) then
    echo " rotate_brain_z $cust_rot_z" >> $out
  endif
  set rgbname = "$outdir/$outstem-$hemi-$patch.patch.flat.rgb"
endif
if ($trans_x != 0) then
  echo " translate_brain_x $trans_x" >> $out
endif
if ($trans_y != 0) then
  echo " translate_brain_y $trans_y" >> $out
endif
if ($trans_z != 0) then
  echo " translate_brain_z $trans_z" >> $out
endif
echo "do_lighting_model -1 -1 -1 -1" '$offset' >> $out
echo "scale_brain $scale" >> $out
echo "redraw" >> $out

echo "set rgb $rgbname" >> $out
if ($saveflag) then
  echo "save_rgb" >> $out
  echo "exit" >> $out
else if ($offscreenflag) then
  echo "exit" >> $out
endif
tksurfer_offscreen -$name $hemi $surf -tcl $out

exit 0;
###############################################

############--------------##################
parse_args:
set cmdline = ($argv);

set name = $argv[1]; shift;
set wfile = $argv[1]; shift;
set hemi = $argv[1]; shift;

while( $#argv != 0 )

  set flag = $argv[1]; shift;
  
  switch($flag)

    case "-polar"
      set polarflag = 1;
      breaksw

    case "-eccen"
      set eccenflag = 1;
      breaksw

    case "-fieldsign"
      set fieldsignflag = 1;
      breaksw

    case "-surf":
      if ( $#argv == 0) goto arg1err;
      set surf = $argv[1]; shift;
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

    case "-angoffset"
      if ( $#argv == 0) goto arg1err;
      set angoffset = $argv[1]; shift;
      breaksw

    case "-cycles"
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

    case "-savergb"
      set saveflag = 1
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

  if ($#name == 0) then
    goto usage_exit;
  endif

  if ($#wfile == 0) then
    goto usage_exit;
  endif

  if ($view != lat && \
      $view != ven && \
      $view != med && \
      $view != pos && \
      $view != dor && \
      $view != cus) then
    echo "ERROR: view must be lat, ven, med, pos, or dor"
    exit 1;
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
  echo "USAGE: surfstats <name> <wfilestem> <hemi> [options]"
  echo " "
  echo "Required Arguments:";
  echo "  name      : subject name"
  echo "  wfile     : file name, excluding infixes and .w extension"
  echo "  hemi      : hemisphere"
  echo " "
  echo "Optional Arguments:";
  echo "  -polar"
  echo "  -eccen"
  echo "  -fieldsign"
  echo "  -surf      surface              : [$surf]"
  echo "  -fthresh   fthresh              : [$fthresh]"
  echo "  -fslope    fslope               : [$fslope]"
  echo "  -fmid      fmid                 : [$fmid]"
  echo "  -fadef     fadef                : [$fadef]"
  echo "  -smooth    smoothsteps          : [$smoothsteps]"
  echo "  -trunc"
  echo "  -angoffset angoffset            : [$angoffset]"
  echo "  -cycles    angcycles            : [$angcycles]"
  echo "  -offset    offset               : [$offset]"
  echo "  -cvfact    cvfact               : [$cvfact]"
  echo "  -scale     zoom factor          : [$scale]"
  echo "  -view      lat,med,ven,pos,dor  : [$view]"
  echo "  -rotx      x rotation           : [$cust_rot_x]"
  echo "  -roty      y rotation           : [$cust_rot_y]"
  echo "  -rotz      z rotation           : [$cust_rot_z]"
  echo "  -transx    x translation        : [$trans_x]"
  echo "  -transy    y translation        : [$trans_y]"
  echo "  -transz    z translation        : [$trans_z]"
  echo "  -flat"
  echo "  -patch     patch                : [$patch]"
  echo "  -scalebar"
  echo "  -colscalebar"
  echo "  -savergb"
  echo "  -offscreen"
  echo "  -outstem   output rgb file stem : [$outstem]"
  echo "  -outdir    output dir           : [$outdir]"
  echo "  -indir     input dir            : [$indir]"
  echo " "
exit 1;
