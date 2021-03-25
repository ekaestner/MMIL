#!/bin/csh -f
#
# render_AtlasTrack.csh -- render images of AtlasTrack fiber streamlines
#
# created:  04/09/12 by Don Hagler
# last mod: 04/10/12 by Don Hagler
#
###############################################################################

# todo: option to directly specify fiber numbers (e.g. -f 101 -f 102 -f 103)
# todo: hard-code colors for each fiber
# todo: create temporary color.dat with colors only for included fnums,
#       excluding missing fibers
# todo: options to direcly specify fname_image and fiberdir
# todo: add parameters in init_parms to parse_args and usage_exit

###############################################################################

## initialize parameters default values ##
goto init_parms;
init_parms_return:

if($#argv < 1) then
  goto usage_exit;
  exit 1;
endif

## parse input arguments ##
goto parse_args;
parse_args_return:

## check parameters ##
goto check_parms;
check_parms_return:

## find background image file ##
goto find_underlay;
find_underlay_return:

## find fiber files ##
goto find_fibers;
find_fibers_return:

## run tractoview ##
goto run_tractoview;
run_tractoview_return:

exit 0;

###############################################################################

###############################################################################
init_parms:

  set indir = ()
  set save_flag = 1
  set fname_tif = 'fibers.tif'
  set fname_log = 'fibers.log'
  set fname_color = 'color.dat'
  set ulay = 'none'
  set ulay_plane = 'sag'
  set view = 'L'
  set xrot = 0
  set yrot = 0
  set zrot = 0
  set xflip = 0
  set yflip = 0
  set zflip = 1
  set image_alpha = 0.4
  set fiber_alpha = 0.4

  ## todo: add to parse_args and usage_exit:
  set fnums = (\
    '101' '102' '103' '104' '105' '106' '107' '108' \
    '109' '110' '115' '116' '117' '118' '119' '120' \
	  '133' '134' '135' '136' '137' '138' '141' '142' \
    '143' '144' '145' '146' '147' '148' '149' '150' \
  )

# fornix
# cingulum (cing)
# cingulum (ph)
# CST
# ATR
# uncinate
# ILF
# IFOF
# SLF

#  set fnums = (\
#    '101' '102' \
#    '103' '104' \
#    '105' '106' \
#    '107' '108' \
#    '109' '110' \
#    '115' '116' \
#    '117' '118' \
#    '119' '120' \
#	  '133' '134' \
#  )

  set ccflag = 1
  set fiber_suffix = '_pthresh0.08_threshFA0.00_minflen12_countatlas_path.grp'
  set fname_ovly = ()
  set ovly_offset = 0
  set ovly_slope = 1
  set ovly_thresh = 0
  set fname_alpha = ()
  set min_fiber_alpha = 0.1
  set rendcor = 1
  set rendsag = 1
  set rendhor = 1
  set forceflag = 1

goto init_parms_return;
###############################################################################

###############################################################################
parse_args:

  set cmdline = ($argv);
  set indir = $argv[1]; shift;

  while( $#argv != 0 )
    set flag = $argv[1]; shift;
    switch($flag)

      case "-save_flag":
        if ( $#argv == 0) goto arg1err;
        set save_flag = $argv[1]; shift;
        breaksw

      case "-fname_tif":
        if ( $#argv == 0) goto arg1err;
        set fname_tif = $argv[1]; shift;
        breaksw

      case "-fname_log":
        if ( $#argv == 0) goto arg1err;
        set fname_log = $argv[1]; shift;
        breaksw

      case "-fname_color":
        if ( $#argv == 0) goto arg1err;
        set fname_color = $argv[1]; shift;
        breaksw

      case "-ulay":
        if ( $#argv == 0) goto arg1err;
        set ulay = $argv[1]; shift;
        breaksw

      case "-ulay_plane":
        if ( $#argv == 0) goto arg1err;
        set ulay_plane = $argv[1]; shift;
        breaksw

      case "-view":
        if ( $#argv == 0) goto arg1err;
        set view = $argv[1]; shift;
        breaksw

      case "-xrot":
        if ( $#argv == 0) goto arg1err;
        set xrot = $argv[1]; shift;
        breaksw

      case "-yrot":
        if ( $#argv == 0) goto arg1err;
        set yrot = $argv[1]; shift;
        breaksw

      case "-zrot":
        if ( $#argv == 0) goto arg1err;
        set zrot = $argv[1]; shift;
        breaksw

      case "-xflip":
        if ( $#argv == 0) goto arg1err;
        set xflip = $argv[1]; shift;
        breaksw

      case "-yflip":
        if ( $#argv == 0) goto arg1err;
        set yflip = $argv[1]; shift;
        breaksw

      case "-zflip":
        if ( $#argv == 0) goto arg1err;
        set zflip = $argv[1]; shift;
        breaksw

      case "-image_alpha":
        if ( $#argv == 0) goto arg1err;
        set image_alpha = $argv[1]; shift;
        breaksw

      case "-fiber_alpha":
        if ( $#argv == 0) goto arg1err;
        set fiber_alpha = $argv[1]; shift;
        breaksw

      case "-fiber_suffix":
        if ( $#argv == 0) goto arg1err;
        set fiber_suffix = $argv[1]; shift;
        breaksw

      default:
        echo ERROR: flag $flag unrecognized. 
        echo $cmdline
        exit 1
        breaksw
    endsw
  end

goto parse_args_return;
###############################################################################

###############################################################################
check_parms:

  if ($#indir == 0) then
    goto usage_exit;
  endif
  if (! -e $indir) then
    echo "ERROR: input directry $indir not found"
    exit 1;
  endif

  if (! -e $fname_color) then
    echo "WARNING: color file $fname_color not found"
    set color_parms = ()
  else
    set color_parms = "-fname_color $fname_color"
  endif

  # todo: do not add if already included in fnums
  # do not add if already included
  # if ("$path" !~ *$dir*) set path = ($dir $path)
  if ($ccflag == 1) then
    set fnums = ($fnums '123')
  else if ($ccflag == 2) then
    set fnums = ($fnums '121' '122')
  endif
  # or, remove duplicates
  # set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`

  if ($save_flag == 1) then
    set save_parms = "-savetif -fname_tif $fname_tif"
  else
    set save_parms = ()
  endif

  set flip_parms = ()
  if ($xflip == 1) then
    set flip_parms = ($flip_parms -xflip)
  else
    set flip_parms = ($flip_parms -noxflip)
  endif
  if ($yflip == 1) then
    set flip_parms = ($flip_parms -yflip)
  else
    set flip_parms = ($flip_parms -noyflip)
  endif
  if ($zflip == 1) then
    set flip_parms = ($flip_parms -zflip)
  else
    set flip_parms = ($flip_parms -nozflip)
  endif

goto check_parms_return;
###############################################################################

###############################################################################
arg1err:

  echo "ERROR: flag $flag requires one argument"

exit 1;
###############################################################################

###############################################################################
usage_exit:

  echo " "
  echo "USAGE: render_AtlasTrack.csh <indir> [options]"
  echo " "
  echo "Required Arguments:"
  echo "  indir : input directory containing AtlasTrack output"
  echo " "
  echo "Optional Arguments:"
  echo "  -save_flag <save_flag>     default = $save_flag"
  echo "     [0|1] whether to save rendered image (otherwise display only)"
  echo "  -fname_tif <fname_tif>     default = $fname_tif"
  echo "     output file name for rendered image (tif format)"
  echo "  -fname_log <fname_log>     default = $fname_log"
  echo "     output file name for log messages (text file)"
  echo "  -fname_color <fname_color> default = $fname_color"
  echo "     name of text file with RGB color codes for each fiber"
  echo "  -ulay <fstem>              default = $ulay"
  echo "     fstem for underlay image (T1, FA, or none)"
  echo "  -ulay_plane <plane>        default = $ulay_plane"
  echo "     plane for underlay image (sag, hor, or cor)"
  echo "  -view <view>               default = $view"
  echo "     image view (L, R, A, P, I, S)"
  echo "  -xrot <xrot>             default = $xrot"
  echo "     [0|1] rotation about x axis"
  echo "  -yrot <yrot>             default = $yrot"
  echo "     [0|1] rotation about y axis"
  echo "  -zrot <zrot>             default = $zrot"
  echo "     [0|1] rotation about z axis"
  echo "  -xflip <xflip>             default = $xflip"
  echo "     [0|1] whether to flip streamlines in x direction"
  echo "  -yflip <yflip>             default = $yflip"
  echo "     [0|1] whether to flip streamlines in y direction"
  echo "  -zflip <zflip>             default = $zflip"
  echo "     [0|1] whether to flip streamlines in z direction"
  echo "  -image_alpha <image_alpha> default = $image_alpha"
  echo "     transparency of image"
  echo "  -fiber_alpha <fiber_alpha> default = $fiber_alpha"
  echo "     transparency of fiber"
  echo " "

exit 1;
###############################################################################

###############################################################################
find_underlay:

  # check underlay file
  if ($ulay == 'none') then
    set ulay_flag = 0
## todo: need to set fname_image to something that exists (like an empty volume?)
  else if($ulay == 'T1') then
    set ulay_flag = 1
    set fname_image = $indir'/'$ulay'_resDTI.mgh'
    if (! -e $fname_image) then
      set fname_image = $indir'/'$ulay'.mgh'
    endif
  else if($ulay == 'FA') then
    set ulay_flag = 1
    set fname_image = $indir'/'$ulay'.mgh'
  else
    echo "ERROR: invalid underlay type $ulay"
    exit 1;
  endif

  if (! -e $fname_image) then
    echo "ERROR: underlay image $fname_image not found"
    exit 1;
  endif

  if ($ulay_flag == 1) then
    if ($ulay_plane == 'sag') then
      set ulay_parms = '-norendhor -norendcor'
    else if ($ulay_plane == 'cor') then
      set ulay_parms = '-norendhor -norendsag'
    else if ($ulay_plane == 'hor') then
      set ulay_parms = '-norendsag -norendcor'
    endif    
  else
    set ulay_parms = '-norendsag -norendhor -norendcor'
  endif

goto find_underlay_return;
###############################################################################

###############################################################################
find_fibers:

  set fiberdir = $indir'/fiber_paths'
  if (! -e $fiberdir) then
    echo "ERROR: fiber paths dir $fiberdir not found"
    exit 1;
  endif

  # create list of fiber files
  set fname_fiber_list = ()
  foreach fnum ($fnums)
    set fname_fiber = $fiberdir/'fiber_'$fnum$fiber_suffix
    if (! -e $fname_fiber) then
      echo "WARNING: file $fname_fiber not found"
      continue
    endif
    set fname_fiber_list = ($fname_fiber_list $fname_fiber)
  end

  if ($#fname_fiber_list == 0) then
    echo "ERROR: no fiber files matching suffix $fiber_suffix in $fiberdir $fiberdir"
    exit 1;
  endif

goto find_fibers_return;
###############################################################################

###############################################################################
run_tractoview:

  # run tractoview
  echo "running tractoview..."
  if (! -e $fname_tif || $forceflag) then
    tractoview \
      -fnameimage $fname_image \
      -fnamefibers $fname_fiber_list \
      -initview $view \
      -xrot $xrot -yrot $yrot -zrot $zrot \
      -image_alpha $image_alpha \
      -fiber_alpha $fiber_alpha \
      $ulay_parms $color_parms $flip_parms $save_parms | tee $fname_log
  endif

goto run_tractoview_return;
###############################################################################

