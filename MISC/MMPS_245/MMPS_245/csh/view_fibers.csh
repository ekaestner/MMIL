#!/bin/csh -f
#
# view_fibers.csh -- wrapper for tractoview
#
# created:  05/18/12 by Don Hagler
# last mod: 05/23/12 by Don Hagler
#
###############################################################################

## initialize parameters default values ##
goto init_parms;
init_parms_return:

if($#argv < 2) then
  goto usage_exit;
  exit 1;
endif

## parse input arguments ##
goto parse_args;
parse_args_return:

## check parameters ##
goto check_parms;
check_parms_return:

## run tractoview ##
goto run_tractoview;
run_tractoview_return:

exit 0;

###############################################################################

###############################################################################
init_parms:

  set fname_image = ()
  set fname_fiber = ()
  set save_flag = 0
  set fname_tif = 'fibers.tif'
  set fname_log = 'fibers.log'
  set fname_color = 'color.dat'
  set ulay_flag = 1
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

  # todo: make these work:
  set fname_ovly = ()
  set ovly_offset = 0
  set ovly_slope = 1
  set ovly_thresh = 0
  set fname_alpha = ()
  set min_fiber_alpha = 0.1
  set rendcor = 1
  set rendsag = 1
  set rendhor = 1
  set forceflag = 0

goto init_parms_return;
###############################################################################

###############################################################################
parse_args:

  set cmdline = ($argv);
  set fname_image = $argv[1]; shift;
  set fname_fiber = $argv[1]; shift;

  set fname_fiber_list = ($fname_fiber)

  while( $#argv != 0 )
    set flag = $argv[1]; shift;
    switch($flag)

      case "-fname_fiber":
        if ( $#argv == 0) goto arg1err;
        set fname_fiber = $argv[1]; shift;
        set fname_fiber_list = ($fname_fiber_list $fname_fiber)
        breaksw

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

      case "-ulay_flag":
        if ( $#argv == 0) goto arg1err;
        set ulay_flag = $argv[1]; shift;
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

  if ($#fname_image == 0) then
    goto usage_exit;
  endif

  if (! -e $fname_image) then
    echo "ERROR: underlay image $fname_image not found"
    exit 1;
  endif

  foreach fname_fiber ($fname_fiber_list)
    if (! -e $fname_fiber) then
      echo "ERROR: file $fname_fiber not found"
      exit 1;
    endif
  end

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

  if (! -e $fname_color) then
    echo "WARNING: color file $fname_color not found"
    set color_parms = ()
  else
    set color_parms = "-fname_color $fname_color"
  endif

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
  echo "USAGE: view_fibers.csh <fname_image> <fname_fiber> [options]"
  echo " "
  echo "Required Arguments:"
  echo "  fname_image : name of file containing underlay volume (mgh format)"
  echo "                NOTE: must be mgh, not mgz"
  echo "  fname_fiber : name of DTI Studio format fiber file (grp format)"
  echo " "
  echo "Optional Arguments:"
  echo "  -fname_fiber <fname_fiber>"
  echo "      specify additional fiber files to view simultaneously"
  echo "  -save_flag <save_flag>     default = $save_flag"
  echo "     [0|1] whether to save rendered image (otherwise display only)"
  echo "  -fname_tif <fname_tif>     default = $fname_tif"
  echo "     output file name for rendered image (tif format)"
  echo "  -fname_log <fname_log>     default = $fname_log"
  echo "     output file name for log messages (text file)"
  echo "  -fname_color <fname_color> default = $fname_color"
  echo "     name of text file with RGB color codes for each fiber"
  echo "  -ulay_flag <ulay_flag>              default = $ulay_flag"
  echo "     [0|1] whether to display underlay image"
  echo "  -ulay_plane <plane>        default = $ulay_plane"
  echo "     plane for underlay image (sag, hor, or cor)"
  echo "  -view <view>               default = $view"
  echo "     image view (L, R, A, P, I, S)"
  echo "  -xrot <xrot>             default = $xrot"
  echo "     rotation about x axis"
  echo "  -yrot <yrot>             default = $yrot"
  echo "     rotation about y axis"
  echo "  -zrot <zrot>             default = $zrot"
  echo "     rotation about z axis"
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
      $ulay_parms $flip_parms $color_parms $save_parms | tee $fname_log
  endif

goto run_tractoview_return;
###############################################################################

