#!/bin/csh -f

##################################################
#
# ts_make_movie -- make movie from freesurfer time-course files
#
#    currently only supports mgh files
#
#  edit log:
#    bcipolli on 2007/06/13
#
###################################################

set progname = ts_make_movie
goto mainprog;


####################################################

initvars:
  ##############################
  # constants

  #  Program names
  
  # Others
  set valid_views = ( lat ven med pos dor cus );
  set valid_monts = ( custom onesub comparesubs groupsummary );
  set space_replacement = llcoolj
  
  ##############################
  #          One-value variables
  set rmintermedflag = 0
  set forceflag = 0
  set remakeonlyflag = 0
  
  set flatflag = 0
  set scalebarflag = 0
  set colscalebarflag = 0
  set polarflag = 0
  set complexflag = 0
  
  set outdir = `pwd`
  set outstem = $progname
  
  set mont = custom
  set nimages = 1
  set nrows = 1
  set ncols = 1
  @ ncells = $nrows * $ncols
  set geometry = 250x250
  
  set bgcolor = \#000000
  set fgcolor = \#ffffff
  set transparent = \#000000

  set title = ""  
  set image_type = "gif";
  
    
  ##############################
  #          Array-based variables
  if ( $?SUBJECTS_DIR ) then
    set sdirs = ( $SUBJECTS_DIR );
  else
    set sdirs = ();
  endif

  set srcfiles = ()
  set checksums = (); # this is computed;
  set subjs = ()
  set hemis = ( rh )
  set surfs = ( inflated )
  set patchs = ( full )
  set fthreshs = ( 0 )
  set fslopes = ( 1 )
  set fmids = ( 1 )
  set sparsesmooths = ( 20 )
  set postsmooths = ( 5 )
  set scales = ( 1.0 )
  set views = ( lat )
  set rotxs = ( 0 )
  set rotys = ( 0 )
  set rotzs = ( 0 )
  set labels = ( " " )
  set offsets = ( 0.2 )
  set cvfacts = ( 1.5 )
  set indirs = ( `pwd` )
  set tfirsts = ( 0 )
  set tlasts = ( -1 )
  set tmaxs = ( -1 )
  set ndigits = ( 3 )
  set positions = ( 1 )
  
  
  goto initvars_return;


####################################################

checkdeps:

  #	Make sure variables we expect to be defined ARE defined
  if ( ! $?IMAGEMAGICK_DIR) goto no_movie_setup;
  if ( ! $?TKSURFER_PROG)   goto no_movie_setup;
  
  #	Just because things are set up doesn't mean
  #	the binaries we want exist.  Do a cursory scan
  #	to make sure all our dependencies are available.
  if ( ! -e "$IMAGEMAGICK_DIR/montage" ) goto no_montage
  if ( ! -e "$IMAGEMAGICK_DIR/convert" ) goto no_convert
  
  if ( "`which fs_surfmgh.csh | grep /`" == "") goto no_fs_surfmgh;
  
  if ( "`which $TKSURFER_PROG | grep /`" == "" ) goto no_tksurfer;
  if ( "`which mri_info | grep /`" == "") goto no_mri_info;

  goto checkdeps_return;

  
####################################################
#
#	Define default array values for a named montage
#
set_montage_values:

  switch ($mont)
    case "custom"
    breaksw
    
    case "onesub"
    breaksw
    
    default
      echo "Unknown montage: $mont";
      exit 10;
    
  endsw
  
  goto set_montage_values_return;
  
  
####################################################

parse_args:
  if ($#argv < 4) goto nea1err;

  set cmdline=($argv);
  
  while( $#argv != 0 )

    set flag = $argv[1]; 
    shift;
    
    switch($flag)
      
      #          Simple flags (no args)
#      case "-polar"
#        set polarflag = 1;
#        set complexflag = 1;
#        breaksw

      case "-nocache"
      case "-rmintermed"
        set rmintermedflag = 1;
        breaksw

      case "-remakeonly"
        set remakeonlyflag = 1;
        breaksw

      case "-force"
        set forceflag = 1;
        breaksw

      
      #          One-value flags
      case "-outdir"
        if ( $#argv == 0) goto arg1err;
        set outdir = $argv[1]; shift;
        breaksw

      
      case "-outstem"
        if ( $#argv == 0) goto arg1err;
        set outstem = $argv[1]; shift;
        breaksw

      case "-mont"
        if ( $#argv == 0) goto arg1err;
        set mont = $argv[1]; shift;
	
        # set the defaults for the specified montage.
        # by doing that now (e.g. as early as possible), then
        # any flags after this one will override the "default"
        # settings for the montage, which is a very powerful thing.
        goto set_montage_values;
        set_montage_values_return:
        breaksw

      case "-nrows"
        if ( $#argv == 0) goto arg1err;
        set nrows = $argv[1]; shift;
        @ ncells = $nrows * $ncols;
        breaksw

      case "-ncols"
        if ( $#argv == 0) goto arg1err;
        set ncols = $argv[1]; shift;
        @ ncells = $nrows * $ncols;
        breaksw

      case "-geometry"
        if ( $#argv == 0) goto arg1err;
        set geometry = $argv[1]; shift;
        breaksw

      case "-bgcolor"
        if ( $#argv == 0) goto arg1err;
        set bgcolor = $argv[1]; shift;
        breaksw


      #   Array-based args
      case "-sdir"
        set sdirs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set sdirs = ( ${sdirs[*]} $argv[1] )
          shift
        end
        if ( $#sdirs == 0 ) goto arg1err;
        if ( $#sdirs > $nimages ) set nimages = $#sdirs;
        breaksw

      
      case "-subj"
        set subjs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set subjs = ( ${subjs[*]} $argv[1] )
          shift
        end
        if ( $#subjs == 0 ) goto arg1err;
	      if ( $#subjs > $nimages ) set nimages = $#subjs;
        breaksw

	case "-hemi"
        set hemis = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set hemis = ( ${hemis[*]} $argv[1] )
          shift
        end
        if ( $#hemis == 0 ) goto arg1err;
		if ( $#hemis > $nimages ) set nimages = $#hemis;
        breaksw
          
      case "-srcfile"
        set srcfiles = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set srcfiles = ( ${srcfiles[*]} $argv[1] )
          shift
        end
        if ( $#srcfiles == 0 ) goto arg1err;
		if ( $#srcfiles > $nimages ) set nimages = $#srcfiles;
        breaksw

      case "-tfirst"
        set tfirsts = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set tfirsts = ( ${tfirsts[*]} $argv[1] )
          shift
	    end
	    if ( $#tfirsts == 0 ) goto arg1err;
		if ( $#tfirsts > $nimages ) set nimages = $#tfirsts;
        breaksw

      case "-tlast"
        set tlasts = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set tlasts = ( ${tlasts[*]} $argv[1] )
          shift
        end
        if ( $#tlasts == 0 ) goto arg1err;
		if ( $#tlasts > $nimages ) set nimages = $#tlasts;
        breaksw

      case "-tmax"
        # push this flag back onto the flag stack
        # with it's actual contentful associate
        #set argv = { '-t
        set tmaxs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          
          set tmaxs = ( ${tmaxs[*]} $argv[1] )
          shift
        end
        if ( $#tmaxs == 0 ) goto arg1err;
		if ( $#tmaxs > $nimages ) set nimages = $#tmaxs;
        breaksw

      case "-surf":
        set surfs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set surfs = ( ${surfs[*]} $argv[1] )
          shift
        end
        if ( $#surfs == 0 ) goto arg1err;
		if ( $#surfs > $nimages ) set nimages = $#surfs;
        breaksw

      case "-fthresh"
        set fthreshs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set fthreshs = ( ${fthreshs[*]} $argv[1] )
          shift
        end
        if ( $#fthreshs == 0 ) goto arg1err;
		if ( $#fthreshs > $nimages ) set nimages = $#fthreshs;
        breaksw

      case "-fslope"
        set fslopes = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set fslopes = ( ${fslopes[*]} $argv[1] )
          shift
        end
        if ( $#fslopes == 0 ) goto arg1err;
		if ( $#fslopes > $nimages ) set nimages = $#fslopes;
        breaksw

      case "-fmid"
        set fmids = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set fmids = ( ${fmids[*]} $argv[1] )
          shift
        end
        if ( $#fmids == 0 ) goto arg1err;
		if ( $#fmids > $nimages ) set nimages = $#fmids;
        breaksw

      case "-postsmooth"
        set postsmooths = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set postsmooths = ( ${postsmooths[*]} $argv[1] )
          shift
        end
        if ( $#postsmooths == 0 ) goto arg1err;
		if ( $#postsmooths > $nimages ) set nimages = $#postsmooths;
        breaksw

      case "-sparsesmooth"
        set sparsesmooths = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set sparsesmooths = ( ${sparsesmooths[*]} $argv[1] )
          shift
        end
        if ( $#sparsesmooths == 0 ) goto arg1err;
		if ( $#sparsesmooths > $nimages ) set nimages = $#sparsesmooths;
        breaksw

      case "-offset"
        set offsets = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set offsets = ( ${offsets[*]} $argv[1] )
          shift
        end
        if ( $#offsets == 0 ) goto arg1err;
		if ( $#offsets > $nimages ) set nimages = $#offsets;
        breaksw

      case "-cvfact"
        set cvfacts = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set cvfacts = ( ${cvfacts[*]} $argv[1] )
          shift
        end
        if ( $#cvfacts == 0 ) goto arg1err;
		if ( $#cvfacts > $nimages ) set nimages = $#cvfacts;
        breaksw

      case "-scale"
        set scales = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set scales = ( ${scales[*]} $argv[1] )
          shift
        end
        if ( $#scales == 0 ) goto arg1err;
		if ( $#scales > $nimages ) set nimages = $#scales;
        breaksw

      case "-view"
        set views = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set views = ( ${views[*]} $argv[1] )
          shift
        end
        if ( $#views == 0 ) goto arg1err;
		if ( $#views > $nimages ) set nimages = $#views;
        breaksw

      case "-rotx"
        set rotxs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set rotxs = ( ${rotxs[*]} $argv[1] )
          shift
        end
        if ( $#rotxs == 0 ) goto arg1err;
		if ( $#rotxs > $nimages ) set nimages = $#rotxs;
        breaksw

      case "-roty"
        set rotys = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set rotys = ( ${rotys[*]} $argv[1] )
          shift
        end
        if ( $#rotys == 0 ) goto arg1err;
 		if ( $#rotys > $nimages ) set nimages = $#rotys;
       breaksw

      case "-rotz"
        set rotzs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set rotzs = ( ${rotzs[*]} $argv[1] )
          shift
        end
        if ( $#rotzs == 0 ) goto arg1err;
		if ( $#rotzs > $nimages ) set nimages = $#rotzs;
        breaksw

      case "-label"
        set labels = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;
          set labels = ( ${labels[*]} `echo $argv[1] | sed "s/\s/$space_replacement/gi"` )
          shift
        end
        if ( $#labels == 0 ) goto arg1err;
        if ( $#labels > $nimages ) set nimages = $#labels;
        breaksw

      case "-patch"
        set patchs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set patchs = ( ${patchs[*]} $argv[1] )
          shift
        end
        if ( $#patchs == 0 ) goto arg1err;
		if ( $#patchs > $nimages ) set nimages = $#patchs;
        breaksw

      case "-indir"
        set indirs = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set indirs = ( ${indirs[*]} $argv[1] )
          shift
        end
        if ( $#indirs == 0 ) goto arg1err;
		if ( $#indirs > $nimages ) set nimages = $#indirs;
        breaksw

      case "-position"
        set positions = ()
        while ( $#argv != 0  )
          if (`echo $argv[1] | grep ^-` != "") break;

          set positions = ( ${positions[*]} $argv[1] )
          shift
        end
        if ( $#positions == 0 ) goto arg1err;
		if ( $#positions > $nimages ) set nimages = $#positions;
        breaksw


#      case "-flat"
#        set flatflag = 1;
#        breaksw

#      case "-scalebar"
#        set scalebarflag = 1
#        breaksw

#      case "-colscalebar"
#        set colscalebarflag = 1
#        breaksw

#      case "-colscalebar"
#        set colscalebarflag = 1
#        breaksw

      
      # error case
      default:
        echo ERROR: Flag $flag unrecognized. 
        echo $cmdline
        exit 1
        breaksw
    endsw
  end
  
  
  goto parse_args_return;


#############################################################

check_params:


  ############################################
  #  Validate required parameters
  
  #
  if ($#subjs == 0) then
    echo "ERROR: subject name not set!"
    goto usage_exit;
  endif

  if ($#srcfiles == 0) then
    echo "ERROR: srcfile not set!"
    goto usage_exit;
  endif

  # check montages
  set found = "";
  foreach goodmont ($valid_monts)
    if ($goodmont == $mont) set found = $goodmont;
  end
  if ($found == "") goto badmont1err;
  
  # check subjects dir
  if ($#sdirs == 0) then
    echo "NO SUBJECTS_DIR set, no -sdir parameter passed.";
    goto usage_exit;
  endif

  #############################################
  #  Validate optional parameters
  
  #  check hemis
  foreach hemi ($hemis)
    if ("$hemi" != "rh" && "$hemi" != "lh") then
      echo "ERROR: hemi must be rh or lh"
      exit 1;
    endif
  end
  
  # check views
  foreach view ($views)
    set found = ""
    foreach goodview ($valid_views)
      if ($goodview == $view) set found = $goodview;
    end;

    if ($found == "") goto badview1err;
  end
  


  ########################################
  #
  #	Validate length of array-based params.
  set allarrs = ( "sdirs" "srcfile" "subj" "hemi" "surf" "patch" "fthresh" "fslope" "fmid" \
                  "sparsesmooth" "postsmooth" "scale" "view" "rotx" "roty" "rotz" "label" \
				  "offset" "cvfact" "indir" "tfirst" "tlast" "tmax" "ndigit" "position" );
  set allarrlens = ( $#sdirs $#srcfiles $#subjs $#hemis $#surfs $#patchs $#fthreshs $#fslopes $#fmids \
                  $#sparsesmooths $#postsmooths $#scales $#views $#rotxs $#rotys $#rotzs $#labels \
				  $#offsets $#cvfacts $#indirs $#tfirsts $#tlasts $#tmaxs $#ndigits $#positions );
  set badarrs = ( );
  set maxdarrs = ( );
  set i = 1
  while ($i <= $#allarrs)
    if ($allarrlens[$i] != 1 && $allarrlens[$i] != $nimages) then
      set badarrs = ( ${badarrs[*]} $allarrs[$i] )
    else if ($allarrlens[$i] == $nimages) then
      set maxdarrs = ( ${maxdarrs[*]} $allarrs[$i] )
    endif
	
    @ i++
  end

  #  We found args with inconsistent #s of values, so
  #  throw an error message
  if ($#badarrs != 0) goto arg2err


  ########################################
  #
  #  Make sure the # of rows & columns can 
  #  support the total # of images displayed.
  
  #	Calculate a grid if no grid was specified.
  if ( $ncells == 1 ) then
     set nrows = `echo "scale=0;sqrt($nimages)" | bc -l`;
     set rem   = `echo "scale=0;$nimages%$nrows" | bc -l`;
     set ncols = `echo "scale=0;$nimages/$nrows" | bc -l`;
     
     if ($rem == 1) set ncols = `echo "scale=0;$ncols+1" | bc -l`;
     @ ncells = $nrows * $ncols;
  endif

  if ( $ncells < $nimages )  goto badgrid1err;


  ########################################
  #
  #  Make sure all the grid positions 
  #  will fit within the current grid
  #
  foreach position ($positions)
    if ( $position > $ncells ) then
      goto badgrid2err;
    endif
  end
  

  goto check_params_return;


#############################################################

nea1err: 
  echo "ERROR: required params not specified";
  goto usage_exit;

arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1
 
arg2err:
  echo "ERROR: different flags have different #s of args:"
  echo "     BAD ( 1<args<max[$nimages] ): ${badarrs[*]}";
  echo "     OK: (   args=max[$nimages] ): ${maxdarrs[*]}";
  echo " ";
  echo "     ** each flag must have either 1 value or (# images) values **"
  echo " ";
  exit 2
  
badgrid1err:
  echo "ERROR: $nrows x $ncols grid cannot support $nimages images.";
  echo "       please increase your grid size, or remove your nrows / ncols parameters"
  echo "       to allow $progname to calculate a grid for you.";
  exit 3
  
badgrid2err:
  echo "ERROR: a $nrows x $ncols grid cannot support the following grid positions:"
  foreach position ($positions)
    if ($position > $ncells) then
      echo "       position $position";
    endif
  end
  echo " ";
  exit 4
  
badsrc1err:
  echo "ERROR: source file(s) could not be found:"
  set i = 1
  while ($i <= $nimages)
    if (! -e "$indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh") then
      echo "       $indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh";
    endif
    @ i++
  end
  exit 5
  
badframes1err:
  echo "ERROR: # of frames is not consistent across images:"
  set i = 1
  while ($i <= $nimages)
    @ nframes = 1 + $tlasts[$i] - $tfirsts[$i];
    echo "       image ${i}:  $nframes frames";
    
    @ i++
  end
  exit 6
  
nosdir1err:
  echo "ERROR: SUBJECTS_DIR not set and no 'sdir' flag... quitting"
  exit 7
  
badsdir1err:
  echo "ERROR: SUBJECTS_DIR not found: "
  set i = 1
  while ($i <= $nimages)
    if ( -e "$sdirs[$i]" ) continue;
    echo "       sdir ${i}:  $sdirs[$i]";
    @ i++
  end
  exit 7
  

badview1err:
  echo "ERROR: all views must be one of [$valid_views[*]]"
  foreach view ($views)
    set found = ""
    foreach goodview ($valid_views)
      if ($goodview == $view) set found = $goodview;
    end;

    if ($found == "") echo "  ** view = $view is invalid!";
    
  end
  exit 8

badmont1err:
  echo "ERROR: montage must be one of [$valid_monts[*]]"
  echo "       specified montage: [$mont]";
  exit 9;
  


no_movie_setup:
  echo "ERROR: looks like you didn't run the SetUpMovieMaker.csh script"
  echo "  before running this program.";
  exit 21;
  
no_fs_surfmgh:
  echo "ERROR: fs_surfmgh.csh was not found in you path, "
  echo "       nor was the common path found ($fs_surfmgh_path)."
  echo "       Please find the directory containing fs_surfmgh.csh and add it to your path."
  exit 10;

no_montage:
no_convert:
no_imagemagik:
  echo "ERROR: some imageMagick tools were not found in you path;"
  echo "       this program requires BOTH the montage AND convert tools to run."
  echo "       Please find the directory containing these tools and add it to your path."
  exit 11;

no_mri_info:
  echo "ERROR: couldn't locate mri_info."
  echo " ";
  echo "You need to resolve this issue before running the movie renderer.";
  exit 15;
  
no_tksurfer
  echo "ERROR: couldn't locate any version of tksurfer; both not in your "
  echo "       path and also not in the expected mmildev shared path.";
  echo " ";
  echo "You need to resolve this issue before running the movie renderer.";
  exit 16;
  
#############################################################
#
# prepare_vars:
#
#    massage internal variables before main program execution,
#    to allow main program logic to stand out more clearly
#

prepare_vars:

  
  ########################################
  #
  #  
  #  At this point, we know that all arrays have either one value or 
  #  (nimages) values.  Simply make all arrays have (nimages) values
  #  so that we don't have to worry about it in the future.
  #
  set i = 2
  set p = 1
  while ($i <= $nimages)
    if ($#sdirs < $i)    set sdirs = ( ${sdirs[*]} ${sdirs[$p]} )
    if ($#srcfiles < $i) set srcfiles = ( ${srcfiles[*]} ${srcfiles[$p]} )
    if ($#subjs < $i)    set subjs = ( ${subjs[*]} ${subjs[$p]} )
    if ($#hemis < $i)    set hemis = ( ${hemis[*]} ${hemis[$p]} )
    if ($#surfs < $i)    set surfs = ( ${surfs[*]} ${surfs[$p]} )
    if ($#patchs < $i)   set patchs = ( ${patchs[*]} ${patchs[$p]} )
    if ($#fthreshs < $i) set fthreshs = ( ${fthreshs[*]} ${fthreshs[$p]} )
    if ($#fslopes < $i)  set fslopes = ( ${fslopes[*]} ${fslopes[$p]} )
    if ($#fmids < $i)    set fmids = ( ${fmids[*]} ${fmids[$p]} )
    if ($#sparsesmooths < $i) set sparsesmooths = ( ${sparsesmooths[*]} ${sparsesmooths[$p]} )
    if ($#postsmooths < $i)   set postsmooths = ( ${postsmooths[*]} ${postsmooths[$p]} )
    if ($#scales < $i)   set scales = ( ${scales[*]} ${scales[$p]} )
    if ($#views < $i)    set views = ( ${views[*]} ${views[$p]} )
    if ($#rotxs < $i)    set rotxs = ( ${rotxs[*]} ${rotxs[$p]} )
    if ($#rotys < $i)    set rotys = ( ${rotys[*]} ${rotys[$p]} )
    if ($#rotzs < $i)    set rotzs = ( ${rotzs[*]} ${rotzs[$p]} )
    if ($#labels < $i)   set labels = ( ${labels[*]} ${labels[$p]} )
    if ($#offsets < $i)  set offsets = ( ${offsets[*]} ${offsets[$p]} )
    if ($#cvfacts < $i)  set cvfacts = ( ${cvfacts[*]} ${cvfacts[$p]} )
    if ($#indirs < $i)   set indirs = ( ${indirs[*]} ${indirs[$p]} )
    if ($#tfirsts < $i)  set tfirsts = ( ${tfirsts[*]} ${tfirsts[$p]} )
    if ($#tlasts < $i)   set tlasts = ( ${tlasts[*]} ${tlasts[$p]} )
    if ($#tmaxs < $i)    set tmaxs = ( ${tmaxs[*]} ${tmaxs[$p]} )
    if ($#ndigits < $i)  set ndigits = ( ${ndigits[*]} ${ndigits[$p]} )

    if ($#positions < $i) set positions = ( ${positions[*]} $i )

    @ i++
    @ p++
  end
  

  #######################################
  #  Validate the input subjects
  #

  set i = 1
  while ($i <= $nimages)

    if ("$sdirs[$i]" == "") then
      goto nosdir1err;
    else if ( ! -e "$sdirs[$i]" ) then
      goto badsdir1err;
    else if (! -e "$sdirs[$i]/$subjs[$i]") then
      echo "$progname : ### Subject not found... quitting: $sdirs[$i]/$subjs[$i]"
      exit 1
    endif
    
    @ i++
  end

  ########################################
  #
  #  For input files that don't have one boundary
  #	 of frames set, then read the boundary from disk
  #
  #  Don't forget to limit this auto-frame-size reading
  #  by tmax
  #
  set i = 1
  while ($i <= $nimages)
  
    # read the total # of frames from disk  
    if ($tlasts[$i] == -1) then
      set nframes = `mri_info $indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh --nframes`;
      
      if ("$nframes" == "") then
        if ( ! -e "$indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh" ) then
	  echo " " ;
	  echo "Source file does not exist: $indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh";
	  exit 3;
	  
	else
	  echo " " ;
	  echo "Exiting due to possibly corrupt source file: $indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh";
	  echo " ";
	  exit 3;
	endif
      endif

      @ tlasts[$i] = $nframes - 1;
    endif
    

    # limit the last frame by tmax
    # (-1 represents 'no limit')
    if ($tlasts[$i] > $tmaxs[$i] && $tmaxs[$i] != -1) then
      set tlasts[$i] = $tmaxs[$i]
    endif
 
    @ i++
  end  



  ########################################
  #
  #  Compute the number of digits needed 
  #  to fully represent all timepoints.  This
  #  allows us to pad timepoints so that ouput
  #  filenames sorted by string value will also
  #  be sorted by timepoint.
  #
  set ndigits = ()
  foreach tlast ($tlasts)
    if ($tlast > 9999) then
      set ndigits = (${ndigits[*]} 5 )
    else if ($tlast > 999) then
      set ndigits = (${ndigits[*]} 4 )
    else if ($tlast > 99) then
      set ndigits = (${ndigits[*]} 3 )
    else if ($tlast > 9) then
      set ndigits = (${ndigits[*]} 2 )
    else
      set ndigits = (${ndigits[*]} 1 )
    endif
  end    

    
  ########################################
  #
  #  Make sure that all images in the movie
  #  have the same total # of frames
  set i = 2
  @ nframes = 1 + $tlasts[1] - $tfirsts[1]
  while ($i <= $nimages)
    @ nframes2 = 1 + $tlasts[$i] - $tfirsts[$i];
    
    if ( "$nframes" != "$nframes2" ) then
      goto badframes1err;
    endif
    
    @ i++
  end


  #################################
  #
  #	Finally, we're ready to generate 
  #     input and output files, and validate
  #     them.
  #  
  set i = 1
  echo "Checking src files and computing checksums...";
  while ($i <= $nimages)
    if ( ! -e $indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh ) then
      goto badsrc1err;
    else
      #set cs = `md5sum "$indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh" | sed -r "s/\s.*//g"`;
      set cs = `date -r "$indirs[$i]/$srcfiles[$i]-$hemis[$i].mgh" "+%Y%m%d-%k%M%S%N"`;
      set checksums = ( ${checksums[*]} $cs );
    endif

    echo "  + computed checksum [$i] of [$nimages] ($checksums[$i])"
    @ i++
  end

  goto prepare_vars_return;
  
  
#############################################################

generatempg:

  #	Set the subjects dir;
  #	WE MUST REMEMBER TO RESET!
  if ( $?SUBJECTS_DIR ) set old_subj_dir = "$SUBJECTS_DIR";
  
  ############################
  #
  # IMAGE-GENERATION LOOP
  #
  # we're going to create all necessary
  # images, calling fs_surfmgh only if/when necessary.
  echo "Searching for image files that need creation..."
  
  #	Loop over each image set
  set i = 1
  while ($i <= $nimages)
    echo "rendering image [$i] of [$nimages], frames [$tfirsts[$i]] - [$tlasts[$i]]";
    
    set rgbdir = `echo $srcfiles[$i] | sed "s/[^\/]*\///g"`
    set rgbdir = $outdir/rgb/$rgbdir;
    set image_stem = $subjs[$i]-$checksums[$i]

    mkdir -p $rgbdir
      
    set rargs = ""
    if ($forceflag)       set rargs = "$rargs -force"
    if ($flatflag)        set rargs = "$rargs -flat -patch $patchs[$i]"
    if ($scalebarflag)    set rargs = "$rargs -scalebar"
    if ($colscalebarflag) set rargs = "$rargs -colscalebar"
    if ($polarflag)       set rargs = "$rargs -polar"

    setenv SUBJECTS_DIR $sdirs[$i];

    fs_surfmgh.csh $subjs[$i] \
		$srcfiles[$i] \
		$hemis[$i] \
		-surf $surfs[$i] \
		-tfirst $tfirsts[$i] \
		-tlast $tlasts[$i] \
		-fthresh $fthreshs[$i] \
		-fmid $fmids[$i] \
		-offset $offsets[$i] \
		-cvfact $cvfacts[$i] \
		-scale $scales[$i] \
		-view $views[$i] \
		-rotx $rotxs[$i] \
		-roty $rotys[$i] \
		-rotz $rotzs[$i] \
		-smooth $postsmooths[$i] \
		-outstem $image_stem \
		-outdir $rgbdir \
		-indir $indirs[$i] \
		$rargs \
		-rmintermed \
		-savergb \
		-offscreen \
		-silent

    if ($?old_subj_dir) setenv SUBJECTS_DIR $old_subj_dir
			
    @ i++
  end
    
      
  ############################
  # MONTAGE LOOP
  #  loop over each timepoint.
  #   - montage the files into a single frame
 
 
  #  Create an output file
  set datetime = `date "+%F.%H-%M-%S"`;
  set outfile_mpg = $outstem.$datetime.$geometry.mpg
  if ( -e "$outdir/$outfile_mpg" ) rm -f "$outdir/$outfile_mpg";

  set montaged_files = ();
  @ nslices = 1 + $tlasts[1] - $tfirsts[1]
  
  echo "Creating $nslices montages ...";
  set montagedir = $outdir/montages
  mkdir -p $montagedir
  
  set n = 1
  while ($n <= $nslices)

    set npadded    = `count -dig $ndigits[1] $n $n`
    set montaged_files = ( ${montaged_files[*]} $montagedir/$outstem-$geometry-$npadded.$image_type );

    echo "Montage [$n] of [$nslices]...";
    
    if ( ! $remakeonlyflag || ! -e ${montaged_files[$#montaged_files]} ) then

    
      #######################################
      #	Loop over each image set
      # to add the image to the list that needs to be
      # montaged
      set image_files = ()
      set i = 1
      while ($i <= $nimages)
        @ t = $n + $tfirsts[$i] - 1;

        set rgbdir = `echo $srcfiles[$i] | sed "s/[^\/]*\///g"`
        set rgbdir = $outdir/rgb/$rgbdir;
      
        set tpadded    = `count -dig $ndigits[$i] $t $t`
        set image_stem = $subjs[$i]-$checksums[$i]
        set image_file = $image_stem-$tpadded-$hemis[$i]-$surfs[$i]-$views[$i].rgb;
  
        set image_files = ( ${image_files[*]} $rgbdir/$image_file );
      
        @ i++
      end


      ####################################   
      #	Now we have a set of $nimages images and a grid of $positions.
      #   We need to reorder images to follow the $positions grid, and to
      #   insert blank images for 'skipped' grid positions.
      set gridded_images = ();

      # loop through every grid position
      set g = 1
      while ($g <= $ncells)

        # search for an image specified for this grid location
        set im = 1
        while ($im <= $#positions)
          if ($positions[$im] == $g) break;
          @ im++
        end
      
        # no image specified for this grid location, so use "blank"
        if ($im <= $#positions) then
	  set gridded_images = ( ${gridded_images[*]} \
	                         -label $labels[$im] \
				 -background $bgcolor \
				 -geometry $geometry \
				 -stroke $fgcolor \
				 -transparent $transparent \
				 $image_files[$im] );
        else
          if ( ! -e "$outdir/montages/_blank.gif" ) then
	    echo "" | "$IMAGEMAGICK_DIR/convert" -background $bgcolor -page 1x1 text:- $outdir/montages/_blank.gif
	  endif
          set gridded_images = ( ${gridded_images[*]} +label $outdir/montages/_blank.gif );
        endif
       
        @ g++
      end
  
      # | sed "s/-label /-label \\'/gi" | sed "s/ -background/\\' -background/gi"
      # Make the montage
      "$IMAGEMAGICK_DIR/montage" -tile ${ncols}x${nrows} \
	      `echo "${gridded_images[*]}" | sed "s/$space_replacement/_/gi"` \
	      ${montaged_files[$#montaged_files]};
    endif

    # Add it to the mpg
    #if ($n == 1) then
    #  convert ${montaged_files[$#montaged_files]} "mpg:$outdir/$outfile_mpg"
    #else
    #  convert -adjoin "$outdir/$outfile_mpg" ${montaged_files[$#montaged_files]} mpg:$outdir/$outfile_mpg
    #endif
    
    @ n++
  end


  # make mpg
  if ( ! -e "$outdir/$outfile_mpg" ) then
    echo "$progname : concatenating [$nslices] images to create mpg..."
    "$IMAGEMAGICK_DIR/convert" -adjoin ${montaged_files[*]} "mpg:$outdir/$outfile_mpg"
  endif
				
  # remove intermediates
  if ( $rmintermedflag ) then
    echo "$progname : removing intermediate files"
    rm -rf $outdir/rgb
    rm -rf $outdir/montages
  endif

  
  #	Set the subjects dir;
  #	WE MUST REMEMBER TO RESET!
  if ($?old_subj_dir) setenv SUBJECTS_DIR $old_subj_dir
  
  goto generatempg_return;


#############################################################

    # Every array variable is either 
    #   (1) length 1 (e.g. use the same value for all images
    #          OR
    #   (2) length $nimages (e.g. use the specified value for each image)
    #
    #           anything different than these will cause a run-time error
  
usage_exit:
  echo " "
  echo "USAGE: $progname -subj <names> -srcfile <srcfiles> [options]"
  echo " "
  echo "Required Arguments:";
  echo "  -subj         subject name           : SUBJECTS_DIR must be set correctly"
  echo "  -srcfile      stc, mgh               : full file name (minus hemi and ext)"
  echo " "
  echo "Optional Arguments:";
  #echo "  -mont         pre-named montage          : [$mont]"
  echo "  [Program Flow args:]";
  echo "  -force        don't use any cached images";
  echo "  -nocache      don't save any intermediate images";
#  echo "  -rmintermed   remove intermediate files (wt, w, rgb)"
  echo "  -remakeonly   flag to remake movie with cached images"
  echo " ";
  echo "  [I/O args:]"
  echo "  -indir        input dir              : [$indirs]"
  echo "  -outdir       output dir             : [$outdir]"
  echo "  -outstem      output stem to prepend : [$outstem]"
  echo "  -sdir         SUBJECTS_DIR (required "
  echo "                when ! -e SUBJECTS_DIR : [$sdirs]"
  echo " "
  echo "  [Data args:]"
  echo "  -hemi         hemisphere             : [$hemis]"
  echo "  -polar        look for _r and _i stc "
  echo "                files, render polar    : [$polarflag]"
  echo "  -tfirst       first time point       : [$tfirsts]"
  echo "  -tlast        last time point        : [$tlasts]"
  echo "  -tmax         max time point         : [$tmaxs]"
  echo " "
  echo "  [Image args:]";
  echo "  -surf         surface                : [$surfs]"
  echo "  -fthresh      fthresh                : [$fthreshs]"
  echo "  -fslope       fslope                 : [$fslopes]"
  echo "  -fmid         fmid                   : [$fmids]"
  echo "  -sparsesmooth sparse smooth iters    : [$sparsesmooths]"
  echo "  -postsmooth   regular smooth iters   : [$postsmooths]"
  echo "  -offset       offset                 : [$offsets]"
  echo "  -cvfact       cvfact                 : [$cvfacts]"
  echo "  -scale        zoom factor            : [$scales]"
  echo "  -view         lat,med,ven,pos,dor    : [$views]"
  echo "  -rotx         x rotation             : [$rotxs]"
  echo "  -roty         y rotation             : [$rotys]"
  echo "  -rotz         z rotation             : [$rotzs]"
  echo "  -label        text label             : [$labels]";
#  echo "  -patch        patch        : [$patchs]"
#  echo "  -flat         use flat patch"
#  echo "  -scalebar"
#  echo "  -colscalebar"
  echo " "
  echo "  Mpg args:"
  echo "  -nrows        # rows in montage      : [$nrows]"
  echo "  -ncols        # cols in montage      : [$ncols]"
  echo "  -position     grid cell #s (l/r t/b) : [$positions]"
  echo "  -geometry     image max width/height : [$geometry]"
  echo "  -bgcolor      background fill color  : [$bgcolor]";
  echo " "

  exit 1;

#############################################################

mainprog:
  goto initvars;
  initvars_return:

  goto checkdeps;
  checkdeps_return:
  
  ## parse and check params ##
  goto parse_args;
  parse_args_return:

  goto check_params;
  check_params_return:

  goto prepare_vars;
  prepare_vars_return:

  goto generatempg;
  generatempg_return:

  goto end;
  
  
#############################################################

end:
  exit 0
